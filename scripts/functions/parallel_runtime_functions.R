# 批量任务并行运行公共函数
#
# 本文件只保存与“批量任务加速运行”有关的通用函数。
# 分析脚本只需要提供任务编号、任务函数和总任务数，即可获得稳定的多进程运行、
# 终端进度条、qs2线程配置、BiocParallel配置和运行时间统计。


start_runtime_timer <- function() {
  # 返回当前时间，供脚本末尾计算总运行时间。
  Sys.time()
}

format_runtime_seconds <- function(runtime_seconds) {
  # 将秒数格式化为终端容易阅读的“分钟+秒”形式。
  runtime_seconds <- as.numeric(runtime_seconds)
  sprintf("%.2f min (%.1f sec)", runtime_seconds / 60, runtime_seconds)
}

print_runtime_summary <- function(start_time, label = "Total runtime") {
  # 按统一格式打印从start_time到当前时刻的总运行时间。
  runtime_seconds <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  )
  cat(label, ": ", format_runtime_seconds(runtime_seconds), "\n", sep = "")
  invisible(runtime_seconds)
}

get_available_worker_count <- function() {
  # logical=TRUE可充分调用Apple Silicon的性能核心与效率核心。
  max(1L, parallel::detectCores(logical = TRUE))
}

make_parallel_execution_strategy <- function(total_tasks, max_workers = NULL) {
  # 默认策略：
  # 1. 任务数充足时优先做外层任务级并行，让CPU持续满负荷；
  # 2. 任务数较少时，单个任务自动获得更多内部线程；
  # 3. 外层并行开启时，任务内部的二级并行保持为1，避免进程过度抢占。
  total_tasks <- max(as.integer(total_tasks), 1L)
  if (is.null(max_workers)) {
    max_workers <- get_available_worker_count()
  }
  max_workers <- max(as.integer(max_workers), 1L)

  task_workers <- if (.Platform$OS.type != "windows") {
    min(max_workers, total_tasks)
  } else {
    1L
  }

  inner_workers <- max(1L, floor(max_workers / task_workers))
  nested_workers <- if (task_workers > 1L) {
    1L
  } else {
    max_workers
  }

  list(
    total_tasks = total_tasks,
    max_workers = max_workers,
    task_workers = task_workers,
    inner_workers = inner_workers,
    qs2_threads_per_task = inner_workers,
    nested_workers = nested_workers
  )
}

configure_parallel_runtime <- function(
    task_workers,
    inner_workers = 1L,
    qs2_threads = inner_workers) {
  # 统一配置R parallel、qs2和BiocParallel。
  # 当外层任务并行时，每个子任务继承这里的qs2线程数，避免每个任务都抢占全部核心。
  task_workers <- max(as.integer(task_workers), 1L)
  inner_workers <- max(as.integer(inner_workers), 1L)
  qs2_threads <- max(as.integer(qs2_threads), 1L)

  options(mc.cores = task_workers)

  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qopt("nthreads", qs2_threads)
  }

  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(
      BiocParallel::MulticoreParam(workers = inner_workers),
      default = TRUE
    )
  }

  invisible(list(
    task_workers = task_workers,
    inner_workers = inner_workers,
    qs2_threads = qs2_threads
  ))
}

.parallel_progress_state <- new.env(parent = emptyenv())

is_progress_terminal <- function() {
  # Rscript在VSCode/终端中运行时通常为TTY，可用单行覆盖式进度条；
  # 若输出被日志系统捕获，则降频打印，避免形成大量重复进度行。
  tryCatch(isatty(stdout()), error = function(error) FALSE)
}

reset_progress_state <- function() {
  .parallel_progress_state$last_plain_percent <- -Inf
  invisible(TRUE)
}

print_parallel_execution_strategy <- function(
    strategy,
    inner_label = "Inner workers per task",
    nested_label = "Nested workers") {
  # 统一打印批量任务并行策略，便于运行前快速确认CPU使用方式。
  print_parallel_metric <- function(label, value) {
    cat(sprintf("%-28s %s\n", paste0(label, ":"), value))
  }

  cat("\nParallel execution strategy:\n")
  print_parallel_metric("Total tasks", strategy$total_tasks)
  print_parallel_metric("Available workers", strategy$max_workers)
  print_parallel_metric("Task-level workers", strategy$task_workers)
  print_parallel_metric(inner_label, strategy$inner_workers)
  print_parallel_metric("qs2 nthreads per task", strategy$qs2_threads_per_task)
  print_parallel_metric(nested_label, strategy$nested_workers)
  invisible(strategy)
}

format_elapsed_time <- function(seconds) {
  # 把秒数压缩成HH:MM:SS或MM:SS，方便进度条实时展示。
  seconds <- max(0, as.integer(round(as.numeric(seconds))))
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  seconds <- seconds %% 60

  if (hours > 0) {
    return(sprintf("%02d:%02d:%02d", hours, minutes, seconds))
  }

  sprintf("%02d:%02d", minutes, seconds)
}

make_progress_line <- function(
    completed,
    total,
    start_time,
    width = 30L,
    label = "Progress") {
  # 构造单行进度条文本。只返回字符串，不直接写终端，便于复用和测试。
  completed <- min(max(as.integer(completed), 0L), as.integer(total))
  total <- max(as.integer(total), 1L)
  width <- max(as.integer(width), 10L)

  progress_fraction <- completed / total
  filled_width <- floor(progress_fraction * width)
  empty_width <- width - filled_width

  bar <- paste0(
    strrep("\u2588", filled_width),
    strrep("\u2591", empty_width)
  )

  elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  eta_text <- "--:--"
  if (completed > 0L && completed < total) {
    eta_seconds <- elapsed_seconds / completed * (total - completed)
    eta_text <- format_elapsed_time(eta_seconds)
  }
  if (completed >= total) {
    eta_text <- "00:00"
  }

  sprintf(
    "%s %6.1f%% |%s| %d/%d | elapsed %s | ETA %s",
    label,
    progress_fraction * 100,
    bar,
    completed,
    total,
    format_elapsed_time(elapsed_seconds),
    eta_text
  )
}

render_progress_line <- function(line, force = FALSE) {
  # 使用回车和ANSI清行在同一行刷新，避免终端刷屏。
  # 如果终端不支持ANSI，回车仍可让大多数终端覆盖当前行。
  max_width <- max(getOption("width", 120L) - 1L, 60L)
  if (nchar(line) > max_width) {
    line <- substr(line, 1L, max_width)
  }

  if (!is_progress_terminal()) {
    percent_match <- regexpr("[0-9]+\\.[0-9]+%", line)
    progress_percent <- NA_real_
    if (percent_match > 0L) {
      progress_percent <- as.numeric(
        sub("%", "", regmatches(line, percent_match))
      )
    }

    last_percent <- .parallel_progress_state$last_plain_percent
    if (is.null(last_percent)) {
      last_percent <- -Inf
    }

    should_print <- force ||
      is.na(progress_percent) ||
      progress_percent <= 0 ||
      progress_percent >= 100 ||
      progress_percent - last_percent >= 10

    if (should_print) {
      cat(line, "\n", sep = "")
      .parallel_progress_state$last_plain_percent <- progress_percent
      flush.console()
    }
    return(invisible(line))
  }

  cat("\r\033[2K", line, sep = "")
  flush.console()
  invisible(line)
}

finish_progress_line <- function(line = NULL) {
  if (!is.null(line)) {
    render_progress_line(line, force = TRUE)
  }
  if (is_progress_terminal()) {
    cat("\n")
  }
  flush.console()
}

execute_parallel_task <- function(task_id, task_function, suppress_output = TRUE) {
  # 批量子任务默认不直接向终端写普通输出，避免打断主进程单行进度条。
  # 错误仍通过try-error对象返回，主进程随后统一汇总失败任务。
  if (!suppress_output) {
    return(try(task_function(task_id), silent = TRUE))
  }

  task_result <- NULL
  invisible(capture.output({
    task_result <- try(
      suppressMessages(task_function(task_id)),
      silent = TRUE
    )
  }))
  task_result
}

setup_parallel_strategy <- function(
    total_tasks,
    max_workers = NULL,
    inner_label = "Inner workers per task",
    nested_label = "Nested workers",
    print_strategy = TRUE) {
  # 一站式生成并配置批量任务并行策略。
  # 返回值中的task_workers用于外层run_parallel_tasks_with_progress；
  # inner_workers和qs2_threads_per_task已在本函数内配置到R运行环境。
  strategy <- make_parallel_execution_strategy(
    total_tasks = total_tasks,
    max_workers = max_workers
  )

  configure_parallel_runtime(
    task_workers = strategy$task_workers,
    inner_workers = strategy$inner_workers,
    qs2_threads = strategy$qs2_threads_per_task
  )

  if (print_strategy) {
    print_parallel_execution_strategy(
      strategy = strategy,
      inner_label = inner_label,
      nested_label = nested_label
    )
  }

  strategy
}

run_parallel_tasks_with_progress <- function(
    task_ids,
    task_function,
    workers,
    suppress_task_output = TRUE,
    progress_label = "Progress") {
  # 主进程维护进度条，子进程负责实际批量任务。
  # 与parallel::mclapply相比，这里可以实时回收已完成任务并刷新终端进度。
  total_task_count <- length(task_ids)
  if (total_task_count == 0L) {
    return(list())
  }

  progress_start_time <- Sys.time()
  reset_progress_state()
  render_progress_line(
    make_progress_line(0L, total_task_count, progress_start_time, label = progress_label),
    force = TRUE
  )
  progress_finished <- FALSE
  on.exit({
    if (!progress_finished) {
      finish_progress_line()
    }
  }, add = TRUE)

  results <- vector("list", total_task_count)
  names(results) <- as.character(task_ids)

  if (workers <= 1L || .Platform$OS.type == "windows") {
    for (task_position in seq_along(task_ids)) {
      task_id <- task_ids[task_position]
      results[[as.character(task_id)]] <- execute_parallel_task(
        task_id = task_id,
        task_function = task_function,
        suppress_output = suppress_task_output
      )
      render_progress_line(
        make_progress_line(
          task_position,
          total_task_count,
          progress_start_time,
          label = progress_label
        ),
        force = task_position >= total_task_count
      )
    }

    finish_progress_line()
    progress_finished <- TRUE
    return(results)
  }

  workers <- max(1L, min(as.integer(workers), total_task_count))
  next_task_position <- 1L
  completed_tasks <- 0L
  active_jobs <- list()
  active_task_by_pid <- list()

  launch_next_task <- function() {
    task_id <- task_ids[next_task_position]
    job <- parallel::mcparallel(
      execute_parallel_task(
        task_id = task_id,
        task_function = task_function,
        suppress_output = suppress_task_output
      ),
      silent = TRUE
    )
    pid <- as.character(job$pid)

    active_jobs[[pid]] <<- job
    active_task_by_pid[[pid]] <<- task_id
    next_task_position <<- next_task_position + 1L
  }

  while (
    next_task_position <= total_task_count &&
      length(active_jobs) < workers
  ) {
    launch_next_task()
  }

  while (completed_tasks < total_task_count) {
    ready_results <- parallel::mccollect(
      active_jobs,
      wait = FALSE,
      timeout = 0.5
    )

    if (is.null(ready_results) || length(ready_results) == 0) {
      Sys.sleep(0.2)
      next
    }

    for (pid in names(ready_results)) {
      task_id <- active_task_by_pid[[pid]]
      results[[as.character(task_id)]] <- ready_results[[pid]]
      active_jobs[[pid]] <- NULL
      active_task_by_pid[[pid]] <- NULL
      completed_tasks <- completed_tasks + 1L
      render_progress_line(
        make_progress_line(
          completed_tasks,
          total_task_count,
          progress_start_time,
          label = progress_label
        ),
        force = completed_tasks >= total_task_count
      )

      while (
        next_task_position <= total_task_count &&
          length(active_jobs) < workers
      ) {
        launch_next_task()
      }
    }
  }

  finish_progress_line()
  progress_finished <- TRUE
  results
}

run_indexed_tasks_with_progress <- function(
    total_tasks,
    task_function,
    workers,
    suppress_task_output = TRUE,
    progress_label = "Progress") {
  # 适合大多数批量脚本：任务天然按1:n编号。
  run_parallel_tasks_with_progress(
    task_ids = seq_len(total_tasks),
    task_function = task_function,
    workers = workers,
    suppress_task_output = suppress_task_output,
    progress_label = progress_label
  )
}

stop_on_parallel_errors <- function(results, task_ids = names(results), label = "parallel tasks") {
  # 统一检查并行任务返回值中的try-error，便于各脚本保持相同的失败提示。
  task_errors <- vapply(results, function(x) {
    is.null(x) || inherits(x, "try-error")
  }, logical(1))
  if (!any(task_errors)) {
    return(invisible(FALSE))
  }

  failed_tasks <- task_ids[task_errors]
  stop(
    "Some ",
    label,
    " failed: ",
    paste(failed_tasks, collapse = ", ")
  )
}
