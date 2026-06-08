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

normalize_parallel_backend <- function(backend) {
  backend <- tolower(trimws(as.character(backend)[1]))
  if (!nzchar(backend) || is.na(backend)) {
    backend <- "auto"
  }
  if (!backend %in% c("auto", "fork", "psock", "serial")) {
    backend <- "auto"
  }
  backend
}

get_parallel_backend_preference <- function() {
  normalize_parallel_backend(
    getOption(
      "parallel_runtime_backend",
      Sys.getenv("PARALLEL_RUNTIME_BACKEND", unset = "auto")
    )
  )
}

is_parallel_fork_allowed <- function() {
  # 交互式source()时，macOS/RStudio/部分R终端中的fork并行容易触发进程级崩溃。
  # quickanalysis默认在交互式会话禁用fork，并自动切换为PSOCK多进程。
  if (isTRUE(getOption("parallel_runtime_disable_fork", FALSE))) {
    return(FALSE)
  }

  if (.Platform$OS.type == "windows") {
    return(FALSE)
  }

  if (interactive() &&
      !isTRUE(getOption("parallel_runtime_allow_interactive_fork", FALSE))) {
    return(FALSE)
  }

  TRUE
}

choose_parallel_backend <- function(max_workers) {
  max_workers <- max(as.integer(max_workers), 1L)
  preference <- get_parallel_backend_preference()
  fork_allowed <- is_parallel_fork_allowed()

  if (max_workers <= 1L || identical(preference, "serial")) {
    return("serial")
  }

  if (identical(preference, "fork")) {
    if (fork_allowed) {
      return("fork")
    }
    return("psock")
  }

  if (identical(preference, "psock")) {
    return("psock")
  }

  if (fork_allowed) {
    return("fork")
  }

  "psock"
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

  backend <- choose_parallel_backend(max_workers)
  parallel_capable <- backend %in% c("fork", "psock")

  task_workers <- if (parallel_capable) {
    min(max_workers, total_tasks)
  } else {
    1L
  }

  inner_workers <- if (parallel_capable) {
    max(1L, floor(max_workers / task_workers))
  } else {
    1L
  }
  nested_workers <- if (parallel_capable && task_workers <= 1L) {
    max_workers
  } else {
    1L
  }

  list(
    backend = backend,
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
    qs2_threads = inner_workers,
    backend = getOption("parallel_runtime_effective_backend", choose_parallel_backend(task_workers))) {
  # 统一配置R parallel、qs2和BiocParallel。
  # 当外层任务并行时，每个子任务继承这里的qs2线程数，避免每个任务都抢占全部核心。
  task_workers <- max(as.integer(task_workers), 1L)
  inner_workers <- max(as.integer(inner_workers), 1L)
  qs2_threads <- max(as.integer(qs2_threads), 1L)
  backend <- normalize_parallel_backend(backend)
  if (backend == "auto") {
    backend <- choose_parallel_backend(task_workers)
  }

  options(mc.cores = task_workers)
  options(parallel_runtime_effective_backend = backend)

  configure_parallel_math_threads(qs2_threads)

  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qopt("nthreads", qs2_threads)
  }

  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    if (identical(backend, "fork") && is_parallel_fork_allowed() && inner_workers > 1L) {
      BiocParallel::register(
        BiocParallel::MulticoreParam(workers = inner_workers),
        default = TRUE
      )
    } else if (identical(backend, "psock") && inner_workers > 1L) {
      BiocParallel::register(
        BiocParallel::SnowParam(workers = inner_workers, type = "SOCK"),
        default = TRUE
      )
    } else {
      BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
    }
  }

  invisible(list(
    backend = backend,
    task_workers = task_workers,
    inner_workers = inner_workers,
    qs2_threads = qs2_threads
  ))
}

configure_parallel_math_threads <- function(thread_count) {
  # 控制BLAS/OpenMP/Accelerate线程数。外层任务并行时每个worker用少量线程，
  # 单任务场景则把线程数让给内部GSEA/矩阵计算，避免过度抢核或闲核。
  thread_count <- max(as.integer(thread_count), 1L)
  Sys.setenv(
    OMP_NUM_THREADS = as.character(thread_count),
    OPENBLAS_NUM_THREADS = as.character(thread_count),
    MKL_NUM_THREADS = as.character(thread_count),
    VECLIB_MAXIMUM_THREADS = as.character(thread_count),
    NUMEXPR_NUM_THREADS = as.character(thread_count)
  )

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(thread_count), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(thread_count), silent = TRUE)
  }

  invisible(thread_count)
}

.parallel_progress_state <- new.env(parent = emptyenv())

is_progress_terminal <- function() {
  # Rscript在VSCode/终端中运行时通常为TTY，可用单行覆盖式进度条；
  # 若输出被日志系统捕获，则降频打印，避免形成大量重复进度行。
  if (isTRUE(getOption("parallel_runtime_force_single_line_progress", FALSE))) {
    return(TRUE)
  }

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
  print_parallel_metric("Backend", strategy$backend)
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

get_parallel_worker_packages <- function() {
  configured <- getOption("parallel_runtime_worker_packages", NULL)
  if (!is.null(configured)) {
    return(unique(as.character(configured)))
  }

  attached_packages <- grep("^package:", search(), value = TRUE)
  attached_packages <- sub("^package:", "", attached_packages)
  base_packages <- c(
    "base", "compiler", "datasets", "graphics", "grDevices",
    "grid", "methods", "parallel", "splines", "stats", "stats4",
    "tools", "utils"
  )
  unique(setdiff(attached_packages, base_packages))
}

get_parallel_export_names <- function() {
  configured <- getOption("parallel_runtime_export_names", NULL)
  if (!is.null(configured)) {
    return(intersect(unique(as.character(configured)), ls(envir = .GlobalEnv, all.names = TRUE)))
  }

  global_names <- ls(envir = .GlobalEnv, all.names = TRUE)
  exclude_patterns <- getOption(
    "parallel_runtime_export_exclude_patterns",
    c(
      "^\\.Random\\.seed$",
      "^analysis_results$",
      "^analysis_summaries$",
      "^expression_input_list$",
      "^gsea_result$",
      "^volcano_result$",
      "^summary_records$",
      "^summary_list$",
      "^raw_task_results$",
      "^normalized_task_results$"
    )
  )
  if (length(exclude_patterns) > 0L) {
    should_exclude <- vapply(global_names, function(name) {
      any(grepl(paste(exclude_patterns, collapse = "|"), name))
    }, logical(1))
    global_names <- global_names[!should_exclude]
  }

  max_export_bytes <- as.numeric(
    getOption("parallel_runtime_max_export_object_size", 256 * 1024^2)
  )
  keep <- vapply(global_names, function(name) {
    object <- try(get(name, envir = .GlobalEnv), silent = TRUE)
    if (inherits(object, "try-error")) {
      return(FALSE)
    }
    is.function(object) || as.numeric(utils::object.size(object)) <= max_export_bytes
  }, logical(1))

  global_names[keep]
}

bootstrap_psock_worker <- function(packages, qs2_threads, inner_workers) {
  options(stringsAsFactors = FALSE)
  configure_parallel_math_threads(qs2_threads)

  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qopt("nthreads", qs2_threads)
  }

  if (length(packages) > 0L) {
    invisible(lapply(
      packages,
      function(package) {
        suppressPackageStartupMessages(
          require(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        )
      }
    ))
  }

  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    if (inner_workers > 1L) {
      BiocParallel::register(
        BiocParallel::SnowParam(workers = inner_workers, type = "SOCK"),
        default = TRUE
      )
    } else {
      BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
    }
  }

  TRUE
}

make_psock_cluster <- function(workers, qs2_threads, inner_workers) {
  workers <- max(as.integer(workers), 1L)
  worker_outfile <- getOption(
    "parallel_runtime_worker_outfile",
    if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
  )
  cluster <- parallel::makeCluster(workers, type = "PSOCK", outfile = worker_outfile)

  packages <- get_parallel_worker_packages()
  parallel::clusterExport(
    cluster,
    varlist = c("configure_parallel_math_threads", "bootstrap_psock_worker"),
    envir = environment()
  )
  parallel::clusterCall(
    cluster,
    bootstrap_psock_worker,
    packages = packages,
    qs2_threads = qs2_threads,
    inner_workers = inner_workers
  )

  export_names <- get_parallel_export_names()
  if (length(export_names) > 0L) {
    parallel::clusterExport(cluster, varlist = export_names, envir = .GlobalEnv)
  }

  cluster
}

run_psock_tasks_with_progress <- function(
    task_ids,
    task_function,
    workers,
    suppress_task_output,
    progress_start_time,
    progress_label) {
  workers <- max(1L, min(as.integer(workers), length(task_ids)))
  qs2_threads <- max(as.integer(getOption("parallel_runtime_worker_qs2_threads", 1L)), 1L)
  inner_workers <- max(as.integer(getOption("parallel_runtime_worker_inner_workers", 1L)), 1L)

  cluster <- make_psock_cluster(
    workers = workers,
    qs2_threads = qs2_threads,
    inner_workers = inner_workers
  )
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  worker_task_function <- function(task_id, task_function, suppress_task_output) {
    local_execute_parallel_task <- function(task_id, task_function, suppress_output = TRUE) {
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

    local_execute_parallel_task(
      task_id = task_id,
      task_function = task_function,
      suppress_output = suppress_task_output
    )
  }

  results <- vector("list", length(task_ids))
  names(results) <- as.character(task_ids)
  active_task_by_node <- rep(NA_integer_, length(cluster))
  next_task_position <- 1L
  completed_tasks <- 0L

  launch_task <- function(node_index, task_position) {
    task_id <- task_ids[[task_position]]
    parallel:::sendCall(
      cluster[[node_index]],
      worker_task_function,
      list(task_id, task_function, suppress_task_output)
    )
    active_task_by_node[[node_index]] <<- task_position
  }

  initial_tasks <- min(workers, length(task_ids))
  for (node_index in seq_len(initial_tasks)) {
    launch_task(node_index, next_task_position)
    next_task_position <- next_task_position + 1L
  }

  while (completed_tasks < length(task_ids)) {
    received <- parallel:::recvOneResult(cluster)
    node_index <- received$node
    task_position <- active_task_by_node[[node_index]]
    task_id <- task_ids[[task_position]]

    results[[as.character(task_id)]] <- received$value
    active_task_by_node[[node_index]] <- NA_integer_
    completed_tasks <- completed_tasks + 1L

    render_progress_line(
      make_progress_line(
        completed_tasks,
        length(task_ids),
        progress_start_time,
        label = progress_label
      ),
      force = completed_tasks >= length(task_ids)
    )

    if (next_task_position <= length(task_ids)) {
      launch_task(node_index, next_task_position)
      next_task_position <- next_task_position + 1L
    }
  }

  results
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
    qs2_threads = strategy$qs2_threads_per_task,
    backend = strategy$backend
  )
  options(
    parallel_runtime_worker_qs2_threads = strategy$qs2_threads_per_task,
    parallel_runtime_worker_inner_workers = strategy$inner_workers
  )

  if (print_strategy) {
    print_strategy <- !isTRUE(getOption("parallel_runtime_quiet_strategy", FALSE))
  }

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

  backend <- normalize_parallel_backend(
    getOption("parallel_runtime_effective_backend", choose_parallel_backend(workers))
  )
  if (backend == "auto") {
    backend <- choose_parallel_backend(workers)
  }
  if (backend == "fork" && !is_parallel_fork_allowed()) {
    backend <- "psock"
  }

  if (workers <= 1L || identical(backend, "serial")) {
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

  if (identical(backend, "psock") || .Platform$OS.type == "windows") {
    results <- run_psock_tasks_with_progress(
      task_ids = task_ids,
      task_function = task_function,
      workers = workers,
      suppress_task_output = suppress_task_output,
      progress_start_time = progress_start_time,
      progress_label = progress_label
    )
    finish_progress_line()
    progress_finished <- TRUE
    return(results)
  }

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
