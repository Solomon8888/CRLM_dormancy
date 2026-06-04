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

print_parallel_execution_strategy <- function(
    strategy,
    inner_label = "Inner workers per task",
    nested_label = "Nested workers") {
  # 统一打印批量任务并行策略，便于运行前快速确认CPU使用方式。
  cat("\nParallel execution strategy:\n")
  cat("Total tasks:              ", strategy$total_tasks, "\n", sep = "")
  cat("Available workers:        ", strategy$max_workers, "\n", sep = "")
  cat("Task-level workers:       ", strategy$task_workers, "\n", sep = "")
  cat(inner_label, ":      ", strategy$inner_workers, "\n", sep = "")
  cat("qs2 nthreads per task:    ", strategy$qs2_threads_per_task, "\n", sep = "")
  cat(nested_label, ":       ", strategy$nested_workers, "\n", sep = "")
  invisible(strategy)
}

run_parallel_tasks_with_progress <- function(task_ids, task_function, workers) {
  # 主进程维护进度条，子进程负责实际批量任务。
  # 与parallel::mclapply相比，这里可以实时回收已完成任务并刷新终端进度。
  total_task_count <- length(task_ids)
  if (total_task_count == 0L) {
    return(list())
  }

  progress_bar <- utils::txtProgressBar(
    min = 0,
    max = total_task_count,
    style = 3
  )
  on.exit(close(progress_bar), add = TRUE)

  results <- vector("list", total_task_count)
  names(results) <- as.character(task_ids)

  if (workers <= 1L || .Platform$OS.type == "windows") {
    for (task_position in seq_along(task_ids)) {
      task_id <- task_ids[task_position]
      results[[as.character(task_id)]] <- try(task_function(task_id), silent = TRUE)
      utils::setTxtProgressBar(progress_bar, task_position)
    }

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
      try(task_function(task_id), silent = TRUE),
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
      utils::setTxtProgressBar(progress_bar, completed_tasks)

      while (
        next_task_position <= total_task_count &&
          length(active_jobs) < workers
      ) {
        launch_next_task()
      }
    }
  }

  results
}
