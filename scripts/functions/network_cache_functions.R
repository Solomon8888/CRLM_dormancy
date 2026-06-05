# 网络资源与远程API结果缓存函数
#
# 设计目标：
# - 缓存统一放在data/reference/下；
# - 默认7天有效，过期后自动重新获取；
# - 缓存对象带metadata校验，避免参数或输入基因改变后误用旧缓存。

is_reference_cache_fresh <- function(cache_file, max_age_days = 7) {
  if (!file.exists(cache_file)) {
    return(FALSE)
  }

  cache_age_days <- as.numeric(
    difftime(Sys.time(), file.info(cache_file)$mtime, units = "days")
  )
  !is.na(cache_age_days) && cache_age_days <= max_age_days
}

make_gene_cache_metadata <- function(method, input, genes, extra = list()) {
  list(
    method = method,
    input_type = input$input_type,
    input_name = input$input_name,
    genes_key = paste(sort(unique(as.character(genes))), collapse = "\t"),
    extra = extra
  )
}

make_resource_cache_metadata <- function(resource_name, species, extra = list()) {
  list(
    resource = resource_name,
    species = tolower(as.character(species)),
    extra = extra
  )
}

read_reference_cache <- function(
    cache_file,
    expected_metadata,
    max_age_days = 7,
    use_cache = TRUE) {
  if (!use_cache || !is_reference_cache_fresh(cache_file, max_age_days = max_age_days)) {
    return(list(found = FALSE, result = NULL))
  }

  cache_object <- tryCatch(
    readRDS(cache_file),
    error = function(error) NULL
  )
  if (
    is.null(cache_object) ||
      !is.list(cache_object) ||
      !all(c("metadata", "result") %in% names(cache_object))
  ) {
    return(list(found = FALSE, result = NULL))
  }

  if (!identical(cache_object$metadata, expected_metadata)) {
    return(list(found = FALSE, result = NULL))
  }

  list(found = TRUE, result = cache_object$result)
}

write_reference_cache <- function(cache_file, metadata, result) {
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    object = list(metadata = metadata, result = result),
    file = cache_file
  )
  invisible(cache_file)
}

# TF脚本历史函数名兼容；后续新脚本可以直接使用上面的通用命名。
is_tf_cache_fresh <- is_reference_cache_fresh
make_tf_gene_cache_metadata <- make_gene_cache_metadata
make_tf_resource_cache_metadata <- make_resource_cache_metadata
read_tf_cache <- read_reference_cache
write_tf_cache <- write_reference_cache

configure_omnipathr_runtime <- function(
    cache_dir,
    log_dir,
    console_loglevel = "WARN",
    loglevel = "INFO") {
  # OmnipathR默认会在当前工作目录生成omnipathr-log。
  # 这里显式指定缓存和日志目录，确保项目根目录只保留代码、数据和正式结果。
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  options(
    omnipathr.cachedir = cache_dir,
    omnipathr.logdir = log_dir,
    omnipathr.console_loglevel = console_loglevel,
    omnipathr.loglevel = loglevel
  )

  if (requireNamespace("OmnipathR", quietly = TRUE)) {
    init_log <- tryCatch(
      get("omnipath_init_log", envir = asNamespace("OmnipathR")),
      error = function(error) NULL
    )
    if (is.function(init_log)) {
      try(init_log(), silent = TRUE)
    }
  }

  invisible(list(cache_dir = cache_dir, log_dir = log_dir))
}
