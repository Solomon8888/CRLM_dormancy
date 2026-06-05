# 结果表格统一保存与读取函数
#
# 当前约定：
# - 分析脚本只保存平铺 CSV，例如 output_dir/result.csv；
# - 不再由 R 脚本额外生成 md/tex 预览文件；
# - Beamer/README 需要的 tex 或 markdown 预览，由对应构建脚本基于 CSV 统一生成；
# - 读取函数兼容历史 output_dir/csv/result.csv 结构，避免旧结果或未重跑结果无法读取。

normalize_report_csv_path <- function(csv_file) {
  csv_file <- as.character(csv_file)
  parent_dir <- basename(dirname(csv_file))

  if (parent_dir %in% c("csv", "md", "tex")) {
    return(file.path(
      dirname(dirname(csv_file)),
      paste0(tools::file_path_sans_ext(basename(csv_file)), ".csv")
    ))
  }

  csv_file
}

legacy_report_csv_path <- function(csv_file) {
  csv_file <- normalize_report_csv_path(csv_file)
  file.path(
    dirname(csv_file),
    "csv",
    paste0(tools::file_path_sans_ext(basename(csv_file)), ".csv")
  )
}

resolve_report_csv_file <- function(csv_file) {
  # 优先读取当前平铺CSV；若尚未重跑脚本，则兼容旧版csv/子目录CSV。
  flat_file <- normalize_report_csv_path(csv_file)
  legacy_file <- legacy_report_csv_path(csv_file)

  if (file.exists(flat_file)) {
    return(flat_file)
  }
  if (file.exists(legacy_file)) {
    return(legacy_file)
  }

  flat_file
}

read_report_csv <- function(csv_file, ...) {
  read.csv(
    resolve_report_csv_file(csv_file),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    ...
  )
}

prefer_report_csv_files <- function(files, group_function) {
  # 同一结果若同时存在平铺CSV和历史csv/子目录CSV，优先使用平铺CSV。
  if (length(files) == 0L) {
    return(files)
  }

  groups <- vapply(files, group_function, character(1))
  is_legacy_layout <- basename(dirname(files)) == "csv"
  order_index <- order(groups, is_legacy_layout, files)
  files <- files[order_index]
  groups <- groups[order_index]

  files[!duplicated(groups)]
}

format_report_csv_values <- function(dat) {
  # 仅用于写出CSV前的展示格式整理，不改变分析脚本中的原始运算对象。
  # 所有数值列都按普通十进制字符写出，避免p值、q值、Entrez等字段出现科学计数法。
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)

  for (column_name in names(dat)) {
    if (is.numeric(dat[[column_name]])) {
      dat[[column_name]] <- ifelse(
        is.na(dat[[column_name]]),
        NA,
        format(
          dat[[column_name]],
          scientific = FALSE,
          trim = TRUE,
          digits = 16
        )
      )
    }
  }

  dat
}

write_csv_with_report_previews <- function(dat, csv_file, n_rows = 21, ...) {
  # 保留历史函数名，避免逐个分析脚本改调用；现在仅写CSV。
  output_file <- normalize_report_csv_path(csv_file)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(format_report_csv_values(dat), output_file, row.names = FALSE, ...)
  invisible(output_file)
}

write_report_table_formats <- function(dat, output_dir, file_stem, n_rows = 21, ...) {
  write_csv_with_report_previews(
    dat = dat,
    csv_file = file.path(output_dir, paste0(file_stem, ".csv")),
    n_rows = n_rows,
    ...
  )
}
