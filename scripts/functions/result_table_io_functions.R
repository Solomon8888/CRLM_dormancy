# 结果表格统一保存函数
#
# 约定：
# - 分析脚本仍然传入形如 output_dir/result.csv 的目标路径；
# - 本函数实际保存到 output_dir/csv/result.csv、output_dir/md/result.md、
#   output_dir/tex/result.tex；
# - 这样所有表格类结果都按文件格式分目录保存，便于Beamer、README和后续脚本调用。

escape_latex_cell <- function(x) {
  if (is.na(x) || x == "") {
    return("--")
  }

  chars <- strsplit(as.character(x), "", fixed = TRUE)[[1]]
  replacements <- c(
    "\\" = "\\textbackslash{}",
    "&" = "\\&",
    "%" = "\\%",
    "$" = "\\$",
    "#" = "\\#",
    "_" = "\\_",
    "{" = "\\{",
    "}" = "\\}",
    "~" = "\\textasciitilde{}",
    "^" = "\\textasciicircum{}"
  )

  chars <- ifelse(chars %in% names(replacements), replacements[chars], chars)
  paste(chars, collapse = "")
}

escape_latex_vector <- function(x) {
  vapply(as.character(x), escape_latex_cell, character(1), USE.NAMES = FALSE)
}

escape_markdown_cell <- function(x) {
  if (is.na(x) || x == "") {
    return("--")
  }

  x <- as.character(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\\|", "\\\\|", x)
  x <- gsub("\r\n|\n|\r", "<br>", x)
  x
}

escape_markdown_vector <- function(x) {
  vapply(as.character(x), escape_markdown_cell, character(1), USE.NAMES = FALSE)
}

get_report_table_paths <- function(csv_file) {
  # 接受旧式 output_dir/name.csv，也接受已经位于csv/md/tex目录中的路径。
  csv_file <- as.character(csv_file)
  table_root <- dirname(csv_file)
  table_format_dir <- basename(table_root)

  if (table_format_dir %in% c("csv", "md", "tex")) {
    table_root <- dirname(table_root)
  }

  file_stem <- tools::file_path_sans_ext(basename(csv_file))
  list(
    root = table_root,
    csv = file.path(table_root, "csv", paste0(file_stem, ".csv")),
    md = file.path(table_root, "md", paste0(file_stem, ".md")),
    tex = file.path(table_root, "tex", paste0(file_stem, ".tex"))
  )
}

resolve_report_csv_file <- function(csv_file) {
  # 后续脚本读取表格时统一调用；优先使用新目录结构，兼容旧版平铺CSV。
  paths <- get_report_table_paths(csv_file)
  if (file.exists(paths$csv)) {
    return(paths$csv)
  }
  if (file.exists(csv_file)) {
    return(csv_file)
  }
  paths$csv
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
  # 同一结果若同时存在旧版平铺CSV和新版csv/子目录CSV，优先使用新版文件。
  if (length(files) == 0L) {
    return(files)
  }

  groups <- vapply(files, group_function, character(1))
  is_new_layout <- basename(dirname(files)) == "csv"
  order_index <- order(groups, !is_new_layout, files)
  files <- files[order_index]
  groups <- groups[order_index]
  files[!duplicated(groups)]
}

write_latex_table_preview <- function(csv_file, tex_file = NULL, n_rows = 21) {
  if (is.null(tex_file)) {
    tex_file <- get_report_table_paths(csv_file)$tex
  }

  dat <- read.csv(
    csv_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )

  if (nrow(dat) > n_rows) {
    dat <- dat[seq_len(n_rows), , drop = FALSE]
  }

  dir.create(dirname(tex_file), recursive = TRUE, showWarnings = FALSE)
  col_spec <- paste(rep("l", ncol(dat)), collapse = "")
  header <- paste(escape_latex_vector(colnames(dat)), collapse = " & ")
  rows <- if (nrow(dat) > 0) {
    apply(dat, 1, function(row) {
      paste(escape_latex_vector(row), collapse = " & ")
    })
  } else {
    character(0)
  }

  table_lines <- c(
    "\\begingroup",
    "\\tiny",
    "\\setlength{\\tabcolsep}{2pt}",
    "\\renewcommand{\\arraystretch}{1.18}",
    sprintf("\\begin{tabular}{@{}%s@{}}", col_spec),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    if (length(rows) > 0) paste0(rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\endgroup"
  )

  writeLines(table_lines, tex_file, useBytes = TRUE)
  tex_file
}

write_markdown_table_preview <- function(csv_file, md_file = NULL, n_rows = 21) {
  if (is.null(md_file)) {
    md_file <- get_report_table_paths(csv_file)$md
  }

  dat <- read.csv(
    csv_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )

  if (nrow(dat) > n_rows) {
    dat <- dat[seq_len(n_rows), , drop = FALSE]
  }

  dir.create(dirname(md_file), recursive = TRUE, showWarnings = FALSE)
  header <- paste(escape_markdown_vector(colnames(dat)), collapse = " | ")
  separator <- paste(rep("---", ncol(dat)), collapse = " | ")
  rows <- if (nrow(dat) > 0) {
    apply(dat, 1, function(row) {
      paste(escape_markdown_vector(row), collapse = " | ")
    })
  } else {
    character(0)
  }

  table_lines <- c(
    paste0("| ", header, " |"),
    paste0("| ", separator, " |"),
    paste0("| ", rows, " |")
  )

  writeLines(table_lines, md_file, useBytes = TRUE)
  md_file
}

write_csv_with_report_previews <- function(dat, csv_file, n_rows = 21, ...) {
  paths <- get_report_table_paths(csv_file)
  dir.create(dirname(paths$csv), recursive = TRUE, showWarnings = FALSE)

  write.csv(dat, paths$csv, row.names = FALSE, ...)
  write_latex_table_preview(paths$csv, tex_file = paths$tex, n_rows = n_rows)
  write_markdown_table_preview(paths$csv, md_file = paths$md, n_rows = n_rows)

  invisible(paths$csv)
}

write_report_table_formats <- function(dat, output_dir, file_stem, n_rows = 21, ...) {
  write_csv_with_report_previews(
    dat = dat,
    csv_file = file.path(output_dir, paste0(file_stem, ".csv")),
    n_rows = n_rows,
    ...
  )
}
