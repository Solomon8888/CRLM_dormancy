# Common table helpers for report-ready result previews.

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

write_latex_table_preview <- function(csv_file, tex_file = NULL, n_rows = 21) {
  if (is.null(tex_file)) {
    tex_file <- sub("[.]csv$", ".tex", csv_file)
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

write_csv_with_latex_preview <- function(dat, csv_file, n_rows = 21, ...) {
  write.csv(dat, csv_file, row.names = FALSE, ...)
  write_latex_table_preview(csv_file, n_rows = n_rows)
}
