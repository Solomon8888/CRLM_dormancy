# 绘图脚本公共函数
#
# 本文件只保存多个绘图脚本都会用到的基础函数。
# 具体图片的颜色、字体大小、边框粗细、PDF尺寸等样式参数仍放在各绘图脚本头部，
# 这样既方便统一维护公共逻辑，也避免不同图形之间互相影响。


sanitize_file_name <- function(x, default = "analysis") {
  # 将分析名或分组名转换成适合文件夹/文件名使用的字符串。
  # 只保留字母、数字、点、下划线和短横线，避免不同操作系统下路径出错。
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- default
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == ""] <- default
  x
}

read_scalar_config <- function(script_file, variable_name, default_value) {
  # 静态读取另一个R脚本里的单值配置。
  # 这里不source目标脚本，避免为了读取阈值而意外重复运行差异分析。
  if (!file.exists(script_file)) {
    return(default_value)
  }

  script_lines <- readLines(script_file, warn = FALSE)
  matched_line <- grep(
    paste0("^\\s*", variable_name, "\\s*<-"),
    script_lines,
    value = TRUE
  )

  if (length(matched_line) == 0) {
    return(default_value)
  }

  value_text <- sub("#.*$", "", matched_line[1])
  value_text <- sub(
    paste0("^\\s*", variable_name, "\\s*<-\\s*"),
    "",
    value_text
  )
  value_text <- trimws(value_text)

  if (is.character(default_value)) {
    value_text <- gsub('^"|"$', "", value_text)
    return(value_text)
  }

  as.numeric(value_text)
}

wrap_label <- function(x, width = 45) {
  # 按固定字符数对长标签换行，主要用于热图样本名。
  x <- as.character(x)

  vapply(x, function(label) {
    n <- nchar(label)
    if (n <= width) return(label)

    starts <- seq(1, n, by = width)
    parts <- substring(label, starts, pmin(starts + width - 1, n))
    paste(parts, collapse = "\n")
  }, character(1))
}

wrap_label_by_underscore <- function(x, width = 10) {
  # 优先按下划线拆分标签；单段仍过长时再按固定宽度换行。
  x <- as.character(x)

  vapply(x, function(label) {
    parts <- strsplit(label, "_", fixed = TRUE)[[1]]

    if (length(parts) == 1) {
      return(wrap_label(label, width))
    }

    tokens <- paste0(parts, "_")
    tokens[length(tokens)] <- parts[length(parts)]

    lines <- character(0)
    current_line <- ""

    for (token in tokens) {
      candidate <- paste0(current_line, token)

      if (nchar(candidate) <= width || current_line == "") {
        current_line <- candidate
      } else {
        lines <- c(lines, current_line)
        current_line <- token
      }
    }

    paste(c(lines, current_line), collapse = "\n")
  }, character(1))
}

get_label_line_count <- function(label_text) {
  # 统计换行标签的实际行数，用于动态调整标签框或PDF尺寸。
  length(strsplit(label_text, "\n", fixed = TRUE)[[1]])
}

get_deg_file_info <- function(table_root, deg_dir_name = "DEG") {
  # 查找01号差异分析脚本输出的tables/<analysis_name>/DEG/all_genes.csv。
  all_gene_files <- list.files(
    table_root,
    pattern = "^all_genes[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  all_gene_files <- all_gene_files[
    basename(dirname(all_gene_files)) == deg_dir_name
  ]

  stopifnot(length(all_gene_files) > 0)

  file_info <- data.frame(
    Analysis_Name = basename(dirname(dirname(all_gene_files))),
    All_Genes_File = all_gene_files,
    stringsAsFactors = FALSE
  )

  file_info <- file_info[order(file_info$Analysis_Name), , drop = FALSE]
  rownames(file_info) <- NULL

  if (any(duplicated(file_info$Analysis_Name))) {
    duplicated_names <- unique(file_info$Analysis_Name[
      duplicated(file_info$Analysis_Name)
    ])
    stop(
      "More than one all_genes.csv file was found for: ",
      paste(duplicated_names, collapse = ", ")
    )
  }

  file_info
}

prepare_volcano_data <- function(
    dat,
    analysis_name = NULL,
    p_value_column,
    p_value_cutoff,
    logfc_cutoff,
    ns_label = "Not significant",
    regulation_levels = c(ns_label, "Down", "Up")) {
  # 将差异分析结果整理成火山图使用的数据。
  # 显著性判定保持为：P值小于阈值且abs(logFC)大于阈值。
  stopifnot("logFC" %in% colnames(dat))
  stopifnot(p_value_column %in% colnames(dat))

  dat$logFC <- as.numeric(dat$logFC)
  dat[[p_value_column]] <- as.numeric(dat[[p_value_column]])

  valid_index <- is.finite(dat$logFC) &
    is.finite(dat[[p_value_column]]) &
    !is.na(dat[[p_value_column]])

  dat <- dat[valid_index, , drop = FALSE]
  stopifnot(nrow(dat) > 0)

  positive_p <- dat[[p_value_column]][dat[[p_value_column]] > 0]
  stopifnot(length(positive_p) > 0)

  min_positive_p <- min(positive_p, na.rm = TRUE)
  safe_p <- dat[[p_value_column]]
  safe_p[safe_p <= 0] <- min_positive_p * 0.1

  dat$Neg_Log10_P <- -log10(safe_p)
  if (!is.null(analysis_name)) {
    dat$Analysis_Name <- analysis_name
  }

  dat$Regulation <- ns_label
  dat$Regulation[
    dat$logFC > logfc_cutoff &
      dat[[p_value_column]] < p_value_cutoff
  ] <- "Up"
  dat$Regulation[
    dat$logFC < -logfc_cutoff &
      dat[[p_value_column]] < p_value_cutoff
  ] <- "Down"

  dat$Regulation <- factor(dat$Regulation, levels = regulation_levels)
  dat[order(dat$Regulation), , drop = FALSE]
}

has_custom_label_genes <- function(custom_label_genes) {
  # 自定义标注基因可为character向量，也可为按分析名命名的list。
  if (is.list(custom_label_genes)) {
    return(length(unlist(custom_label_genes, use.names = FALSE)) > 0)
  }

  length(custom_label_genes) > 0
}

get_custom_genes_for_analysis <- function(analysis_name, custom_label_genes) {
  # 从自定义配置中提取当前分析需要标注的基因。
  # list中的all表示所有分析都标注。
  if (!has_custom_label_genes(custom_label_genes)) {
    return(character(0))
  }

  if (is.list(custom_label_genes)) {
    global_genes <- if ("all" %in% names(custom_label_genes)) {
      custom_label_genes[["all"]]
    } else {
      character(0)
    }
    analysis_genes <- if (analysis_name %in% names(custom_label_genes)) {
      custom_label_genes[[analysis_name]]
    } else {
      character(0)
    }

    return(unique(trimws(as.character(c(global_genes, analysis_genes)))))
  }

  unique(trimws(as.character(custom_label_genes)))
}

match_custom_genes <- function(dat, custom_genes, match_columns) {
  # 按多个候选列匹配自定义基因；每个输入基因最多保留第一个匹配行。
  match_columns <- intersect(match_columns, colnames(dat))
  if (length(match_columns) == 0) {
    return(dat[0, , drop = FALSE])
  }

  custom_genes <- unique(trimws(as.character(custom_genes)))
  custom_genes <- custom_genes[!is.na(custom_genes) & custom_genes != ""]
  if (length(custom_genes) == 0) {
    return(dat[0, , drop = FALSE])
  }

  matched_rows <- integer(0)
  for (custom_gene in custom_genes) {
    matched_index <- integer(0)
    for (match_column in match_columns) {
      column_values <- trimws(as.character(dat[[match_column]]))
      matched_index <- which(column_values == custom_gene)
      if (length(matched_index) > 0) {
        break
      }
    }

    if (length(matched_index) > 0) {
      matched_rows <- c(matched_rows, matched_index[1])
    }
  }

  matched_rows <- unique(matched_rows)
  if (length(matched_rows) == 0) {
    return(dat[0, , drop = FALSE])
  }

  dat[matched_rows, , drop = FALSE]
}

get_custom_gene_label_data <- function(
    plot_data,
    symbol_column,
    custom_label_genes,
    match_columns) {
  # 只标注用户指定的基因；不在当前分析中的基因会被自然跳过。
  if (!symbol_column %in% colnames(plot_data)) {
    return(plot_data[0, , drop = FALSE])
  }

  label_data_list <- lapply(unique(plot_data$Analysis_Name), function(analysis_name) {
    custom_genes <- get_custom_genes_for_analysis(analysis_name, custom_label_genes)
    custom_genes <- custom_genes[!is.na(custom_genes) & custom_genes != ""]

    if (length(custom_genes) == 0) {
      return(plot_data[0, , drop = FALSE])
    }

    dat <- plot_data[
      plot_data$Analysis_Name == analysis_name &
        plot_data$Regulation %in% c("Up", "Down"),
      ,
      drop = FALSE
    ]

    dat <- match_custom_genes(dat, custom_genes, match_columns)
    dat$Gene_Label <- trimws(as.character(dat[[symbol_column]]))
    dat[!is.na(dat$Gene_Label) & dat$Gene_Label != "", , drop = FALSE]
  })

  label_data <- do.call(rbind, label_data_list)
  rownames(label_data) <- NULL
  label_data
}

get_top_gene_label_data <- function(
    plot_data,
    symbol_column,
    p_value_column,
    top_up_n = 5,
    top_down_n = 5) {
  # 每个分析设计分别取Up和Down中P值最小的基因用于图中标注。
  if (!symbol_column %in% colnames(plot_data)) {
    return(plot_data[0, , drop = FALSE])
  }

  label_data_list <- lapply(unique(plot_data$Analysis_Name), function(analysis_name) {
    dat <- plot_data[
      plot_data$Analysis_Name == analysis_name &
        plot_data$Regulation %in% c("Up", "Down"),
      ,
      drop = FALSE
    ]

    dat$Gene_Label <- trimws(as.character(dat[[symbol_column]]))
    dat <- dat[!is.na(dat$Gene_Label) & dat$Gene_Label != "", , drop = FALSE]

    up_dat <- dat[dat$Regulation == "Up", , drop = FALSE]
    down_dat <- dat[dat$Regulation == "Down", , drop = FALSE]

    if (nrow(up_dat) > 0) {
      up_dat <- up_dat[
        order(up_dat[[p_value_column]], -abs(up_dat$logFC)),
        ,
        drop = FALSE
      ]
      up_dat <- up_dat[seq_len(min(top_up_n, nrow(up_dat))), , drop = FALSE]
    }

    if (nrow(down_dat) > 0) {
      down_dat <- down_dat[
        order(down_dat[[p_value_column]], -abs(down_dat$logFC)),
        ,
        drop = FALSE
      ]
      down_dat <- down_dat[seq_len(min(top_down_n, nrow(down_dat))), , drop = FALSE]
    }

    rbind(up_dat, down_dat)
  })

  label_data <- do.call(rbind, label_data_list)
  rownames(label_data) <- NULL
  label_data
}

count_status <- function(status_counts, status_name) {
  # 安全读取table中的计数；缺失类别返回0，避免summary里出现NA。
  value <- status_counts[status_name]
  if (is.na(value)) {
    return(0L)
  }

  as.integer(value)
}
