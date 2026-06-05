# 绘图脚本公共配置和函数
#
# 本文件保存00、03、04号绘图脚本共用的视觉风格配置和基础函数。
# 具体图片的样本选择、分析组合、PDF尺寸、专属排版参数仍放在各绘图脚本头部。


# 0. 共用绘图风格配置 ----------------------------------------------------------

# 全部绘图文字统一使用Helvetica黑色粗体，便于不同图之间保持一致。
TEXT_FONT_FAMILY <- "Helvetica"
TEXT_FONT_FACE <- "bold"
TEXT_COLOR <- "black"
BASE_FONT_SIZE <- 12

# ggplot类图形的坐标轴/边框线宽；ComplexHeatmap热图会在脚本内单独放大。
AXIS_LINE_WIDTH <- 0.8

# 传统火山图和多组火山图共用的显著基因点样式。
UP_COLOR <- "#D73027"
DOWN_COLOR <- "#2166AC"
POINT_SIZE <- 3.2
POINT_ALPHA <- 0.60

# 绘图脚本同步输出PNG，方便插入Markdown文档。
PNG_DPI <- 300

# 图中展示基因名时优先使用的列；若该列不存在则不会标注基因名。
TOP_GENE_SYMBOL_COLUMN <- "Symbol"

# 自定义标注基因可按这些列匹配；图中展示文字仍优先使用TOP_GENE_SYMBOL_COLUMN。
CUSTOM_LABEL_MATCH_COLUMNS <- c("Symbol", "Feature_ID", "GeneID", "Ensembl", "Entrez")

# 未配置CUSTOM_LABEL_GENES时，默认每个分析标注Up 5个和Down 5个Top基因。
TOP_UP_LABEL_N <- 5
TOP_DOWN_LABEL_N <- 5

# Top基因文字和引线样式；03号和04号火山图共用。
TOP_GENE_LABEL_FONT_SIZE <- 3.4
TOP_GENE_LABEL_FONT_FACE <- TEXT_FONT_FACE
TOP_GENE_LABEL_BOX_PADDING <- 0.30
TOP_GENE_LABEL_POINT_PADDING <- 0.20
TOP_GENE_LABEL_SEGMENT_WIDTH <- 0.28
TOP_GENE_LABEL_MAX_OVERLAPS <- Inf
TOP_GENE_LABEL_FORCE <- 2.0
TOP_GENE_LABEL_FORCE_PULL <- 0.25

# 标注文字沿用Up/Down颜色，但略微加深，避免在PDF里显得发灰。
TOP_GENE_LABEL_COLOR_DARKEN <- 0.78


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

darken_color <- function(color, fraction = 0.80) {
  # 按原色等比例加深颜色，主要用于彩色文字和边框。
  # fraction越小颜色越深；fraction为1时保持原色。
  rgb_matrix <- grDevices::col2rgb(color) / 255
  rgb_matrix <- pmax(pmin(rgb_matrix * fraction, 1), 0)

  grDevices::rgb(
    red = rgb_matrix[1, ],
    green = rgb_matrix[2, ],
    blue = rgb_matrix[3, ],
    names = names(color)
  )
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

sync_de_thresholds_from_script <- function(
    script_file,
    p_value_column,
    p_value_cutoff,
    logfc_cutoff,
    sync = TRUE) {
  # 从01号差异分析脚本同步火山图使用的P值列名、P值阈值和logFC阈值。
  # sync为FALSE时直接返回当前脚本头部手动配置的阈值。
  if (sync) {
    p_value_column <- read_scalar_config(
      script_file,
      "P_VALUE_COLUMN",
      p_value_column
    )
    p_value_cutoff <- read_scalar_config(
      script_file,
      "P_VALUE_CUTOFF",
      p_value_cutoff
    )
    logfc_cutoff <- read_scalar_config(
      script_file,
      "LOGFC_CUTOFF",
      logfc_cutoff
    )
  }

  stopifnot(is.character(p_value_column))
  stopifnot(length(p_value_column) == 1)
  stopifnot(is.numeric(p_value_cutoff))
  stopifnot(is.numeric(logfc_cutoff))
  stopifnot(p_value_cutoff > 0)
  stopifnot(logfc_cutoff > 0)

  list(
    p_value_column = p_value_column,
    p_value_cutoff = p_value_cutoff,
    logfc_cutoff = logfc_cutoff
  )
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

get_display_labels <- function(sample_info, label_column, fallback_column = "Sample_ID") {
  # 提取图中展示用的样本名；若指定列为空，则回退到Sample_ID。
  stopifnot(label_column %in% colnames(sample_info))
  stopifnot(fallback_column %in% colnames(sample_info))

  display_labels <- trimws(as.character(sample_info[[label_column]]))
  empty_label_index <- display_labels == "" | is.na(display_labels)
  display_labels[empty_label_index] <- sample_info[[fallback_column]][empty_label_index]

  stopifnot(!any(duplicated(display_labels)))
  display_labels
}

get_label_line_count <- function(label_text) {
  # 统计换行标签的实际行数，用于动态调整标签框或PDF尺寸。
  length(strsplit(label_text, "\n", fixed = TRUE)[[1]])
}

get_named_brewer_palette <- function(levels, palette = "Set2") {
  # 为离散分组生成命名颜色；少于3个水平时仍按RColorBrewer要求取3色后截取。
  levels <- sort(unique(as.character(levels)))
  colors <- RColorBrewer::brewer.pal(
    max(3, length(levels)),
    palette
  )[seq_along(levels)]
  names(colors) <- levels

  colors
}

prepare_sample_correlation <- function(
    expr_matrix,
    correlation_method = "pearson",
    clustering_method = "complete") {
  # 对样本表达矩阵计算样本相关性和层次聚类。
  # 输入矩阵要求行为基因、列为样本；内部使用log2(x + 1)并去除零方差基因。
  expr_for_correlation <- log2(expr_matrix + 1)

  gene_sd <- apply(expr_for_correlation, 1, sd, na.rm = TRUE)
  expr_for_correlation <- expr_for_correlation[
    is.finite(gene_sd) & gene_sd > 0,
    ,
    drop = FALSE
  ]
  stopifnot(nrow(expr_for_correlation) > 1)

  cor_matrix <- cor(
    expr_for_correlation,
    method = correlation_method,
    use = "pairwise.complete.obs"
  )
  stopifnot(!any(is.na(cor_matrix)))

  row_distance <- as.dist(1 - cor_matrix)
  sample_hclust <- hclust(row_distance, method = clustering_method)

  list(
    expr_for_correlation = expr_for_correlation,
    cor_matrix = cor_matrix,
    sample_hclust = sample_hclust
  )
}

get_deg_file_info <- function(table_root, deg_dir_name = "DEG") {
  # 查找01号差异分析脚本输出的all_genes.csv。
  # 当前结构为tables/<analysis_name>/DEG/all_genes.csv；
  # 同时兼容历史tables/<analysis_name>/DEG/csv/all_genes.csv。
  all_gene_files <- list.files(
    table_root,
    pattern = "^all_genes[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  all_gene_files <- all_gene_files[
    basename(dirname(all_gene_files)) == deg_dir_name |
      (
        basename(dirname(all_gene_files)) == "csv" &
          basename(dirname(dirname(all_gene_files))) == deg_dir_name
      )
  ]

  get_analysis_name <- function(file_name) {
    if (basename(dirname(file_name)) == "csv") {
      return(basename(dirname(dirname(dirname(file_name)))))
    }

    basename(dirname(dirname(file_name)))
  }

  if (exists("prefer_report_csv_files", mode = "function")) {
    all_gene_files <- prefer_report_csv_files(all_gene_files, get_analysis_name)
  }

  stopifnot(length(all_gene_files) > 0)

  file_info <- data.frame(
    Analysis_Name = vapply(all_gene_files, get_analysis_name, character(1)),
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

get_selected_analysis_names <- function(file_info, analyses_to_plot) {
  # 根据脚本头部配置选择需要绘图的分析设计。
  # analyses_to_plot为"all"时，自动使用全部可用的DEG结果。
  if (identical(analyses_to_plot, "all")) {
    return(file_info$Analysis_Name)
  }

  selected_analyses <- analyses_to_plot
  missing_analyses <- setdiff(selected_analyses, file_info$Analysis_Name)
  if (length(missing_analyses) > 0) {
    stop(
      "No all_genes.csv file was found for: ",
      paste(missing_analyses, collapse = ", ")
    )
  }

  selected_analyses
}

read_deg_result <- function(file_info, analysis_name) {
  # 读取指定分析设计对应的all_genes.csv。
  # file_info来自get_deg_file_info()，可避免在主脚本里反复拼路径。
  file_index <- match(analysis_name, file_info$Analysis_Name)
  if (is.na(file_index)) {
    stop("No all_genes.csv file was found for: ", analysis_name)
  }

  read.csv(
    if (exists("resolve_report_csv_file", mode = "function")) {
      resolve_report_csv_file(file_info$All_Genes_File[file_index])
    } else {
      file_info$All_Genes_File[file_index]
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
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

get_volcano_label_data <- function(
    plot_data,
    custom_label_genes,
    symbol_column,
    match_columns,
    p_value_column,
    top_up_n = 5,
    top_down_n = 5) {
  # 火山图基因标注的统一入口：
  # 1. 如果用户配置了CUSTOM_LABEL_GENES，只标注用户指定基因；
  # 2. 如果未配置，则自动标注每个分析的Top Up/Down基因。
  if (has_custom_label_genes(custom_label_genes)) {
    return(get_custom_gene_label_data(
      plot_data = plot_data,
      symbol_column = symbol_column,
      custom_label_genes = custom_label_genes,
      match_columns = match_columns
    ))
  }

  get_top_gene_label_data(
    plot_data = plot_data,
    symbol_column = symbol_column,
    p_value_column = p_value_column,
    top_up_n = top_up_n,
    top_down_n = top_down_n
  )
}

get_regulation_label_colors <- function(
    label_data,
    up_color,
    down_color,
    darken_fraction = 0.78) {
  # 根据Up/Down状态生成基因标注文字和引线颜色。
  # 颜色沿用点的红蓝配色，但略微加深，便于PDF中阅读。
  if (nrow(label_data) == 0) {
    return(character(0))
  }

  ifelse(
    as.character(label_data$Regulation) == "Up",
    darken_color(up_color, darken_fraction),
    darken_color(down_color, darken_fraction)
  )
}

add_volcano_gene_label_layer <- function(
    plot,
    label_data,
    label_colors,
    use_gg_repel,
    text_family,
    fontface,
    font_size,
    box_padding,
    point_padding,
    segment_width,
    force,
    force_pull,
    max_overlaps,
    seed = 1,
    nudge_y = NULL,
    fallback_vjust = -0.8) {
  # 为传统火山图和多组火山图添加基因名标注。
  # 优先使用ggrepel自动避让；如果未安装ggrepel，则退回普通geom_text。
  if (nrow(label_data) == 0) {
    return(plot)
  }

  if (use_gg_repel) {
    repel_args <- list(
      data = label_data,
      mapping = ggplot2::aes(label = Gene_Label),
      family = text_family,
      fontface = fontface,
      size = font_size,
      box.padding = box_padding,
      point.padding = point_padding,
      segment.size = segment_width,
      segment.alpha = 0.65,
      segment.color = label_colors,
      min.segment.length = 0,
      force = force,
      force_pull = force_pull,
      max.overlaps = max_overlaps,
      seed = seed,
      color = label_colors,
      show.legend = FALSE
    )

    if (!is.null(nudge_y)) {
      repel_args$nudge_y <- nudge_y
    }

    return(plot + do.call(ggrepel::geom_text_repel, repel_args))
  }

  plot + ggplot2::geom_text(
    data = label_data,
    ggplot2::aes(label = Gene_Label),
    family = text_family,
    fontface = fontface,
    size = font_size,
    vjust = fallback_vjust,
    color = label_colors,
    check_overlap = TRUE,
    show.legend = FALSE
  )
}

count_status <- function(status_counts, status_name) {
  # 安全读取table中的计数；缺失类别返回0，避免summary里出现NA。
  value <- status_counts[status_name]
  if (is.na(value)) {
    return(0L)
  }

  as.integer(value)
}

get_plot_output_paths <- function(plot_file) {
  # 接受旧式 output_dir/name.pdf，也接受已经位于pdf/png目录中的路径。
  # 实际保存时统一放到 output_dir/pdf/name.pdf 与 output_dir/png/name.png。
  plot_file <- as.character(plot_file)
  plot_root <- dirname(plot_file)
  plot_format_dir <- basename(plot_root)

  if (plot_format_dir %in% c("pdf", "png")) {
    plot_root <- dirname(plot_root)
  }

  file_stem <- tools::file_path_sans_ext(basename(plot_file))
  list(
    root = plot_root,
    pdf = file.path(plot_root, "pdf", paste0(file_stem, ".pdf")),
    png = file.path(plot_root, "png", paste0(file_stem, ".png"))
  )
}

get_png_file_from_pdf <- function(pdf_file) {
  get_plot_output_paths(pdf_file)$png
}

save_grid_pdf_png <- function(pdf_file, width, height, draw_fun, png_dpi = PNG_DPI) {
  # 使用同一套绘图函数分别输出矢量PDF和高清PNG。
  plot_paths <- get_plot_output_paths(pdf_file)
  dir.create(dirname(plot_paths$pdf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_paths$png), recursive = TRUE, showWarnings = FALSE)

  Cairo::CairoPDF(
    file = plot_paths$pdf,
    width = width,
    height = height,
    bg = "white"
  )
  draw_fun()
  invisible(dev.off())

  Cairo::Cairo(
    file = plot_paths$png,
    type = "png",
    width = width,
    height = height,
    units = "in",
    dpi = png_dpi,
    bg = "white",
    canvas = "white"
  )
  draw_fun()
  invisible(dev.off())

  stopifnot(file.exists(plot_paths$pdf))
  stopifnot(file.exists(plot_paths$png))

  invisible(list(pdf_file = plot_paths$pdf, png_file = plot_paths$png))
}

save_ggplot_pdf_png <- function(plot, pdf_file, width, height, png_dpi = PNG_DPI) {
  # 保存ggplot图形为PDF和PNG。
  save_grid_pdf_png(
    pdf_file = pdf_file,
    width = width,
    height = height,
    png_dpi = png_dpi,
    draw_fun = function() print(plot)
  )
}

save_ggplot_pdf <- function(plot, pdf_file, width, height) {
  # 兼容旧调用名；现在同步输出同名PNG。
  save_ggplot_pdf_png(
    plot = plot,
    pdf_file = pdf_file,
    width = width,
    height = height
  )
}
