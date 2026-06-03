# GSE114012多组差异表达火山图
#
# 读取01号limma脚本生成的all_genes.csv结果，按配置将多个差异分析设计
# 合并到一张多组火山图中。每个差异分析设计占一个分面；
# 横轴为-log10(P值)，纵轴为logFC。图中只展示达到阈值的Up和Down基因；
# NS基因只参与终端汇总，不在图中绘制。阈值与01号脚本保持一致。
#
# 参考图形思路：
# 1. R2Omics多组火山图：按分组展示显著性和logFC分布，显著基因用颜色突出。
# 2. openbiox/Bizard多组火山图：多个比较对在同一图中并排展示。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

DE_SCRIPT_FILE <- "scripts/GSE114012/01_limma_differential_expression.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "multiple_volcano")

# 是否从01号差异分析脚本中同步显著性阈值。
# 建议保持TRUE，保证多组火山图阈值与差异分析结果完全一致。
SYNC_THRESHOLDS_FROM_01_SCRIPT <- TRUE

# 只有在SYNC_THRESHOLDS_FROM_01_SCRIPT为FALSE时，下面三项才作为手动阈值使用。
P_VALUE_COLUMN <- "P.Value"
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

# 多套多组火山图方案。
# 新增方案时，只需要按下面格式增加一行；顺序会影响图中分组排列。
# 每个方案会单独输出一个PDF：
# results/ngs/GSE114012/plots/multiple_volcano/<scheme_name>/multiple_volcano_plot.pdf
MULTIPLE_VOLCANO_SCHEMES <- list(
  CRC_LRC_core = c("DLD1", "HCT15", "HT55", "RKO", "SW48", "SW948", "ALL", "DLD1_HCT15", "DLD1_HCT15_SW48"),
  DLD1_HCT15_SW48 = c("DLD1", "HCT15", "SW48"),
  DLD1_HCT15 = c("DLD1", "HCT15"),
  HT55_SW948 = c("HT55", "SW948"),
  SW948_RKO = c("SW948", "RKO")
)

# 运行哪些多组火山图方案。
# 可设为names(MULTIPLE_VOLCANO_SCHEMES)运行全部；也可只写部分方案名。
SCHEMES_TO_RUN <- names(MULTIPLE_VOLCANO_SCHEMES)
# SCHEMES_TO_RUN <- c("CRC_LRC_core")

# 重跑时清理整个multiple_volcano目录，保证最终只保留当前配置生成的方案图。
CLEAN_MULTIPLE_VOLCANO_ROOT <- TRUE

# 重跑时清理当前方案目录内旧PDF，避免旧文件残留。
OVERWRITE_SCHEME_OUTPUT <- TRUE

# 颜色和点样式与03号火山图保持一致。
UP_COLOR <- "#D73027"
DOWN_COLOR <- "#2166AC"
POINT_SIZE <- 3.2
POINT_ALPHA <- 0.60

# 每组标注Top显著基因。Top排序优先按P值，其次按logFC幅度。
# 若CUSTOM_LABEL_GENES为空，则自动标注每组Up 5个和Down 5个Top基因。
# 若CUSTOM_LABEL_GENES不为空，则只标注这里配置的基因，不再自动标注Top基因。
# 用法1：所有分析都标注同一批基因。
# CUSTOM_LABEL_GENES <- c("MYH15", "CSF1R", "SRSF1")
# 用法2：不同分析标注不同基因；all表示所有分析都标注。
# CUSTOM_LABEL_GENES <- list(all = c("MYH15"), DLD1 = c("ABCG1", "SRSF1"))
CUSTOM_LABEL_GENES <- character(0)

TOP_GENE_SYMBOL_COLUMN <- "Symbol"
# 自定义标注基因可按下面这些列匹配；图中展示文字仍优先使用TOP_GENE_SYMBOL_COLUMN。
CUSTOM_LABEL_MATCH_COLUMNS <- c("Symbol", "Feature_ID", "GeneID", "Ensembl", "Entrez")
TOP_UP_LABEL_N <- 5
TOP_DOWN_LABEL_N <- 5
TOP_GENE_LABEL_FONT_SIZE <- 3.4
TOP_GENE_LABEL_FONT_FACE <- "bold"
TOP_GENE_LABEL_BOX_PADDING <- 0.30
TOP_GENE_LABEL_POINT_PADDING <- 0.20
TOP_GENE_LABEL_SEGMENT_WIDTH <- 0.28
TOP_GENE_LABEL_MAX_OVERLAPS <- Inf
TOP_GENE_LABEL_NUDGE_Y <- 0.80
TOP_GENE_LABEL_FORCE <- 2.0
TOP_GENE_LABEL_FORCE_PULL <- 0.25

# 每组中间的彩色标签；参考CNScolor/sciRcolor中HEX_color[[30]]的顶刊离散配色。
# 未列出的分析会自动分配颜色；若出现重复颜色，脚本会自动替换。
GROUP_LABEL_COLORS <- c(
  "DLD1" = "#c70008",
  "HCT15" = "#eb7400",
  "HT55" = "#006b3f",
  "RKO" = "#0052a2",
  "SW48" = "#0090c1",
  "SW948" = "#b9b453",
  "ALL" = "#805190",
  "DLD1_HCT15" = "#8e5640",
  "DLD1_HCT15_SW48" = "#d6779f"
)
GROUP_LABEL_ALPHA <- 0.24
GROUP_LABEL_WRAP_WIDTH <- 10
GROUP_LABEL_FONT_SIZE <- 5.3
GROUP_LABEL_LINE_HEIGHT <- 1.18
# 组名框位于logFC阈值内侧：当前阈值0.5时为±0.4；阈值1时为±0.9。
GROUP_LABEL_BOX_LOGFC_GAP <- 0.10
GROUP_LABEL_BOX_MIN_FRACTION <- 0.80
GROUP_LABEL_FONT_MIN_SIZE <- 3.9
GROUP_LABEL_FONT_LINE_SHRINK <- 0.55
GROUP_LABEL_LINE_HEIGHT_MIN <- 0.95
GROUP_LABEL_LINE_HEIGHT_SHRINK <- 0.08
GROUP_LABEL_BOX_X_MARGIN_FRACTION <- 0.05
GROUP_LABEL_BORDER_WIDTH <- 0.9

# 边框样式。
AXIS_LINE_WIDTH <- 0.8
PANEL_SPACING_X_MM <- 4.6
LEGEND_TOP_MARGIN_PT <- 68
X_AXIS_PADDING_FRACTION <- 0.12
X_AXIS_PADDING_MIN <- 0.22
X_AXIS_LEFT_MIN_FRACTION <- 0.78
Y_AXIS_PADDING_FRACTION <- 0.07
LABEL_BOX_HEIGHT_PDF_SCALE <- 1.8

# PDF和字体设置。文件名固定为multiple_volcano_plot.pdf；方案名由目录体现。
BASE_PDF_HEIGHT <- 6.0
GROUP_WIDTH_INCH <- 1.78
LEGEND_WIDTH_INCH <- 1.18
MIN_PDF_WIDTH <- 7.2
MAX_PDF_WIDTH <- 20.0
MAX_PDF_HEIGHT <- 8.2
TEXT_FONT_FAMILY <- "Helvetica"
BASE_FONT_SIZE <- 12
TEXT_FONT_FACE <- "bold"
TEXT_COLOR <- "black"
LEGEND_TEXT_SIZE <- 13.5
LEGEND_POINT_SIZE_SCALE <- 1.45
LEGEND_KEY_SIZE_MM <- 6.6

options(width = 200)


# 1. 加载包 -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(Cairo)
})

USE_GG_REPEL <- requireNamespace("ggrepel", quietly = TRUE)


# 2. 常用函数 -----------------------------------------------------------------

sanitize_file_name <- function(x) {
  # 文件夹名保留字母、数字、下划线、点和短横线，避免不同系统下路径出错。
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

read_scalar_config <- function(script_file, variable_name, default_value) {
  # 静态读取01号脚本顶部配置，不source脚本，避免重复运行差异分析。
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

wrap_label_by_underscore <- function(x, width = 10) {
  # 分组名较长时按下划线拆行，避免底部标签重叠。
  x <- as.character(x)

  vapply(x, function(label) {
    parts <- strsplit(label, "_", fixed = TRUE)[[1]]

    if (length(parts) == 1) {
      return(label)
    }

    lines <- character(0)
    current_line <- ""

    for (part in parts) {
      token <- if (current_line == "") part else paste0("_", part)
      candidate <- paste0(current_line, token)

      if (nchar(candidate) <= width || current_line == "") {
        current_line <- candidate
      } else {
        lines <- c(lines, current_line)
        current_line <- part
      }
    }

    paste(c(lines, current_line), collapse = "\n")
  }, character(1))
}

get_group_label_colors <- function(selected_analyses) {
  # 优先使用脚本头部指定的颜色；新增分析名时自动补充一组可区分颜色。
  group_colors <- GROUP_LABEL_COLORS[selected_analyses]
  missing_index <- is.na(group_colors)
  duplicated_index <- duplicated(group_colors) & !is.na(group_colors)

  if (any(missing_index | duplicated_index)) {
    replacement_colors <- grDevices::hcl.colors(
      sum(missing_index | duplicated_index),
      palette = "Dark 3"
    )
    group_colors[missing_index | duplicated_index] <- replacement_colors
  }

  if (any(duplicated(group_colors))) {
    group_colors <- grDevices::hcl.colors(
      length(selected_analyses),
      palette = "Dark 3"
    )
  }

  names(group_colors) <- selected_analyses
  group_colors
}

get_label_line_count <- function(label_text) {
  length(strsplit(label_text, "\n", fixed = TRUE)[[1]])
}

get_group_label_box_y_limits <- function() {
  # 组名框放在显著logFC阈值内侧：阈值为0.5时为±0.4，阈值为1时为±0.9。
  box_half_height <- max(
    LOGFC_CUTOFF - GROUP_LABEL_BOX_LOGFC_GAP,
    LOGFC_CUTOFF * GROUP_LABEL_BOX_MIN_FRACTION
  )
  box_half_height <- min(box_half_height, LOGFC_CUTOFF)

  c(-box_half_height, box_half_height)
}

get_group_label_text_style <- function(max_line_count) {
  # 多行组名通过字体和行距轻微动态收紧，以适配阈值内侧的彩色外框。
  line_extra <- max(max_line_count - 1, 0)

  list(
    font_size = max(
      GROUP_LABEL_FONT_SIZE - line_extra * GROUP_LABEL_FONT_LINE_SHRINK,
      GROUP_LABEL_FONT_MIN_SIZE
    ),
    line_height = max(
      GROUP_LABEL_LINE_HEIGHT - line_extra * GROUP_LABEL_LINE_HEIGHT_SHRINK,
      GROUP_LABEL_LINE_HEIGHT_MIN
    )
  )
}

has_custom_label_genes <- function() {
  # character向量或list均可；为空时自动使用top基因。
  if (is.list(CUSTOM_LABEL_GENES)) {
    return(length(unlist(CUSTOM_LABEL_GENES, use.names = FALSE)) > 0)
  }

  length(CUSTOM_LABEL_GENES) > 0
}

get_custom_genes_for_analysis <- function(analysis_name) {
  if (!has_custom_label_genes()) {
    return(character(0))
  }

  if (is.list(CUSTOM_LABEL_GENES)) {
    global_genes <- if ("all" %in% names(CUSTOM_LABEL_GENES)) {
      CUSTOM_LABEL_GENES[["all"]]
    } else {
      character(0)
    }
    analysis_genes <- if (analysis_name %in% names(CUSTOM_LABEL_GENES)) {
      CUSTOM_LABEL_GENES[[analysis_name]]
    } else {
      character(0)
    }

    return(unique(trimws(as.character(c(global_genes, analysis_genes)))))
  }

  unique(trimws(as.character(CUSTOM_LABEL_GENES)))
}

match_custom_genes <- function(dat, custom_genes) {
  # 自定义基因可按Symbol、Feature_ID等配置列匹配；展示文本仍使用TOP_GENE_SYMBOL_COLUMN。
  match_columns <- intersect(CUSTOM_LABEL_MATCH_COLUMNS, colnames(dat))
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

get_custom_gene_label_data <- function(plot_data) {
  # 只标注用户指定的基因；不在当前分析中的基因会被自然跳过。
  if (!TOP_GENE_SYMBOL_COLUMN %in% colnames(plot_data)) {
    return(plot_data[0, , drop = FALSE])
  }

  label_data_list <- lapply(unique(plot_data$Analysis_Name), function(analysis_name) {
    custom_genes <- get_custom_genes_for_analysis(analysis_name)
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

    dat$Gene_Label <- trimws(as.character(dat[[TOP_GENE_SYMBOL_COLUMN]]))
    dat <- match_custom_genes(dat, custom_genes)
    dat$Gene_Label <- trimws(as.character(dat[[TOP_GENE_SYMBOL_COLUMN]]))
    dat[!is.na(dat$Gene_Label) & dat$Gene_Label != "", , drop = FALSE]
  })

  label_data <- do.call(rbind, label_data_list)
  rownames(label_data) <- NULL
  label_data
}

get_top_gene_label_data <- function(plot_data) {
  # 每个分析设计分别取Up 5个和Down 5个最显著基因用于图中标注。
  if (!TOP_GENE_SYMBOL_COLUMN %in% colnames(plot_data)) {
    return(plot_data[0, , drop = FALSE])
  }

  label_data_list <- lapply(unique(plot_data$Analysis_Name), function(analysis_name) {
    dat <- plot_data[
      plot_data$Analysis_Name == analysis_name &
        plot_data$Regulation %in% c("Up", "Down"),
      ,
      drop = FALSE
    ]

    dat$Gene_Label <- trimws(as.character(dat[[TOP_GENE_SYMBOL_COLUMN]]))
    dat <- dat[!is.na(dat$Gene_Label) & dat$Gene_Label != "", , drop = FALSE]

    up_dat <- dat[dat$Regulation == "Up", , drop = FALSE]
    down_dat <- dat[dat$Regulation == "Down", , drop = FALSE]

    if (nrow(up_dat) > 0) {
      up_dat <- up_dat[
        order(up_dat[[P_VALUE_COLUMN]], -abs(up_dat$logFC)),
        ,
        drop = FALSE
      ]
      up_dat <- up_dat[seq_len(min(TOP_UP_LABEL_N, nrow(up_dat))), , drop = FALSE]
    }

    if (nrow(down_dat) > 0) {
      down_dat <- down_dat[
        order(down_dat[[P_VALUE_COLUMN]], -abs(down_dat$logFC)),
        ,
        drop = FALSE
      ]
      down_dat <- down_dat[seq_len(min(TOP_DOWN_LABEL_N, nrow(down_dat))), , drop = FALSE]
    }

    rbind(up_dat, down_dat)
  })

  label_data <- do.call(rbind, label_data_list)
  rownames(label_data) <- NULL
  label_data
}

get_deg_file_info <- function(table_root) {
  # 01号脚本将结果保存为tables/<analysis_name>/DEG/all_genes.csv。
  all_gene_files <- list.files(
    table_root,
    pattern = "^all_genes[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  all_gene_files <- all_gene_files[
    basename(dirname(all_gene_files)) == "DEG"
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

prepare_group_data <- function(dat, analysis_name) {
  # 显著性判定与01号脚本保持一致：P值小于阈值且abs(logFC)大于阈值。
  stopifnot("logFC" %in% colnames(dat))
  stopifnot(P_VALUE_COLUMN %in% colnames(dat))

  dat$logFC <- as.numeric(dat$logFC)
  dat[[P_VALUE_COLUMN]] <- as.numeric(dat[[P_VALUE_COLUMN]])

  valid_index <- is.finite(dat$logFC) &
    is.finite(dat[[P_VALUE_COLUMN]]) &
    !is.na(dat[[P_VALUE_COLUMN]])

  dat <- dat[valid_index, , drop = FALSE]
  stopifnot(nrow(dat) > 0)

  positive_p <- dat[[P_VALUE_COLUMN]][dat[[P_VALUE_COLUMN]] > 0]
  stopifnot(length(positive_p) > 0)

  min_positive_p <- min(positive_p, na.rm = TRUE)
  safe_p <- dat[[P_VALUE_COLUMN]]
  safe_p[safe_p <= 0] <- min_positive_p * 0.1

  dat$Neg_Log10_P <- -log10(safe_p)
  dat$Analysis_Name <- analysis_name
  dat$Regulation <- "NS"
  dat$Regulation[
    dat$logFC > LOGFC_CUTOFF &
      dat[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF
  ] <- "Up"
  dat$Regulation[
    dat$logFC < -LOGFC_CUTOFF &
      dat[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF
  ] <- "Down"

  dat$Regulation <- factor(
    dat$Regulation,
    levels = c("NS", "Down", "Up")
  )

  dat
}

get_group_layout_data <- function(plot_data, selected_analyses) {
  # 计算每个分面自己的横轴范围，并给中心组名框预留统一高度。
  # 组名框在横轴两侧保留固定比例空白，保证各分面标签框的物理宽度一致。
  group_labels <- wrap_label_by_underscore(
    selected_analyses,
    width = GROUP_LABEL_WRAP_WIDTH
  )
  names(group_labels) <- selected_analyses

  threshold_x <- -log10(P_VALUE_CUTOFF)
  max_line_count <- max(vapply(
    group_labels,
    get_label_line_count,
    integer(1)
  ))

  label_style <- get_group_label_text_style(max_line_count)
  label_box_y_limits <- get_group_label_box_y_limits()
  box_height <- diff(label_box_y_limits)

  do.call(rbind, lapply(selected_analyses, function(analysis_name) {
    dat <- plot_data[
      plot_data$Analysis_Name == analysis_name &
        plot_data$Regulation %in% c("Up", "Down"),
      ,
      drop = FALSE
    ]

    if (nrow(dat) == 0) {
      dat <- plot_data[plot_data$Analysis_Name == analysis_name, , drop = FALSE]
    }

    x_max <- max(
      dat$Neg_Log10_P,
      threshold_x,
      na.rm = TRUE
    )

    x_range <- max(x_max - threshold_x, 0.5)
    x_padding <- max(
      x_range * X_AXIS_PADDING_FRACTION,
      X_AXIS_PADDING_MIN
    )
    x_limit_min <- max(
      threshold_x * X_AXIS_LEFT_MIN_FRACTION,
      threshold_x - x_padding
    )
    x_limit_max <- x_max + x_padding
    box_x_margin <- (x_limit_max - x_limit_min) *
      GROUP_LABEL_BOX_X_MARGIN_FRACTION

    data.frame(
      Analysis_Name = analysis_name,
      Label = group_labels[analysis_name],
      Label_X = (x_limit_min + box_x_margin + x_limit_max - box_x_margin) / 2,
      Label_Y = 0,
      Box_X_Min = x_limit_min + box_x_margin,
      Box_X_Max = x_limit_max - box_x_margin,
      Box_Y_Min = label_box_y_limits[1],
      Box_Y_Max = label_box_y_limits[2],
      X_Min = x_limit_min,
      X_Max = x_limit_max,
      X_Data_Min = threshold_x,
      X_Data_Max = x_max,
      X_Padding = x_padding,
      Box_X_Margin = box_x_margin,
      Label_Box_Height = box_height,
      Label_Font_Size = label_style$font_size,
      Label_Line_Height = label_style$line_height,
      stringsAsFactors = FALSE
    )
  }))
}

apply_group_label_display_y <- function(plot_data, group_layout) {
  # 图中纵坐标直接使用真实logFC；PDF高度会按组名框动态拉长。
  plot_data$Plot_LogFC <- plot_data$logFC

  plot_data
}

get_axis_info <- function(plot_data, group_layout) {
  # 纵坐标保持真实logFC；组名框只参与坐标范围和PDF高度计算。
  plotted_data <- plot_data[
    plot_data$Regulation %in% c("Up", "Down"),
    ,
    drop = FALSE
  ]

  true_y_min <- min(
    plotted_data$logFC,
    -LOGFC_CUTOFF,
    na.rm = TRUE
  )
  true_y_max <- max(
    plotted_data$logFC,
    LOGFC_CUTOFF,
    na.rm = TRUE
  )
  true_y_span <- max(true_y_max - true_y_min, 1)
  true_y_padding <- true_y_span * Y_AXIS_PADDING_FRACTION
  true_limits <- c(true_y_min - true_y_padding, true_y_max + true_y_padding)

  true_breaks <- pretty(true_limits, n = 7)
  true_breaks <- true_breaks[
    true_breaks >= true_limits[1] &
      true_breaks <= true_limits[2]
  ]
  true_breaks <- unique(c(true_breaks, -LOGFC_CUTOFF, 0, LOGFC_CUTOFF))
  true_breaks <- true_breaks[
    true_breaks >= true_limits[1] &
      true_breaks <= true_limits[2]
  ]
  true_breaks <- sort(unique(round(true_breaks, 6)))

  list(
    true_limits = true_limits,
    display_limits = true_limits,
    breaks = true_breaks,
    labels = format(true_breaks, trim = TRUE, scientific = FALSE)
  )
}

get_pdf_size <- function(group_count, display_y_limits, group_layout) {
  # 宽度随组数增加；高度按logFC范围和组名框高度等比例动态拉长。
  y_span <- diff(display_y_limits)
  label_box_height <- max(group_layout$Label_Box_Height, na.rm = TRUE)
  label_height_extra <- max(
    label_box_height - LOGFC_CUTOFF,
    0
  ) * LABEL_BOX_HEIGHT_PDF_SCALE

  pdf_height <- BASE_PDF_HEIGHT +
    min(max((y_span - 5) * 0.08, 0), 1.2) +
    label_height_extra
  pdf_height <- min(pdf_height, MAX_PDF_HEIGHT)

  pdf_width <- GROUP_WIDTH_INCH * group_count + LEGEND_WIDTH_INCH + 1.2
  pdf_width <- max(pdf_width, MIN_PDF_WIDTH)
  pdf_width <- max(pdf_width, pdf_height * 1.12)
  pdf_width <- min(pdf_width, MAX_PDF_WIDTH)

  list(
    width = pdf_width,
    height = pdf_height
  )
}

make_multiple_volcano_plot <- function(plot_data, group_layout, selected_analyses, axis_info) {
  plot_data$Analysis_Name <- factor(
    plot_data$Analysis_Name,
    levels = selected_analyses
  )
  group_layout$Analysis_Name <- factor(
    group_layout$Analysis_Name,
    levels = selected_analyses
  )

  point_data <- plot_data[
    plot_data$Regulation %in% c("Up", "Down"),
    ,
    drop = FALSE
  ]
  stopifnot(nrow(point_data) > 0)

  top_label_data <- if (has_custom_label_genes()) {
    get_custom_gene_label_data(plot_data)
  } else {
    get_top_gene_label_data(plot_data)
  }
  if (nrow(top_label_data) > 0) {
    top_label_data$Analysis_Name <- factor(
      top_label_data$Analysis_Name,
      levels = selected_analyses
    )
    top_label_data$Label_Nudge_Y <- ifelse(
      top_label_data$Regulation == "Up",
      TOP_GENE_LABEL_NUDGE_Y,
      -TOP_GENE_LABEL_NUDGE_Y
    )
  }

  x_axis_anchor_data <- data.frame(
    Analysis_Name = factor(
      rep(as.character(group_layout$Analysis_Name), each = 2),
      levels = selected_analyses
    ),
    Neg_Log10_P = c(
      as.vector(rbind(group_layout$X_Min, group_layout$X_Max))
    ),
    Plot_LogFC = 0,
    stringsAsFactors = FALSE
  )

  group_label_colors <- get_group_label_colors(selected_analyses)

  point_colors <- c(
    "Down" = DOWN_COLOR,
    "Up" = UP_COLOR
  )

  volcano_plot <- ggplot(
    point_data,
    aes(x = Neg_Log10_P, y = Plot_LogFC, color = Regulation)
  ) +
    geom_point(
      size = POINT_SIZE,
      alpha = POINT_ALPHA,
      shape = 16,
      stroke = 0
    ) +
    geom_blank(
      data = x_axis_anchor_data,
      aes(x = Neg_Log10_P, y = Plot_LogFC),
      inherit.aes = FALSE
    ) +
    facet_grid(
      cols = vars(Analysis_Name),
      scales = "free_x",
      space = "fixed"
    ) +
    scale_color_manual(
      values = point_colors,
      breaks = c("Up", "Down"),
      labels = c("Sig_Up", "Sig_Down")
    ) +
    scale_x_continuous(
      breaks = function(x) {
        x_breaks <- pretty(x, n = 4)
        x_breaks[x_breaks >= 0]
      },
      expand = expansion(mult = 0, add = 0)
    ) +
    scale_y_continuous(
      breaks = axis_info$breaks,
      labels = axis_info$labels,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(
      ylim = axis_info$display_limits,
      clip = "off"
    ) +
    labs(
      x = paste0("-log10(", P_VALUE_COLUMN, ")"),
      y = "log2 fold change",
      color = NULL
    ) +
    theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    theme(
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.22),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = TEXT_COLOR,
        fill = NA,
        linewidth = AXIS_LINE_WIDTH
      ),
      panel.spacing.x = grid::unit(PANEL_SPACING_X_MM, "mm"),
      axis.line = element_line(color = TEXT_COLOR, linewidth = AXIS_LINE_WIDTH),
      axis.text.x = element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.text.y = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.position = "right",
      legend.justification = "top",
      legend.box.just = "top",
      legend.text = element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        size = LEGEND_TEXT_SIZE
      ),
      legend.key = element_blank(),
      legend.key.height = grid::unit(LEGEND_KEY_SIZE_MM, "mm"),
      legend.key.width = grid::unit(LEGEND_KEY_SIZE_MM, "mm"),
      legend.box.spacing = grid::unit(8, "pt"),
      legend.margin = margin(LEGEND_TOP_MARGIN_PT, 0, 0, 4, unit = "pt"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      text = element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        family = TEXT_FONT_FAMILY
      ),
      plot.margin = margin(10, 12, 10, 10, unit = "pt")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = POINT_SIZE * LEGEND_POINT_SIZE_SCALE,
          alpha = POINT_ALPHA
        )
      )
    )

  for (analysis_name in selected_analyses) {
    current_label <- group_layout[
      group_layout$Analysis_Name == analysis_name,
      ,
      drop = FALSE
    ]

    volcano_plot <- volcano_plot +
      geom_rect(
        data = current_label,
        aes(
          xmin = Box_X_Min,
          xmax = Box_X_Max,
          ymin = Box_Y_Min,
          ymax = Box_Y_Max
        ),
        inherit.aes = FALSE,
        color = group_label_colors[analysis_name],
        fill = adjustcolor(
          group_label_colors[analysis_name],
          alpha.f = GROUP_LABEL_ALPHA
        ),
        linewidth = GROUP_LABEL_BORDER_WIDTH
      ) +
      geom_text(
        data = current_label,
        aes(x = Label_X, y = Label_Y, label = Label),
        inherit.aes = FALSE,
        family = TEXT_FONT_FAMILY,
        fontface = TEXT_FONT_FACE,
        size = current_label$Label_Font_Size[1],
        lineheight = current_label$Label_Line_Height[1],
        color = group_label_colors[analysis_name]
      )
  }

  if (nrow(top_label_data) > 0 && USE_GG_REPEL) {
    volcano_plot <- volcano_plot +
      ggrepel::geom_text_repel(
        data = top_label_data,
        aes(label = Gene_Label),
        family = TEXT_FONT_FAMILY,
        fontface = TOP_GENE_LABEL_FONT_FACE,
        size = TOP_GENE_LABEL_FONT_SIZE,
        box.padding = TOP_GENE_LABEL_BOX_PADDING,
        point.padding = TOP_GENE_LABEL_POINT_PADDING,
        segment.size = TOP_GENE_LABEL_SEGMENT_WIDTH,
        segment.alpha = 0.65,
        min.segment.length = 0,
        nudge_y = top_label_data$Label_Nudge_Y,
        force = TOP_GENE_LABEL_FORCE,
        force_pull = TOP_GENE_LABEL_FORCE_PULL,
        max.overlaps = TOP_GENE_LABEL_MAX_OVERLAPS,
        seed = 1,
        show.legend = FALSE
      )
  }

  if (nrow(top_label_data) > 0 && !USE_GG_REPEL) {
    volcano_plot <- volcano_plot +
      geom_text(
        data = top_label_data,
        aes(label = Gene_Label),
        family = TEXT_FONT_FAMILY,
        fontface = TOP_GENE_LABEL_FONT_FACE,
        size = TOP_GENE_LABEL_FONT_SIZE,
        vjust = ifelse(top_label_data$Regulation == "Up", -0.7, 1.2),
        check_overlap = TRUE,
        show.legend = FALSE
      )
  }

  volcano_plot
}

run_multiple_volcano_scheme <- function(scheme_name, selected_analyses, file_info) {
  # 单套方案的完整流程：读取数据、绘图、保存PDF、返回终端汇总表。
  stopifnot(is.character(scheme_name))
  stopifnot(length(scheme_name) == 1)
  stopifnot(scheme_name != "")
  stopifnot(length(selected_analyses) >= 2)

  selected_analyses <- unique(selected_analyses)

  missing_analyses <- setdiff(selected_analyses, file_info$Analysis_Name)
  if (length(missing_analyses) > 0) {
    stop(
      "No all_genes.csv file was found for scheme ",
      scheme_name,
      ": ",
      paste(missing_analyses, collapse = ", ")
    )
  }

  selected_file_info <- file_info[
    match(selected_analyses, file_info$Analysis_Name),
    ,
    drop = FALSE
  ]

  plot_data_list <- vector("list", length(selected_analyses))
  names(plot_data_list) <- selected_analyses

  for (analysis_index in seq_along(selected_analyses)) {
    analysis_name <- selected_analyses[analysis_index]
    all_genes_file <- selected_file_info$All_Genes_File[
      selected_file_info$Analysis_Name == analysis_name
    ]

    dat <- read.csv(
      all_genes_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    plot_data_list[[analysis_index]] <- prepare_group_data(
      dat = dat,
      analysis_name = analysis_name
    )
  }

  plot_data <- do.call(rbind, plot_data_list)
  plot_data <- plot_data[order(plot_data$Regulation), , drop = FALSE]

  plot_data_for_axis <- plot_data[
    plot_data$Regulation %in% c("Up", "Down"),
    ,
    drop = FALSE
  ]
  if (nrow(plot_data_for_axis) == 0) {
    stop("No significant Up/Down genes were found for scheme: ", scheme_name)
  }

  group_layout <- get_group_layout_data(
    plot_data = plot_data,
    selected_analyses = selected_analyses
  )
  plot_data <- apply_group_label_display_y(
    plot_data = plot_data,
    group_layout = group_layout
  )

  axis_info <- get_axis_info(
    plot_data = plot_data,
    group_layout = group_layout
  )
  pdf_size <- get_pdf_size(
    group_count = length(selected_analyses),
    display_y_limits = axis_info$display_limits,
    group_layout = group_layout
  )

  multiple_volcano_plot <- make_multiple_volcano_plot(
    plot_data = plot_data,
    group_layout = group_layout,
    selected_analyses = selected_analyses,
    axis_info = axis_info
  )

  output_dir <- file.path(PLOT_ROOT, sanitize_file_name(scheme_name))
  if (OVERWRITE_SCHEME_OUTPUT && dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pdf_file <- file.path(output_dir, "multiple_volcano_plot.pdf")

  Cairo::CairoPDF(
    file = pdf_file,
    width = pdf_size$width,
    height = pdf_size$height,
    bg = "white"
  )
  print(multiple_volcano_plot)
  invisible(dev.off())

  stopifnot(file.exists(pdf_file))

  summary_table <- do.call(rbind, lapply(selected_analyses, function(analysis_name) {
    dat <- plot_data[plot_data$Analysis_Name == analysis_name, , drop = FALSE]
    plotted_dat <- dat[dat$Regulation %in% c("Up", "Down"), , drop = FALSE]
    status_counts <- table(dat$Regulation)
    count_status <- function(status_name) {
      value <- status_counts[status_name]
      if (is.na(value)) {
        return(0L)
      }

      as.integer(value)
    }

    layout_row <- group_layout[
      group_layout$Analysis_Name == analysis_name,
      ,
      drop = FALSE
    ]

    data.frame(
      Plot_Name = scheme_name,
      Analysis_Name = analysis_name,
      Total_Genes = nrow(dat),
      Genes_Plotted = nrow(plotted_dat),
      Up = count_status("Up"),
      Down = count_status("Down"),
      NS = count_status("NS"),
      X_Min = round(layout_row$X_Min, 2),
      X_Max = round(layout_row$X_Max, 2),
      Y_Min = round(axis_info$true_limits[1], 2),
      Y_Max = round(axis_info$true_limits[2], 2),
      PDF_Width = round(pdf_size$width, 2),
      PDF_Height = round(pdf_size$height, 2),
      PDF_File = pdf_file,
      stringsAsFactors = FALSE
    )
  }))

  summary_table
}


# 3. 同步阈值并查找输入文件 ----------------------------------------------------

if (SYNC_THRESHOLDS_FROM_01_SCRIPT) {
  P_VALUE_COLUMN <- read_scalar_config(
    DE_SCRIPT_FILE,
    "P_VALUE_COLUMN",
    P_VALUE_COLUMN
  )
  P_VALUE_CUTOFF <- read_scalar_config(
    DE_SCRIPT_FILE,
    "P_VALUE_CUTOFF",
    P_VALUE_CUTOFF
  )
  LOGFC_CUTOFF <- read_scalar_config(
    DE_SCRIPT_FILE,
    "LOGFC_CUTOFF",
    LOGFC_CUTOFF
  )
}

stopifnot(is.character(P_VALUE_COLUMN))
stopifnot(is.numeric(P_VALUE_CUTOFF))
stopifnot(is.numeric(LOGFC_CUTOFF))
stopifnot(P_VALUE_CUTOFF > 0)
stopifnot(LOGFC_CUTOFF > 0)

file_info <- get_deg_file_info(TABLE_ROOT)

stopifnot(is.list(MULTIPLE_VOLCANO_SCHEMES))
stopifnot(length(MULTIPLE_VOLCANO_SCHEMES) > 0)
stopifnot(is.character(SCHEMES_TO_RUN))
stopifnot(length(SCHEMES_TO_RUN) > 0)

missing_schemes <- setdiff(SCHEMES_TO_RUN, names(MULTIPLE_VOLCANO_SCHEMES))
if (length(missing_schemes) > 0) {
  stop(
    "The following schemes are not defined in MULTIPLE_VOLCANO_SCHEMES: ",
    paste(missing_schemes, collapse = ", ")
  )
}


# 4. 绘制并保存多组火山图 ------------------------------------------------------

cat("\nRunning multiple volcano plot generation...\n")
cat("Schemes: ", paste(SCHEMES_TO_RUN, collapse = ", "), "\n", sep = "")
cat("P value column: ", P_VALUE_COLUMN, "\n", sep = "")
cat("P value cutoff: ", P_VALUE_CUTOFF, "\n", sep = "")
cat("logFC cutoff: ", LOGFC_CUTOFF, "\n", sep = "")

if (CLEAN_MULTIPLE_VOLCANO_ROOT && dir.exists(PLOT_ROOT)) {
  unlink(PLOT_ROOT, recursive = TRUE)
}
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)

summary_list <- vector("list", length(SCHEMES_TO_RUN))
names(summary_list) <- SCHEMES_TO_RUN

for (scheme_index in seq_along(SCHEMES_TO_RUN)) {
  scheme_name <- SCHEMES_TO_RUN[scheme_index]
  selected_analyses <- MULTIPLE_VOLCANO_SCHEMES[[scheme_name]]

  cat(
    "[",
    scheme_index,
    "/",
    length(SCHEMES_TO_RUN),
    "] ",
    scheme_name,
    ": ",
    paste(selected_analyses, collapse = ", "),
    "\n",
    sep = ""
  )

  summary_list[[scheme_index]] <- run_multiple_volcano_scheme(
    scheme_name = scheme_name,
    selected_analyses = selected_analyses,
    file_info = file_info
  )
}

summary_table <- do.call(rbind, summary_list)


# 5. 终端快速汇总 --------------------------------------------------------------

cat("\nMultiple volcano plot summary:\n")
print(
  summary_table[
    ,
    c(
      "Plot_Name", "Analysis_Name", "Total_Genes", "Genes_Plotted",
      "Up", "Down", "NS", "X_Min", "X_Max", "Y_Min", "Y_Max",
      "PDF_Width", "PDF_Height"
    )
  ],
  row.names = FALSE
)

cat("\nMultiple volcano plot was saved in:\n")
saved_files <- unique(summary_table$PDF_File)
cat(paste(saved_files, collapse = "\n"), "\n", sep = "")
cat("\nMultiple volcano plot generation finished.\n")
