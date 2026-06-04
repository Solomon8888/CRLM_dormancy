# GSE114012多组差异表达火山图
#
# 读取01号limma脚本生成的all_genes.csv结果，按配置将多个差异分析设计
# 合并到一张多组火山图中。每个差异分析设计占一个分面；
# 横轴为-log10(P值)，纵轴为logFC。图中只展示达到阈值的Up和Down基因；
# NS基因只参与终端汇总，不在图中绘制。阈值与01号脚本保持一致。


# 0. 可修改配置 ---------------------------------------------------------------

# 当前脚本只服务于GSE114012这个NGS数据集；目录结构也依赖这两个字段。
DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

# 01号脚本用于同步P值列名和显著性阈值；
# PLOTTING_FUNCTION_FILE保存跨绘图脚本共用的风格配置和基础函数。
DE_SCRIPT_FILE <- "scripts/GSE114012/01_limma_differential_expression.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"

# 输入目录为01号脚本输出的tables/<analysis_name>/DEG/all_genes.csv。
# 输出目录为plots/multiple_volcano/<scheme_name>/。
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
# 每个方案会单独输出PDF和PNG：
# results/ngs/GSE114012/plots/multiple_volcano/<scheme_name>/multiple_volcano_plot.*
MULTIPLE_VOLCANO_SCHEMES <- list(
  ALL = c("DLD1", "HCT15", "HT55", "RKO", "SW48", "SW948"),
  DLD1_HCT15_SW48 = c("DLD1", "HCT15", "SW48"),
  DLD1_HCT15 = c("DLD1", "HCT15"),
  SW48_RKO = c("SW48", "RKO"),
  HT55_SW948 = c("HT55", "SW948"),
  COMBINATION_1 = c("DLD1_HCT15_SW48", "DLD1_HCT15")
)

# 运行哪些多组火山图方案。
# 可设为names(MULTIPLE_VOLCANO_SCHEMES)运行全部；也可只写部分方案名。
SCHEMES_TO_RUN <- names(MULTIPLE_VOLCANO_SCHEMES)
# SCHEMES_TO_RUN <- c("CRC_LRC_core")

# 重跑时清理整个multiple_volcano目录，保证最终只保留当前配置生成的方案图。
CLEAN_MULTIPLE_VOLCANO_ROOT <- TRUE

# 重跑时清理当前方案目录内旧图片，避免旧文件残留。
OVERWRITE_SCHEME_OUTPUT <- TRUE

# 每组标注Top显著基因。Top排序优先按P值，其次按logFC幅度。
# 若CUSTOM_LABEL_GENES为空，则自动标注每组Up 5个和Down 5个Top基因。
# 若CUSTOM_LABEL_GENES不为空，则只标注这里配置的基因，不再自动标注Top基因。
# 用法1：所有分析都标注同一批基因。
# CUSTOM_LABEL_GENES <- c("MYH15", "CSF1R", "SRSF1")
# 用法2：不同分析标注不同基因；all表示所有分析都标注。
# CUSTOM_LABEL_GENES <- list(all = c("MYH15"), DLD1 = c("ABCG1", "SRSF1"))
CUSTOM_LABEL_GENES <- character(0)

# 多组火山图中Top基因上下避让距离。其他Top基因标注风格来自公共配置文件。
TOP_GENE_LABEL_NUDGE_Y <- 0.80

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
GROUP_LABEL_TEXT_DARKEN <- 0.78
GROUP_LABEL_BORDER_DARKEN <- 0.88

# 组名框位于logFC阈值内侧：当前阈值0.5时为±0.4；阈值1时为±0.9。
GROUP_LABEL_BOX_LOGFC_GAP <- 0.10
GROUP_LABEL_BOX_MIN_FRACTION <- 0.80

# 多行组名时，字体和行距会轻微收紧，保证文字被彩色外框完整覆盖。
GROUP_LABEL_FONT_MIN_SIZE <- 3.9
GROUP_LABEL_FONT_LINE_SHRINK <- 0.55
GROUP_LABEL_LINE_HEIGHT_MIN <- 0.95
GROUP_LABEL_LINE_HEIGHT_SHRINK <- 0.08

# 彩色组名框横向宽度略小于单个分面宽度，并保持各分面一致。
GROUP_LABEL_BOX_X_MARGIN_FRACTION <- 0.05
GROUP_LABEL_BORDER_WIDTH <- 0.9

# 分面边框、分面间距、图例位置和坐标留白。
PANEL_SPACING_X_MM <- 4.6
LEGEND_TOP_MARGIN_PT <- 68
X_AXIS_PADDING_FRACTION <- 0.12
X_AXIS_PADDING_MIN <- 0.22
X_AXIS_LEFT_MIN_FRACTION <- 0.78
Y_AXIS_PADDING_FRACTION <- 0.07

# 组名框越高，PDF纵向空间会等比例略微增加，避免压缩真实logFC坐标区域。
LABEL_BOX_HEIGHT_PDF_SCALE <- 1.8

# PDF/PNG和字体设置。文件名固定为multiple_volcano_plot.*；方案名由目录体现。
BASE_PDF_HEIGHT <- 6.0
GROUP_WIDTH_INCH <- 1.78
LEGEND_WIDTH_INCH <- 1.18
MIN_PDF_WIDTH <- 7.2
MAX_PDF_WIDTH <- 20.0
MAX_PDF_HEIGHT <- 8.2

# 右侧图例的文字和圆点大小。
LEGEND_TEXT_SIZE <- 13.5
LEGEND_POINT_SIZE_SCALE <- 1.45
LEGEND_KEY_SIZE_MM <- 6.6

options(width = 200)


# 1. 加载包 -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

source(PLOTTING_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

USE_GG_REPEL <- requireNamespace("ggrepel", quietly = TRUE)


# 2. 常用函数 -----------------------------------------------------------------

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

  top_label_data <- get_volcano_label_data(
    plot_data = plot_data,
    custom_label_genes = CUSTOM_LABEL_GENES,
    symbol_column = TOP_GENE_SYMBOL_COLUMN,
    match_columns = CUSTOM_LABEL_MATCH_COLUMNS,
    p_value_column = P_VALUE_COLUMN,
    top_up_n = TOP_UP_LABEL_N,
    top_down_n = TOP_DOWN_LABEL_N
  )
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
  top_label_colors <- get_regulation_label_colors(
    label_data = top_label_data,
    up_color = UP_COLOR,
    down_color = DOWN_COLOR,
    darken_fraction = TOP_GENE_LABEL_COLOR_DARKEN
  )

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
  group_label_text_colors <- darken_color(
    group_label_colors,
    fraction = GROUP_LABEL_TEXT_DARKEN
  )
  group_label_border_colors <- darken_color(
    group_label_colors,
    fraction = GROUP_LABEL_BORDER_DARKEN
  )

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
        color = group_label_border_colors[analysis_name],
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
        color = group_label_text_colors[analysis_name]
      )
  }

  add_volcano_gene_label_layer(
    plot = volcano_plot,
    label_data = top_label_data,
    label_colors = top_label_colors,
    use_gg_repel = USE_GG_REPEL,
    text_family = TEXT_FONT_FAMILY,
    fontface = TOP_GENE_LABEL_FONT_FACE,
    font_size = TOP_GENE_LABEL_FONT_SIZE,
    box_padding = TOP_GENE_LABEL_BOX_PADDING,
    point_padding = TOP_GENE_LABEL_POINT_PADDING,
    segment_width = TOP_GENE_LABEL_SEGMENT_WIDTH,
    force = TOP_GENE_LABEL_FORCE,
    force_pull = TOP_GENE_LABEL_FORCE_PULL,
    max_overlaps = TOP_GENE_LABEL_MAX_OVERLAPS,
    nudge_y = top_label_data$Label_Nudge_Y,
    fallback_vjust = ifelse(top_label_data$Regulation == "Up", -0.7, 1.2)
  )
}

run_multiple_volcano_scheme <- function(scheme_name, selected_analyses, file_info) {
  # 单套方案的完整流程：读取数据、绘图、保存图片、返回终端汇总表。
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

  plot_data_list <- vector("list", length(selected_analyses))
  names(plot_data_list) <- selected_analyses

  for (analysis_index in seq_along(selected_analyses)) {
    analysis_name <- selected_analyses[analysis_index]
    dat <- read_deg_result(file_info, analysis_name)

    plot_data_list[[analysis_index]] <- prepare_volcano_data(
      dat = dat,
      analysis_name = analysis_name,
      p_value_column = P_VALUE_COLUMN,
      p_value_cutoff = P_VALUE_CUTOFF,
      logfc_cutoff = LOGFC_CUTOFF,
      ns_label = "NS",
      regulation_levels = c("NS", "Down", "Up")
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

  output_files <- save_ggplot_pdf_png(
    plot = multiple_volcano_plot,
    pdf_file = pdf_file,
    width = pdf_size$width,
    height = pdf_size$height
  )

  summary_table <- do.call(rbind, lapply(selected_analyses, function(analysis_name) {
    dat <- plot_data[plot_data$Analysis_Name == analysis_name, , drop = FALSE]
    plotted_dat <- dat[dat$Regulation %in% c("Up", "Down"), , drop = FALSE]
    status_counts <- table(dat$Regulation)

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
      Up = count_status(status_counts, "Up"),
      Down = count_status(status_counts, "Down"),
      NS = count_status(status_counts, "NS"),
      X_Min = round(layout_row$X_Min, 2),
      X_Max = round(layout_row$X_Max, 2),
      Y_Min = round(axis_info$true_limits[1], 2),
      Y_Max = round(axis_info$true_limits[2], 2),
      PDF_Width = round(pdf_size$width, 2),
      PDF_Height = round(pdf_size$height, 2),
      PDF_File = output_files$pdf_file,
      PNG_File = output_files$png_file,
      stringsAsFactors = FALSE
    )
  }))

  summary_table
}


# 3. 同步阈值并查找输入文件 ----------------------------------------------------

threshold_config <- sync_de_thresholds_from_script(
  script_file = DE_SCRIPT_FILE,
  p_value_column = P_VALUE_COLUMN,
  p_value_cutoff = P_VALUE_CUTOFF,
  logfc_cutoff = LOGFC_CUTOFF,
  sync = SYNC_THRESHOLDS_FROM_01_SCRIPT
)
P_VALUE_COLUMN <- threshold_config$p_value_column
P_VALUE_CUTOFF <- threshold_config$p_value_cutoff
LOGFC_CUTOFF <- threshold_config$logfc_cutoff

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

run_one_multiple_volcano_scheme <- function(scheme_index) {
  scheme_name <- SCHEMES_TO_RUN[scheme_index]
  selected_analyses <- MULTIPLE_VOLCANO_SCHEMES[[scheme_name]]

  run_multiple_volcano_scheme(
    scheme_name = scheme_name,
    selected_analyses = selected_analyses,
    file_info = file_info
  )
}

parallel_strategy <- setup_parallel_strategy(
  total_tasks = length(SCHEMES_TO_RUN),
  inner_label = "Multiple volcano inner workers",
  nested_label = "Nested workers"
)

summary_list <- run_indexed_tasks_with_progress(
  total_tasks = length(SCHEMES_TO_RUN),
  workers = parallel_strategy$task_workers,
  task_function = run_one_multiple_volcano_scheme
)
stop_on_parallel_errors(summary_list, task_ids = SCHEMES_TO_RUN, label = "multiple volcano schemes")
names(summary_list) <- SCHEMES_TO_RUN

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

cat("\nMultiple volcano PDF plots were saved in:\n")
saved_files <- unique(summary_table$PDF_File)
cat(paste(saved_files, collapse = "\n"), "\n", sep = "")
cat("\nMultiple volcano PNG plots were saved in:\n")
saved_png_files <- unique(summary_table$PNG_File)
cat(paste(saved_png_files, collapse = "\n"), "\n", sep = "")
cat("\nMultiple volcano plot generation finished.\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
