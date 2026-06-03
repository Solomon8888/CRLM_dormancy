# GSE114012差异表达火山图
#
# 读取01号limma脚本生成的all_genes.csv结果，为每个差异分析设计绘制火山图。
# 横坐标为logFC，纵坐标为-log10(P值)；上调基因标红，下调基因标蓝，
# 未达到阈值的基因标灰。输出PDF为矢量格式，适合后续排版和发表使用。


# 0. 可修改配置 ---------------------------------------------------------------

# 当前脚本只服务于GSE114012这个NGS数据集；目录结构也依赖这两个字段。
DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

# 01号脚本用于同步P值列名和显著性阈值；公共绘图函数脚本用于读取配置和查找DEG结果。
DE_SCRIPT_FILE <- "scripts/GSE114012/01_limma_differential_expression.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"

# 输入目录为01号脚本输出的tables/<analysis_name>/DEG/all_genes.csv。
# 输出目录为plots/volcano/<analysis_name>/volcano_plot.pdf。
RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "volcano")

# 是否从01号差异分析脚本中同步显著性阈值。
# 建议保持TRUE，保证火山图阈值与差异分析结果完全一致。
SYNC_THRESHOLDS_FROM_01_SCRIPT <- TRUE

# 只有在SYNC_THRESHOLDS_FROM_01_SCRIPT为FALSE时，下面三项才作为手动阈值使用。
P_VALUE_COLUMN <- "P.Value"
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

# 需要绘制哪些分析。设为"all"时自动绘制全部tables/<analysis_name>/DEG/all_genes.csv。
ANALYSES_TO_PLOT <- "all"
# ANALYSES_TO_PLOT <- c("DLD1", "HCT15")

# 重跑时清理当前分析火山图目录内旧PDF，避免旧文件残留。
CLEAN_VOLCANO_OUTPUT_DIR <- TRUE

# 火山图颜色和点样式。Up/Down颜色与04号多组火山图保持一致。
UP_COLOR <- "#D73027"
DOWN_COLOR <- "#2166AC"
NOT_SIGNIFICANT_COLOR <- "#B8B8B8"
POINT_SIZE <- 3.2
POINT_ALPHA <- 0.60

# 阈值线和坐标轴样式。阈值线只用于03号传统火山图；04号多组火山图不展示阈值虚线。
THRESHOLD_LINE_COLOR <- "#333333"
THRESHOLD_LINE_WIDTH <- 0.45
THRESHOLD_LINE_TYPE <- "dashed"
AXIS_LINE_WIDTH <- 0.8

# PDF和字体设置。文件名固定为volcano_plot.pdf；分析名由目录体现。
# 火山图本体保持正方形；整体宽度只因右侧图例略大于高度。
BASE_PDF_HEIGHT <- 6.2
MAX_EXTRA_PDF_HEIGHT <- 0.8
LEGEND_WIDTH_INCH <- 0.95
RIGHT_LEGEND_GAP_INCH <- 0.18
MAX_PDF_WIDTH_HEIGHT_RATIO <- 1.28
PANEL_HEIGHT_WIDTH_RATIO <- 1.0

# 全部文字统一为Helvetica黑色粗体，便于和00/04号绘图脚本保持一致。
TEXT_FONT_FAMILY <- "Helvetica"
BASE_FONT_SIZE <- 12
TEXT_FONT_FACE <- "bold"

options(width = 200)


# 1. 加载包 -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(Cairo)
})

source(PLOTTING_FUNCTION_FILE)

# 2. 绘图函数 -----------------------------------------------------------------

get_axis_limits <- function(plot_data) {
  # 横坐标左右对称，保证0点居中；纵坐标按当前分析的P值范围自动放大。
  x_abs <- max(abs(plot_data$logFC), LOGFC_CUTOFF, na.rm = TRUE)
  x_limit <- ceiling(x_abs * 1.08 * 2) / 2
  x_limit <- max(x_limit, 1)

  y_threshold <- -log10(P_VALUE_CUTOFF)
  y_limit <- ceiling(max(plot_data$Neg_Log10_P, y_threshold, na.rm = TRUE) * 1.08)
  y_limit <- max(y_limit, 2)

  list(
    x = c(-x_limit, x_limit),
    y = c(0, y_limit)
  )
}

get_pdf_size <- function(axis_limits) {
  # 纵向空间随-log10(P值)范围略微增加；横向只为右侧图例预留少量空间。
  y_span <- diff(axis_limits$y)
  extra_height <- min(max((y_span - 6) * 0.04, 0), MAX_EXTRA_PDF_HEIGHT)
  pdf_height <- BASE_PDF_HEIGHT + extra_height

  pdf_width <- pdf_height + LEGEND_WIDTH_INCH + RIGHT_LEGEND_GAP_INCH
  pdf_width <- min(pdf_width, pdf_height * MAX_PDF_WIDTH_HEIGHT_RATIO)
  pdf_width <- max(pdf_width, pdf_height * 1.08)

  list(
    width = pdf_width,
    height = pdf_height
  )
}

make_volcano_plot <- function(plot_data, axis_limits) {
  ggplot(
    plot_data,
    aes(x = logFC, y = Neg_Log10_P, color = Regulation)
  ) +
    geom_point(
      size = POINT_SIZE,
      alpha = POINT_ALPHA,
      shape = 16,
      stroke = 0
    ) +
    geom_vline(
      xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
      linewidth = THRESHOLD_LINE_WIDTH,
      linetype = THRESHOLD_LINE_TYPE,
      color = THRESHOLD_LINE_COLOR
    ) +
    geom_hline(
      yintercept = -log10(P_VALUE_CUTOFF),
      linewidth = THRESHOLD_LINE_WIDTH,
      linetype = THRESHOLD_LINE_TYPE,
      color = THRESHOLD_LINE_COLOR
    ) +
    scale_color_manual(
      values = c(
        "Not significant" = NOT_SIGNIFICANT_COLOR,
        "Down" = DOWN_COLOR,
        "Up" = UP_COLOR
      ),
      breaks = c("Up", "Down", "Not significant"),
      labels = c("Up", "Down", "NS")
    ) +
    scale_x_continuous(
      limits = axis_limits$x,
      breaks = pretty(axis_limits$x, n = 7),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = axis_limits$y,
      breaks = pretty(axis_limits$y, n = 6),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "log2 fold change",
      y = paste0("-log10(", P_VALUE_COLUMN, ")"),
      color = NULL
    ) +
    theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    theme(
      panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = AXIS_LINE_WIDTH),
      axis.line = element_line(color = "black", linewidth = AXIS_LINE_WIDTH),
      axis.text = element_text(color = "black", face = TEXT_FONT_FACE),
      axis.title = element_text(color = "black", face = TEXT_FONT_FACE),
      aspect.ratio = PANEL_HEIGHT_WIDTH_RATIO,
      legend.position = "right",
      legend.text = element_text(color = "black", face = TEXT_FONT_FACE),
      legend.key = element_blank(),
      legend.key.height = grid::unit(5.5, "mm"),
      legend.key.width = grid::unit(5.5, "mm"),
      legend.box.spacing = grid::unit(8, "pt"),
      legend.margin = margin(0, 0, 0, 4, unit = "pt"),
      strip.text = element_text(color = "black", face = TEXT_FONT_FACE),
      text = element_text(color = "black", face = TEXT_FONT_FACE),
      plot.margin = margin(10, 12, 10, 10, unit = "pt")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = POINT_SIZE * 1.15, alpha = 0.85)
      )
    )
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

if (identical(ANALYSES_TO_PLOT, "all")) {
  selected_analyses <- file_info$Analysis_Name
} else {
  selected_analyses <- ANALYSES_TO_PLOT
}

missing_analyses <- setdiff(selected_analyses, file_info$Analysis_Name)
if (length(missing_analyses) > 0) {
  stop(
    "No all_genes.csv file was found for: ",
    paste(missing_analyses, collapse = ", ")
  )
}

selected_file_info <- file_info[
  match(selected_analyses, file_info$Analysis_Name),
  ,
  drop = FALSE
]


# 4. 绘制并保存火山图 ----------------------------------------------------------

summary_list <- vector("list", length(selected_analyses))

cat("\nRunning volcano plot generation...\n")
cat("P value column: ", P_VALUE_COLUMN, "\n", sep = "")
cat("P value cutoff: ", P_VALUE_CUTOFF, "\n", sep = "")
cat("logFC cutoff: ", LOGFC_CUTOFF, "\n", sep = "")

progress_bar <- utils::txtProgressBar(
  min = 0,
  max = length(selected_analyses),
  style = 3
)

for (i in seq_along(selected_analyses)) {
  analysis_name <- selected_analyses[i]
  all_genes_file <- selected_file_info$All_Genes_File[
    selected_file_info$Analysis_Name == analysis_name
  ]

  dat <- read.csv(
    all_genes_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  plot_data <- prepare_volcano_data(
    dat = dat,
    p_value_column = P_VALUE_COLUMN,
    p_value_cutoff = P_VALUE_CUTOFF,
    logfc_cutoff = LOGFC_CUTOFF,
    ns_label = "Not significant",
    regulation_levels = c("Not significant", "Down", "Up")
  )
  axis_limits <- get_axis_limits(plot_data)
  pdf_size <- get_pdf_size(axis_limits)
  volcano_plot <- make_volcano_plot(plot_data, axis_limits)

  output_dir <- file.path(PLOT_ROOT, sanitize_file_name(analysis_name))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (CLEAN_VOLCANO_OUTPUT_DIR) {
    unlink(list.files(output_dir, pattern = "[.]pdf$", full.names = TRUE))
  }

  pdf_file <- file.path(output_dir, "volcano_plot.pdf")

  Cairo::CairoPDF(
    file = pdf_file,
    width = pdf_size$width,
    height = pdf_size$height,
    bg = "white"
  )
  print(volcano_plot)
  invisible(dev.off())

  stopifnot(file.exists(pdf_file))

  status_counts <- table(plot_data$Regulation)
  summary_list[[i]] <- data.frame(
    Analysis_Name = analysis_name,
    Genes_Plotted = nrow(plot_data),
    Up = count_status(status_counts, "Up"),
    Down = count_status(status_counts, "Down"),
    Not_Significant = count_status(status_counts, "Not significant"),
    X_Min = axis_limits$x[1],
    X_Max = axis_limits$x[2],
    Y_Max = axis_limits$y[2],
    PDF_Width = round(pdf_size$width, 2),
    PDF_Height = round(pdf_size$height, 2),
    PDF_File = pdf_file,
    stringsAsFactors = FALSE
  )

  utils::setTxtProgressBar(progress_bar, i)
}

close(progress_bar)


# 5. 终端快速汇总 --------------------------------------------------------------

summary_table <- do.call(rbind, summary_list)

cat("\nVolcano plot summary:\n")
print(
  summary_table[
    ,
    c(
      "Analysis_Name", "Genes_Plotted", "Up", "Down",
      "Not_Significant", "X_Min", "X_Max", "Y_Max",
      "PDF_Width", "PDF_Height"
    )
  ],
  row.names = FALSE
)

cat("\nVolcano plots were saved in:\n")
cat(file.path(PLOT_ROOT, "<analysis_name>", "volcano_plot.pdf"), "\n", sep = "")
cat("\nVolcano plot generation finished.\n")
