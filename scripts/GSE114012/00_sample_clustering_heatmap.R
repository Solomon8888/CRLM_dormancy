# GSE114012 LRC sample clustering heatmap
#
# Extract all LRC samples from the GSE114012 SummarizedExperiment object and
# draw a sample-sample TPM correlation heatmap with hierarchical clustering.
# Sample labels are taken from the clinical Title column.


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"

SE_RDS_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_clinical_edit.csv"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"

PLOT_DIR <- "results/ngs/GSE114012/plots"

# 可选："coding", "protein", "protein_coding", "non_coding", "all"
GENE_BIOTYPE_FILTER <- "coding"

TPM_ASSAY_NAME <- "tpm"

LRC_COLUMN <- "dormant_lrc_or_cycling_bulk_"
LRC_VALUE <- "LRC"
SAMPLE_LABEL_COLUMN <- "Title"

CORRELATION_METHOD <- "pearson"
CLUSTERING_METHOD <- "complete"

# 样本名较长时换行。列名在45度显示，默认更保守换行以避免重叠。
ROW_LABEL_WIDTH <- 45
COLUMN_LABEL_WIDTH <- 45

# 热图主体中每个小格子的边长。行列样本数一致时，热图主体保持正方形。
CELL_SIZE_MM <- 12

# 聚类树和样本注释条的空间。聚类树只反映TPM表达相关性，不按分组强行排序。
ROW_DEND_WIDTH_MM <- 28
COLUMN_DEND_HEIGHT_MM <- 24
SAMPLE_ANNOTATION_SIZE_MM <- 4
DENDROGRAM_LINE_WIDTH_SCALE <- 0.95

# 加粗倍率。当前按4倍设置，用于热图框线和样本名显示。
BOLDNESS_MULTIPLIER <- 4
BASE_CELL_BORDER_WIDTH <- 0.95
BASE_SAMPLE_FONT_SIZE <- 3.8
BASE_ANNOTATION_FONT_SIZE <- 4.2
BASE_LEGEND_FONT_SIZE <- 4.2

CELL_BORDER_WIDTH <- BASE_CELL_BORDER_WIDTH * BOLDNESS_MULTIPLIER
DENDROGRAM_LINE_WIDTH <- CELL_BORDER_WIDTH * DENDROGRAM_LINE_WIDTH_SCALE

# 文字统一使用黑色粗体。样本名字号不变，只统一字体和加粗方式。
TEXT_FONT_FAMILY <- "Helvetica"
TEXT_FONT_FACE <- "bold"
TEXT_COLOR <- "black"
LABEL_HEATMAP_GAP_SPACES <- 2

# 红白蓝渐变配色。需要换色时只改这里即可。
HEATMAP_COLOR_LOW <- "#0d0dbb7f"   # deep blue
HEATMAP_COLOR_MID <- "#FFFFFF"   # white
HEATMAP_COLOR_HIGH <- "#cd0e0e"  # bright red
CORRELATION_COLOR_MIN <- 0.80
CORRELATION_COLOR_MAX <- 1.00

# PDF大小会根据样本数量和样本名长度自动调整。
MIN_PDF_WIDTH <- 8.0
MIN_PDF_HEIGHT <- 8.0
MAX_PDF_WIDTH <- 24.0
MAX_PDF_HEIGHT <- 22.0

SAMPLE_FONT_SIZE <- BASE_SAMPLE_FONT_SIZE * BOLDNESS_MULTIPLIER
COLUMN_FONT_SCALE <- 1.00
COLUMN_SAMPLE_FONT_SIZE <- SAMPLE_FONT_SIZE * COLUMN_FONT_SCALE
ANNOTATION_FONT_SIZE <- BASE_ANNOTATION_FONT_SIZE * BOLDNESS_MULTIPLIER
LEGEND_FONT_SIZE <- BASE_LEGEND_FONT_SIZE * BOLDNESS_MULTIPLIER

# 图例放在整体图片左侧，并与热图主体保持较大间距。
OUTER_MARGIN_INCH <- 0.35
LEGEND_LEFT_WIDTH_INCH <- 1.65
LEGEND_HEATMAP_GAP_INCH <- 0.45
LEGEND_INNER_MARGIN_INCH <- 0.05
LEGEND_TOP_EXTRA_MM <- 26
LEGEND_ITEM_GAP_MM <- 5
LEGEND_GROUP_GAP_MM <- 18
CORRELATION_LEGEND_HEIGHT_MM <- 50
LEGEND_TOP_MARGIN_INCH <- (
  COLUMN_DEND_HEIGHT_MM + SAMPLE_ANNOTATION_SIZE_MM + LEGEND_TOP_EXTRA_MM
) / 25.4

options(width = 200)


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
  library(Cairo)
})

source(FUNCTION_FILE)


# 2. 读取数据 ------------------------------------------------------------------

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot(SAMPLE_LABEL_COLUMN %in% colnames(clinical_data))
stopifnot(LRC_COLUMN %in% colnames(clinical_data))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

missing_samples <- setdiff(colnames(se), clinical_data$Sample_ID)
stopifnot(length(missing_samples) == 0)

sample_info_all <- clinical_data[
  match(colnames(se), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
rownames(sample_info_all) <- sample_info_all$Sample_ID
stopifnot(all(sample_info_all$Sample_ID == colnames(se)))


# 3. 提取LRC样本 ---------------------------------------------------------------

lrc_status <- trimws(as.character(sample_info_all[[LRC_COLUMN]]))
lrc_status[is.na(lrc_status)] <- ""
lrc_sample_index <- lrc_status == LRC_VALUE

sample_info <- sample_info_all[lrc_sample_index, , drop = FALSE]
stopifnot(nrow(sample_info) >= 2)

sample_labels <- trimws(as.character(sample_info[[SAMPLE_LABEL_COLUMN]]))
sample_labels[sample_labels == "" | is.na(sample_labels)] <- sample_info$Sample_ID[
  sample_labels == "" | is.na(sample_labels)
]

stopifnot(!any(duplicated(sample_labels)))


# 4. 准备TPM表达矩阵 -----------------------------------------------------------

assay_names <- names(SummarizedExperiment::assays(se))
stopifnot(TPM_ASSAY_NAME %in% assay_names)

tpm_all <- as.matrix(SummarizedExperiment::assay(se, TPM_ASSAY_NAME))
tpm_all <- tpm_all[, sample_info$Sample_ID, drop = FALSE]
stopifnot(is.numeric(tpm_all))

feature_id <- rownames(tpm_all)
if (is.null(feature_id)) {
  feature_id <- paste0("Feature_", seq_len(nrow(tpm_all)))
  rownames(tpm_all) <- feature_id
}

gene_annotation <- data.frame(
  Feature_ID = feature_id,
  as.data.frame(rowData(se), stringsAsFactors = FALSE),
  check.names = FALSE
)
rownames(gene_annotation) <- rownames(tpm_all)

gene_biotype_filter <- trimws(tolower(GENE_BIOTYPE_FILTER))
gene_biotype_filter <- gsub("-", "_", gene_biotype_filter)
if (gene_biotype_filter %in% c("protein", "protein_coding")) {
  gene_biotype_filter <- "coding"
}

gene_filter <- filter_genes_by_biotype(
  exprSet = tpm_all,
  gene_annotation = gene_annotation,
  biotype_filter = gene_biotype_filter
)

tpm <- gene_filter$exprSet
expr_for_correlation <- log2(tpm + 1)

gene_sd <- apply(expr_for_correlation, 1, sd, na.rm = TRUE)
expr_for_correlation <- expr_for_correlation[
  is.finite(gene_sd) & gene_sd > 0,
  ,
  drop = FALSE
]
stopifnot(nrow(expr_for_correlation) > 1)


# 5. 计算样本相关性和层次聚类 --------------------------------------------------

cor_matrix <- cor(
  expr_for_correlation,
  method = CORRELATION_METHOD,
  use = "pairwise.complete.obs"
)

stopifnot(!any(is.na(cor_matrix)))

row_distance <- as.dist(1 - cor_matrix)
sample_hclust <- hclust(row_distance, method = CLUSTERING_METHOD)

# 聚类只基于TPM表达谱得到的样本相关性，不按细胞系分组强行排序。


# 6. 样本标签、注释和PDF尺寸 ---------------------------------------------------

wrap_label <- function(x, width = 45) {
  x <- as.character(x)

  vapply(x, function(label) {
    n <- nchar(label)
    if (n <= width) return(label)

    starts <- seq(1, n, by = width)
    parts <- substring(label, starts, pmin(starts + width - 1, n))
    paste(parts, collapse = "\n")
  }, character(1))
}

wrap_label_by_underscore <- function(x, width = 8) {
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

    lines <- c(lines, current_line)
    paste(lines, collapse = "\n")
  }, character(1))
}

plot_row_labels <- wrap_label(sample_labels, ROW_LABEL_WIDTH)
plot_column_labels <- wrap_label_by_underscore(sample_labels, COLUMN_LABEL_WIDTH)

label_padding <- paste(rep(" ", LABEL_HEATMAP_GAP_SPACES), collapse = "")
plot_row_labels <- paste0(label_padding, plot_row_labels)
plot_column_labels <- paste0(label_padding, plot_column_labels)

rownames(cor_matrix) <- sample_info$Sample_ID
colnames(cor_matrix) <- sample_info$Sample_ID

cell_line_levels <- sort(unique(sample_info$cell_line))
cell_line_palette <- RColorBrewer::brewer.pal(
  max(3, length(cell_line_levels)),
  "Set2"
)[seq_along(cell_line_levels)]
names(cell_line_palette) <- cell_line_levels

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Cell_line = sample_info$cell_line,
  col = list(Cell_line = cell_line_palette),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = grid::gpar(
    fontsize = ANNOTATION_FONT_SIZE,
    fontface = TEXT_FONT_FACE,
    fontfamily = TEXT_FONT_FAMILY,
    col = TEXT_COLOR
  ),
  simple_anno_size = grid::unit(SAMPLE_ANNOTATION_SIZE_MM, "mm"),
  gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
  border = TRUE
)

left_annotation <- ComplexHeatmap::rowAnnotation(
  Cell_line = sample_info$cell_line,
  col = list(Cell_line = cell_line_palette),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = grid::unit(SAMPLE_ANNOTATION_SIZE_MM, "mm"),
  gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
  border = TRUE
)

n_samples <- ncol(cor_matrix)
row_label_parts <- unlist(strsplit(plot_row_labels, "\n", fixed = TRUE))
column_label_parts <- strsplit(plot_column_labels, "\n", fixed = TRUE)

max_row_label_chars <- max(nchar(row_label_parts))
max_column_label_chars <- max(nchar(unlist(column_label_parts)))
max_column_label_lines <- max(lengths(column_label_parts))

heatmap_body_inch <- n_samples * CELL_SIZE_MM / 25.4
row_label_space <- max_row_label_chars * SAMPLE_FONT_SIZE * 0.010
col_label_space <- max_column_label_chars * COLUMN_SAMPLE_FONT_SIZE * 0.010 +
  max_column_label_lines * COLUMN_SAMPLE_FONT_SIZE * 0.020

heatmap_panel_width <- heatmap_body_inch + row_label_space +
  ROW_DEND_WIDTH_MM / 25.4 + SAMPLE_ANNOTATION_SIZE_MM / 25.4 + 0.95
heatmap_panel_height <- heatmap_body_inch + col_label_space +
  COLUMN_DEND_HEIGHT_MM / 25.4 + SAMPLE_ANNOTATION_SIZE_MM / 25.4 + 0.95

pdf_width <- OUTER_MARGIN_INCH * 2 +
  LEGEND_LEFT_WIDTH_INCH + LEGEND_HEATMAP_GAP_INCH + heatmap_panel_width
pdf_height <- OUTER_MARGIN_INCH * 2 + heatmap_panel_height

pdf_width <- min(max(pdf_width, MIN_PDF_WIDTH), MAX_PDF_WIDTH)
pdf_height <- min(max(pdf_height, MIN_PDF_HEIGHT), MAX_PDF_HEIGHT)

plot_filter_name <- sanitize_file_name(gene_filter$filter)
pdf_file <- file.path(
  PLOT_DIR,
  paste0(DATASET_ID, "_LRC_", plot_filter_name, "_TPM_sample_clustering_heatmap.pdf")
)


# 7. 绘制聚类热图 --------------------------------------------------------------

color_min <- CORRELATION_COLOR_MIN
color_max <- CORRELATION_COLOR_MAX
color_mid <- (color_min + color_max) / 2

heatmap_colors <- circlize::colorRamp2(
  c(color_min, color_mid, color_max),
  c(HEATMAP_COLOR_LOW, HEATMAP_COLOR_MID, HEATMAP_COLOR_HIGH)
)

heatmap_body_size <- grid::unit(n_samples * CELL_SIZE_MM, "mm")

ht <- ComplexHeatmap::Heatmap(
  cor_matrix,
  name = "Correlation",
  col = heatmap_colors,
  cluster_rows = sample_hclust,
  cluster_columns = sample_hclust,
  top_annotation = top_annotation,
  left_annotation = left_annotation,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_labels = plot_row_labels,
  column_labels = plot_column_labels,
  row_names_side = "right",
  column_names_side = "bottom",
  row_names_max_width = grid::unit(row_label_space, "in"),
  column_names_max_height = grid::unit(col_label_space, "in"),
  row_names_gp = grid::gpar(
    fontsize = SAMPLE_FONT_SIZE,
    fontface = TEXT_FONT_FACE,
    fontfamily = TEXT_FONT_FAMILY,
    col = TEXT_COLOR
  ),
  column_names_gp = grid::gpar(
    fontsize = COLUMN_SAMPLE_FONT_SIZE,
    fontface = TEXT_FONT_FACE,
    fontfamily = TEXT_FONT_FAMILY,
    col = TEXT_COLOR
  ),
  column_names_rot = 45,
  rect_gp = grid::gpar(
    col = "black",
    lwd = CELL_BORDER_WIDTH
  ),
  border = TRUE,
  border_gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
  row_dend_width = grid::unit(ROW_DEND_WIDTH_MM, "mm"),
  column_dend_height = grid::unit(COLUMN_DEND_HEIGHT_MM, "mm"),
  row_dend_gp = grid::gpar(col = "black", lwd = DENDROGRAM_LINE_WIDTH),
  column_dend_gp = grid::gpar(col = "black", lwd = DENDROGRAM_LINE_WIDTH),
  width = heatmap_body_size,
  height = heatmap_body_size,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    title = "Correlation",
    at = round(seq(color_min, color_max, length.out = 5), 2),
    labels_gp = grid::gpar(
      fontsize = LEGEND_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    title_gp = grid::gpar(
      fontsize = LEGEND_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    )
  )
)

cell_line_legend <- ComplexHeatmap::Legend(
  title = "Cell_line",
  at = cell_line_levels,
  type = "grid",
  legend_gp = grid::gpar(
    fill = cell_line_palette[cell_line_levels],
    col = "black",
    lwd = CELL_BORDER_WIDTH
  ),
  labels_gp = grid::gpar(
    fontsize = LEGEND_FONT_SIZE,
    fontface = TEXT_FONT_FACE,
    fontfamily = TEXT_FONT_FAMILY,
    col = TEXT_COLOR
  ),
  title_gp = grid::gpar(
    fontsize = LEGEND_FONT_SIZE,
    fontface = TEXT_FONT_FACE,
    fontfamily = TEXT_FONT_FAMILY,
    col = TEXT_COLOR
  ),
  grid_height = grid::unit(5, "mm"),
  grid_width = grid::unit(5, "mm"),
  gap = grid::unit(LEGEND_ITEM_GAP_MM, "mm"),
  row_gap = grid::unit(LEGEND_ITEM_GAP_MM, "mm"),
  title_gap = grid::unit(LEGEND_ITEM_GAP_MM, "mm")
)

correlation_legend <- ComplexHeatmap::Legend(
  title = "Correlation",
  col_fun = heatmap_colors,
  at = round(seq(color_min, color_max, length.out = 5), 2),
  labels_gp = grid::gpar(
    fontsize = LEGEND_FONT_SIZE,
    fontface = TEXT_FONT_FACE,
    fontfamily = TEXT_FONT_FAMILY,
    col = TEXT_COLOR
  ),
  title_gp = grid::gpar(
    fontsize = LEGEND_FONT_SIZE,
    fontface = TEXT_FONT_FACE,
    fontfamily = TEXT_FONT_FAMILY,
    col = TEXT_COLOR
  ),
  legend_height = grid::unit(CORRELATION_LEGEND_HEIGHT_MM, "mm"),
  grid_width = grid::unit(6, "mm"),
  title_gap = grid::unit(LEGEND_ITEM_GAP_MM, "mm")
)

legend_pack <- ComplexHeatmap::packLegend(
  cell_line_legend,
  correlation_legend,
  direction = "vertical",
  gap = grid::unit(LEGEND_GROUP_GAP_MM, "mm")
)

Cairo::CairoPDF(
  file = pdf_file,
  width = pdf_width,
  height = pdf_height,
  bg = "white"
)

grid::grid.newpage()
grid::pushViewport(grid::viewport(
  layout = grid::grid.layout(
    nrow = 3,
    ncol = 5,
    widths = grid::unit.c(
      grid::unit(OUTER_MARGIN_INCH, "in"),
      grid::unit(LEGEND_LEFT_WIDTH_INCH, "in"),
      grid::unit(LEGEND_HEATMAP_GAP_INCH, "in"),
      grid::unit(heatmap_panel_width, "in"),
      grid::unit(OUTER_MARGIN_INCH, "in")
    ),
    heights = grid::unit.c(
      grid::unit(OUTER_MARGIN_INCH, "in"),
      grid::unit(heatmap_panel_height, "in"),
      grid::unit(OUTER_MARGIN_INCH, "in")
    )
  )
))

grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
ComplexHeatmap::draw(
  legend_pack,
  x = grid::unit(LEGEND_INNER_MARGIN_INCH, "in"),
  y = grid::unit(1, "npc") - grid::unit(LEGEND_TOP_MARGIN_INCH, "in"),
  just = c("left", "top")
)
grid::popViewport()

grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 4))
ComplexHeatmap::draw(
  ht,
  newpage = FALSE,
  show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)
grid::popViewport(2)

invisible(dev.off())


# 8. 输出信息 ------------------------------------------------------------------

cat("\nLRC sample TPM clustering heatmap finished.\n")
cat("LRC samples: ", nrow(sample_info), "\n", sep = "")
cat("Gene biotype filter: ", gene_filter$filter, "\n", sep = "")
cat("TPM assay: ", TPM_ASSAY_NAME, "\n", sep = "")
cat("Genes used after biotype filter: ", gene_filter$selected_gene_count, "\n", sep = "")
cat("Genes used after removing zero-variance genes: ", nrow(expr_for_correlation), "\n", sep = "")
cat("PDF size: ", round(pdf_width, 2), " x ", round(pdf_height, 2), " inches\n", sep = "")
cat("PDF file: ", pdf_file, "\n", sep = "")
