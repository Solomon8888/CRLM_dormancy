# GSE114012 LRC/BULK样本聚类热图
#
# 从当前数据集的SummarizedExperiment对象中分别提取全部LRC和BULK样本，
# 使用TPM表达量计算样本间相关性，并绘制层次聚类热图。
# 样本显示名称来自临床信息表中的Title列。


# 0. 可修改配置 ---------------------------------------------------------------

# 当前脚本只服务于GSE114012这个NGS数据集；目录结构也依赖这两个字段。
DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

# 输入文件：SE对象提供TPM矩阵和基因注释，临床表提供样本Title和LRC/BULK分组。
SE_RDS_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_clinical_edit.csv"

# FUNCTION_FILE用于基因类型筛选；
# PLOTTING_FUNCTION_FILE保存跨绘图脚本共用的风格配置和基础函数。
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"

# 输出根目录。最终PDF保存到plots/sample_clustering_heatmap/<gene_biotype>/<sample_group>/heatmap.pdf。
RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "sample_clustering_heatmap")

# 重跑时清理当前热图输出目录里的旧PDF，避免新旧命名混在一起。
CLEAN_PLOT_OUTPUT_DIR <- TRUE

# 基因类型筛选。可选："coding", "protein", "protein_coding", "non_coding", "all"。
# 这里沿用01号差异分析脚本的基因类型筛选逻辑。
GENE_BIOTYPE_FILTER <- "coding"

# SE对象中TPM矩阵的assay名称。聚类热图使用TPM，不使用原始count。
TPM_ASSAY_NAME <- "tpm"

# 样本筛选列和值。脚本会依次为SAMPLE_GROUPS中的每一组样本绘制热图。
LRC_COLUMN <- "dormant_lrc_or_cycling_bulk_"
LRC_VALUE <- "LRC"
BULK_VALUE <- "BULK"
SAMPLE_GROUPS <- c(
  LRC = LRC_VALUE,
  BULK = BULK_VALUE
)

# 图中展示的样本名来源。若Title为空，会自动回退到Sample_ID。
SAMPLE_LABEL_COLUMN <- "Title"

# 样本相关性和层次聚类方法。聚类只基于TPM表达模式，不按细胞系分组强行排序。
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

# 样本名与热图之间的空格，避免粗体样本名贴到热图本体。
LABEL_HEATMAP_GAP_SPACES <- 2

# 红白蓝渐变配色。需要换色时只改这里即可。
HEATMAP_COLOR_LOW <- "#0d0dbb7f"   # 深蓝
HEATMAP_COLOR_MID <- "#FFFFFF"     # 白色
HEATMAP_COLOR_HIGH <- "#cd0e0e"    # 鲜红
CORRELATION_COLOR_MIN <- 0.80
CORRELATION_COLOR_MAX <- 1.00

# PDF大小会根据样本数量和样本名长度自动调整；上下限用于避免图片过小或过大。
MIN_PDF_WIDTH <- 8.0
MIN_PDF_HEIGHT <- 8.0
MAX_PDF_WIDTH <- 24.0
MAX_PDF_HEIGHT <- 22.0

SAMPLE_FONT_SIZE <- BASE_SAMPLE_FONT_SIZE * BOLDNESS_MULTIPLIER
COLUMN_FONT_SCALE <- 1.00
COLUMN_SAMPLE_FONT_SIZE <- SAMPLE_FONT_SIZE * COLUMN_FONT_SCALE
ANNOTATION_FONT_SIZE <- BASE_ANNOTATION_FONT_SIZE * BOLDNESS_MULTIPLIER
LEGEND_FONT_SIZE <- BASE_LEGEND_FONT_SIZE * BOLDNESS_MULTIPLIER

# 图例放在整体图片左侧，并与热图主体保持适度间距。
# 这些参数只控制整体排版，不改变热图本体的聚类结果。
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
source(PLOTTING_FUNCTION_FILE)


# 2. 读取数据 ------------------------------------------------------------------

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


# 3. 准备TPM表达矩阵 -----------------------------------------------------------

assay_names <- names(SummarizedExperiment::assays(se))
stopifnot(TPM_ASSAY_NAME %in% assay_names)

tpm_all <- as.matrix(SummarizedExperiment::assay(se, TPM_ASSAY_NAME))
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

# 输出目录用基因类型区分；PDF文件名只保留结果类型。
plot_filter_name <- sanitize_file_name(gene_filter$filter)
PLOT_DIR <- file.path(PLOT_ROOT, plot_filter_name)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

if (CLEAN_PLOT_OUTPUT_DIR) {
  unlink(list.files(
    PLOT_DIR,
    pattern = "[.]pdf$",
    recursive = TRUE,
    full.names = TRUE
  ))

  legacy_pdf_files <- list.files(
    PLOT_ROOT,
    pattern = "[.]pdf$",
    full.names = TRUE
  )
  if (length(legacy_pdf_files) > 0) {
    unlink(legacy_pdf_files)
  }
}

tpm_filtered <- gene_filter$exprSet


# 4. 单组样本聚类热图绘制函数 --------------------------------------------------

draw_sample_clustering_heatmap <- function(sample_group_name, sample_group_value) {
  # 对一个样本集合完成筛选、相关性计算、层次聚类、PDF排版和保存。
  # 聚类只基于TPM表达谱得到的样本相关性，不按细胞系分组强行排序。
  sample_status <- trimws(as.character(sample_info_all[[LRC_COLUMN]]))
  sample_status[is.na(sample_status)] <- ""
  sample_index <- sample_status == sample_group_value

  sample_info <- sample_info_all[sample_index, , drop = FALSE]
  stopifnot(nrow(sample_info) >= 2)

  sample_labels <- get_display_labels(
    sample_info = sample_info,
    label_column = SAMPLE_LABEL_COLUMN
  )

  tpm <- tpm_filtered[, sample_info$Sample_ID, drop = FALSE]
  correlation_result <- prepare_sample_correlation(
    expr_matrix = tpm,
    correlation_method = CORRELATION_METHOD,
    clustering_method = CLUSTERING_METHOD
  )
  expr_for_correlation <- correlation_result$expr_for_correlation
  cor_matrix <- correlation_result$cor_matrix
  sample_hclust <- correlation_result$sample_hclust

  plot_row_labels <- wrap_label(sample_labels, ROW_LABEL_WIDTH)
  plot_column_labels <- wrap_label_by_underscore(sample_labels, COLUMN_LABEL_WIDTH)

  label_padding <- paste(rep(" ", LABEL_HEATMAP_GAP_SPACES), collapse = "")
  plot_row_labels <- paste0(label_padding, plot_row_labels)
  plot_column_labels <- paste0(label_padding, plot_column_labels)

  rownames(cor_matrix) <- sample_info$Sample_ID
  colnames(cor_matrix) <- sample_info$Sample_ID

  cell_line_palette <- get_named_brewer_palette(sample_info$cell_line)
  cell_line_levels <- names(cell_line_palette)

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

  output_dir <- file.path(PLOT_DIR, sanitize_file_name(sample_group_name))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(output_dir, "heatmap.pdf")

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
    show_heatmap_legend = FALSE
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
  stopifnot(file.exists(pdf_file))

  data.frame(
    Sample_Group = sample_group_name,
    Group_Value = sample_group_value,
    Samples = nrow(sample_info),
    Genes_After_Biotype_Filter = gene_filter$selected_gene_count,
    Genes_After_Removing_Zero_Variance = nrow(expr_for_correlation),
    PDF_Width = round(pdf_width, 2),
    PDF_Height = round(pdf_height, 2),
    PDF_File = pdf_file,
    stringsAsFactors = FALSE
  )
}


# 5. 依次绘制LRC和BULK样本聚类热图 -------------------------------------------

summary_list <- lapply(names(SAMPLE_GROUPS), function(sample_group_name) {
  draw_sample_clustering_heatmap(
    sample_group_name = sample_group_name,
    sample_group_value = SAMPLE_GROUPS[[sample_group_name]]
  )
})

summary_table <- do.call(rbind, summary_list)
rownames(summary_table) <- NULL


# 6. 输出信息 ------------------------------------------------------------------

cat("\nSample TPM clustering heatmaps finished.\n")
cat("Gene biotype filter: ", gene_filter$filter, "\n", sep = "")
cat("TPM assay: ", TPM_ASSAY_NAME, "\n", sep = "")
cat("\nHeatmap summary:\n")
print(
  summary_table[
    ,
    c(
      "Sample_Group", "Samples",
      "Genes_After_Biotype_Filter", "Genes_After_Removing_Zero_Variance",
      "PDF_Width", "PDF_Height"
    )
  ],
  row.names = FALSE
)

cat("\nPDF files:\n")
cat(paste(summary_table$PDF_File, collapse = "\n"), "\n", sep = "")
