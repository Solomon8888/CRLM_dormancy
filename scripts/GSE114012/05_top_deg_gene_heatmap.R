# GSE114012 Top差异表达基因热图
#
# 自动读取01号脚本输出的每一套significant_genes.csv，
# 提取其中按P值排序的Top差异基因，并使用SE对象中的TPM表达量绘制表达热图。
# 热图参考ComplexHeatmap复杂热图风格：顶部为样本分组和细胞系注释，
# 左侧用细条标记Up/Down方向，主体展示log2(TPM + 1)后的行Z-score。


# 0. 可修改配置 ---------------------------------------------------------------

# 当前脚本只服务于GSE114012这个NGS数据集；目录结构也依赖这两个字段。
DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

# 输入文件：SE对象提供TPM矩阵，临床表提供样本分组、细胞系和Title。
SE_RDS_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_clinical_edit.csv"

# 01号脚本用于同步P值列名；FUNCTION_FILE用于解析差异分析设计；
# PLOTTING_FUNCTION_FILE保存跨绘图脚本共用的风格配置和PDF/PNG输出函数。
DE_SCRIPT_FILE <- "scripts/GSE114012/01_limma_differential_expression.R"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"

# 输入目录为01号脚本输出的tables/<analysis_name>/DEG/significant_genes.csv。
# 输出目录为plots/gene_heatmap/<analysis_name>/gene_heatmap.*。
RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "gene_heatmap")

# 是否从01号差异分析脚本中同步P值列名。
SYNC_P_VALUE_COLUMN_FROM_01_SCRIPT <- TRUE
P_VALUE_COLUMN <- "P.Value"

# 需要绘制哪些分析。设为"all"时自动绘制全部significant_genes.csv。
ANALYSES_TO_PLOT <- "all"
# ANALYSES_TO_PLOT <- c("DLD1", "HCT15")

# 重跑时清理整个gene_heatmap目录，保证最终只保留当前配置生成的图片。
CLEAN_GENE_HEATMAP_ROOT <- TRUE

# SE对象中TPM矩阵的assay名称。基因表达热图使用TPM，不使用原始count。
TPM_ASSAY_NAME <- "tpm"

# 每个分析设计最多展示Top多少个显著差异基因。
# 若某个分析的显著基因不足该数量，则展示全部显著基因。
TOP_GENE_COUNT <- 50

# 基因ID、基因名、logFC列名。01号脚本最终DEG表不再保存Feature_ID。
GENE_ID_COLUMN <- "GeneID"
GENE_SYMBOL_COLUMN <- "Symbol"
LOGFC_COLUMN <- "logFC"

# 样本名显示列；若Title为空，会回退到Sample_ID。
SAMPLE_LABEL_COLUMN <- "Title"
CELL_LINE_COLUMN <- "cell_line"

# 热图数值处理：log2(TPM + 1)后按基因做Z-score，并截断极端值。
Z_SCORE_MIN <- -2
Z_SCORE_MAX <- 2

# 行列聚类设置。行按Up/Down分块但不显示文字标题；列按LRC、BULK分块，组内按真实表达模式聚类。
CLUSTER_ROWS <- FALSE
CLUSTER_COLUMNS <- TRUE
ROW_SPLIT_BY_DIRECTION <- TRUE
COLUMN_SPLIT_BY_GROUP <- TRUE
CLUSTERING_METHOD <- "complete"

# 是否显示基因名和样本名。样本名不换行，斜45度完整显示。
SHOW_ROW_NAMES <- TRUE
SHOW_COLUMN_NAMES <- TRUE
ROW_LABEL_WIDTH <- 24
COLUMN_NAMES_ROT <- 45
ROW_NAME_LEFT_PADDING_SPACES <- 1
COLUMN_NAME_LEFT_PADDING_SPACES <- 2

# 热图颜色：深蓝-白-鲜红，对应低表达Z-score、中位、 高表达Z-score。
HEATMAP_COLOR_LOW <- "#2166AC"
HEATMAP_COLOR_MID <- "#FFFFFF"
HEATMAP_COLOR_HIGH <- "#D73027"

# 样本分组颜色。每个分析的实验组会自动使用红色；BULK使用深灰。
CONTROL_GROUP_COLOR <- "#6A6A6A"
EXPERIMENT_GROUP_COLOR <- "#D73027"
DISPLAY_CONTROL_GROUP <- "BULK"
DISPLAY_EXPERIMENT_GROUP <- "LRC"
DIRECTION_COLORS <- c(
  Up = "#D73027",
  Down = "#2166AC"
)

# 图片尺寸动态参数。样本数决定宽度，基因数决定高度。
# 通过固定单元格宽高比例并动态计算PDF画布，使全部热图保持一致视觉比例。
CELL_WIDTH_MM <- 10.4
CELL_HEIGHT_MM <- 3.8
ROW_DEND_WIDTH_MM <- 18
COLUMN_DEND_HEIGHT_MM <- 9
TOP_ANNOTATION_SIZE_MM <- 1.35
TOP_ANNOTATION_INTERNAL_GAP_MM <- 1.2
TOP_TO_HEATMAP_GAP_MM <- 0.35
DIRECTION_ANNOTATION_SIZE_MM <- 1.6
DIRECTION_TO_HEATMAP_GAP_MM <- 0.35
RIGHT_ROW_NAME_MAX_WIDTH_INCH <- 1.45
OUTER_MARGIN_INCH <- 0.38
LEGEND_WIDTH_INCH <- 1.50
LEGEND_TITLE_GAP_MM <- 3.0
HEATMAP_TO_LEGEND_GAP_MM <- 6.0

MIN_PDF_WIDTH <- 7.0
MAX_PDF_WIDTH <- 32.0
MIN_PDF_HEIGHT <- 6.0
MAX_PDF_HEIGHT <- 24.0

# 字体和线条。沿用公共配置中的Helvetica粗体，并加粗热图框线。
HEATMAP_FONT_MULTIPLIER <- 1.30
ROW_NAME_FONT_SIZE <- 5.8 * HEATMAP_FONT_MULTIPLIER
COLUMN_NAME_FONT_SIZE <- 6.1 * HEATMAP_FONT_MULTIPLIER
ANNOTATION_FONT_SIZE <- 7.0 * HEATMAP_FONT_MULTIPLIER
LEGEND_FONT_SIZE <- 7.2 * HEATMAP_FONT_MULTIPLIER
COLUMN_TITLE_FONT_SIZE <- 10.0 * HEATMAP_FONT_MULTIPLIER
CELL_BORDER_WIDTH <- 1.80
DENDROGRAM_LINE_WIDTH <- 1.90

options(width = 200)


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

source(FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

# 让Top基因热图的红蓝体系明确跟随火山图脚本的公共显著基因配色。
HEATMAP_COLOR_LOW <- DOWN_COLOR
HEATMAP_COLOR_HIGH <- UP_COLOR
EXPERIMENT_GROUP_COLOR <- UP_COLOR
DIRECTION_COLORS <- c(
  Up = UP_COLOR,
  Down = DOWN_COLOR
)


# 2. 常用函数 -----------------------------------------------------------------

get_significant_file_info <- function(table_root, deg_dir_name = "DEG") {
  significant_files <- list.files(
    table_root,
    pattern = "^significant_genes[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  significant_files <- significant_files[
    basename(dirname(significant_files)) == deg_dir_name |
      (
        basename(dirname(significant_files)) == "csv" &
          basename(dirname(dirname(significant_files))) == deg_dir_name
      )
  ]

  get_analysis_name <- function(file_name) {
    if (basename(dirname(file_name)) == "csv") {
      return(basename(dirname(dirname(dirname(file_name)))))
    }

    basename(dirname(dirname(file_name)))
  }

  if (exists("prefer_report_csv_files", mode = "function")) {
    significant_files <- prefer_report_csv_files(significant_files, get_analysis_name)
  }

  stopifnot(length(significant_files) > 0)

  file_info <- data.frame(
    Analysis_Name = vapply(significant_files, get_analysis_name, character(1)),
    Significant_File = significant_files,
    stringsAsFactors = FALSE
  )

  file_info <- file_info[order(file_info$Analysis_Name), , drop = FALSE]
  rownames(file_info) <- NULL

  if (any(duplicated(file_info$Analysis_Name))) {
    duplicated_names <- unique(file_info$Analysis_Name[
      duplicated(file_info$Analysis_Name)
    ])
    stop(
      "More than one significant_genes.csv file was found for: ",
      paste(duplicated_names, collapse = ", ")
    )
  }

  file_info
}

read_significant_result <- function(file_info, analysis_name) {
  file_index <- match(analysis_name, file_info$Analysis_Name)
  if (is.na(file_index)) {
    stop("No significant_genes.csv file was found for: ", analysis_name)
  }

  read.csv(
    if (exists("resolve_report_csv_file", mode = "function")) {
      resolve_report_csv_file(file_info$Significant_File[file_index])
    } else {
      file_info$Significant_File[file_index]
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

select_top_deg_genes <- function(dat) {
  stopifnot(GENE_ID_COLUMN %in% colnames(dat))
  stopifnot(LOGFC_COLUMN %in% colnames(dat))
  stopifnot(P_VALUE_COLUMN %in% colnames(dat))

  dat[[LOGFC_COLUMN]] <- as.numeric(dat[[LOGFC_COLUMN]])
  dat[[P_VALUE_COLUMN]] <- as.numeric(dat[[P_VALUE_COLUMN]])

  valid_index <- !is.na(dat[[GENE_ID_COLUMN]]) &
    dat[[GENE_ID_COLUMN]] != "" &
    is.finite(dat[[LOGFC_COLUMN]]) &
    is.finite(dat[[P_VALUE_COLUMN]])

  dat <- dat[valid_index, , drop = FALSE]
  dat <- dat[
    order(dat[[P_VALUE_COLUMN]], -abs(dat[[LOGFC_COLUMN]])),
    ,
    drop = FALSE
  ]
  dat <- dat[!duplicated(dat[[GENE_ID_COLUMN]]), , drop = FALSE]

  dat[seq_len(min(TOP_GENE_COUNT, nrow(dat))), , drop = FALSE]
}

make_row_labels <- function(top_deg) {
  if (GENE_SYMBOL_COLUMN %in% colnames(top_deg)) {
    row_labels <- trimws(as.character(top_deg[[GENE_SYMBOL_COLUMN]]))
    empty_label <- is.na(row_labels) | row_labels == ""
    row_labels[empty_label] <- top_deg[[GENE_ID_COLUMN]][empty_label]
  } else {
    row_labels <- top_deg[[GENE_ID_COLUMN]]
  }

  row_labels <- make.unique(wrap_label(row_labels, ROW_LABEL_WIDTH))
  paste0(strrep(" ", ROW_NAME_LEFT_PADDING_SPACES), row_labels)
}

prepare_expression_zscore <- function(tpm_matrix, gene_ids, sample_ids) {
  expr <- log2(tpm_matrix[gene_ids, sample_ids, drop = FALSE] + 1)
  expr_z <- t(scale(t(expr)))
  expr_z[!is.finite(expr_z)] <- 0
  expr_z <- pmax(pmin(expr_z, Z_SCORE_MAX), Z_SCORE_MIN)
  expr_z
}

get_analysis_design_row <- function(analysis_name, analysis_designs) {
  design_index <- match(analysis_name, analysis_designs$Analysis_Name)
  if (is.na(design_index)) {
    stop("No clinical analysis design was found for: ", analysis_name)
  }

  analysis_designs[design_index, , drop = FALSE]
}

get_pdf_size <- function(n_genes, n_samples, max_column_label_chars, max_column_label_lines) {
  heatmap_width <- n_samples * CELL_WIDTH_MM / 25.4
  heatmap_height <- n_genes * CELL_HEIGHT_MM / 25.4

  column_label_space <- if (SHOW_COLUMN_NAMES) {
    max_column_label_chars * COLUMN_NAME_FONT_SIZE * 0.010 +
      max_column_label_lines * COLUMN_NAME_FONT_SIZE * 0.018
  } else {
    0
  }

  row_name_space <- if (SHOW_ROW_NAMES) {
    RIGHT_ROW_NAME_MAX_WIDTH_INCH
  } else {
    0
  }

  pdf_width <- OUTER_MARGIN_INCH * 2 +
    DIRECTION_ANNOTATION_SIZE_MM / 25.4 +
    DIRECTION_TO_HEATMAP_GAP_MM / 25.4 +
    ROW_DEND_WIDTH_MM / 25.4 +
    heatmap_width +
    row_name_space +
    HEATMAP_TO_LEGEND_GAP_MM / 25.4 +
    LEGEND_WIDTH_INCH

  pdf_height <- OUTER_MARGIN_INCH * 2 +
    COLUMN_DEND_HEIGHT_MM / 25.4 +
    (TOP_ANNOTATION_SIZE_MM * 2 +
      TOP_ANNOTATION_INTERNAL_GAP_MM +
      TOP_TO_HEATMAP_GAP_MM) / 25.4 +
    heatmap_height +
    column_label_space

  list(
    width = min(max(pdf_width, MIN_PDF_WIDTH), MAX_PDF_WIDTH),
    height = min(max(pdf_height, MIN_PDF_HEIGHT), MAX_PDF_HEIGHT)
  )
}

get_annotation_colors <- function(sample_info, group_column) {
  group_colors <- c(
    EXPERIMENT_GROUP_COLOR,
    CONTROL_GROUP_COLOR
  )
  names(group_colors) <- c(DISPLAY_EXPERIMENT_GROUP, DISPLAY_CONTROL_GROUP)
  group_colors <- group_colors[unique(as.character(sample_info[[group_column]]))]

  cell_line_colors <- get_named_brewer_palette(sample_info[[CELL_LINE_COLUMN]])

  list(
    group_colors = group_colors,
    cell_line_colors = cell_line_colors
  )
}

make_top_annotation <- function(sample_info, group_column, annotation_colors) {
  ComplexHeatmap::HeatmapAnnotation(
    Group = sample_info[[group_column]],
    Group_Cell_line_gap = ComplexHeatmap::anno_empty(
      height = grid::unit(TOP_ANNOTATION_INTERNAL_GAP_MM, "mm"),
      border = FALSE
    ),
    Cell_line = sample_info[[CELL_LINE_COLUMN]],
    Heatmap_gap = ComplexHeatmap::anno_empty(
      height = grid::unit(TOP_TO_HEATMAP_GAP_MM, "mm"),
      border = FALSE
    ),
    col = list(
      Group = annotation_colors$group_colors,
      Cell_line = annotation_colors$cell_line_colors
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(TOP_ANNOTATION_SIZE_MM, "mm"),
    gap = grid::unit(0, "mm"),
    gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
    border = TRUE
  )
}

make_direction_annotation <- function(direction) {
  ComplexHeatmap::rowAnnotation(
    Direction = direction,
    Heatmap_gap = ComplexHeatmap::anno_empty(
      width = grid::unit(DIRECTION_TO_HEATMAP_GAP_MM, "mm"),
      border = FALSE
    ),
    col = list(Direction = DIRECTION_COLORS),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(DIRECTION_ANNOTATION_SIZE_MM, "mm"),
    gap = grid::unit(0, "mm"),
    gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
    border = TRUE
  )
}

make_heatmap_legends <- function(annotation_colors, heatmap_colors) {
  legend_text_gp <- grid::gpar(
    fontsize = LEGEND_FONT_SIZE,
    fontface = TEXT_FONT_FACE,
    fontfamily = TEXT_FONT_FAMILY,
    col = TEXT_COLOR
  )

  list(
    ComplexHeatmap::Legend(
      title = "Group",
      at = names(annotation_colors$group_colors),
      labels = names(annotation_colors$group_colors),
      type = "grid",
      legend_gp = grid::gpar(fill = annotation_colors$group_colors, col = "black"),
      title_gp = legend_text_gp,
      labels_gp = legend_text_gp,
      title_gap = grid::unit(LEGEND_TITLE_GAP_MM, "mm"),
      direction = "vertical"
    ),
    ComplexHeatmap::Legend(
      title = "Cell_Line",
      at = names(annotation_colors$cell_line_colors),
      labels = names(annotation_colors$cell_line_colors),
      type = "grid",
      legend_gp = grid::gpar(fill = annotation_colors$cell_line_colors, col = "black"),
      title_gp = legend_text_gp,
      labels_gp = legend_text_gp,
      title_gap = grid::unit(LEGEND_TITLE_GAP_MM, "mm"),
      direction = "vertical"
    ),
    ComplexHeatmap::Legend(
      title = "Direction",
      at = c("Up", "Down"),
      labels = c("Up", "Down"),
      type = "grid",
      legend_gp = grid::gpar(
        fill = DIRECTION_COLORS[c("Up", "Down")],
        col = "black"
      ),
      title_gp = legend_text_gp,
      labels_gp = legend_text_gp,
      title_gap = grid::unit(LEGEND_TITLE_GAP_MM, "mm"),
      direction = "vertical"
    ),
    ComplexHeatmap::Legend(
      title = "Z-Score",
      at = pretty(c(Z_SCORE_MIN, Z_SCORE_MAX), n = 5),
      col_fun = heatmap_colors,
      border = "black",
      title_gp = legend_text_gp,
      labels_gp = legend_text_gp,
      title_gap = grid::unit(LEGEND_TITLE_GAP_MM, "mm"),
      direction = "vertical"
    )
  )
}

draw_top_deg_gene_heatmap <- function(analysis_name, file_info, analysis_designs, tpm_matrix, sample_info_all) {
  significant_result <- read_significant_result(file_info, analysis_name)
  top_deg <- select_top_deg_genes(significant_result)
  stopifnot(nrow(top_deg) > 0)

  gene_ids <- top_deg[[GENE_ID_COLUMN]]
  gene_ids <- gene_ids[gene_ids %in% rownames(tpm_matrix)]
  top_deg <- top_deg[top_deg[[GENE_ID_COLUMN]] %in% gene_ids, , drop = FALSE]
  gene_ids <- top_deg[[GENE_ID_COLUMN]]
  stopifnot(length(gene_ids) > 0)

  analysis_design <- get_analysis_design_row(analysis_name, analysis_designs)
  design_samples <- prepare_design_samples(
    sample_info = sample_info_all,
    group_column_index = analysis_design$Column_Index,
    experiment_group = analysis_design$Experiment_Group
  )

  sample_info <- design_samples$sample_info
  group_list <- design_samples$group_list
  control_group <- design_samples$control_group

  sample_info$Display_Group <- ifelse(
    as.character(group_list) == control_group,
    DISPLAY_CONTROL_GROUP,
    DISPLAY_EXPERIMENT_GROUP
  )
  sample_info$Display_Group <- factor(
    sample_info$Display_Group,
    levels = c(DISPLAY_EXPERIMENT_GROUP, DISPLAY_CONTROL_GROUP)
  )

  sample_order <- order(sample_info$Display_Group)
  sample_info <- sample_info[sample_order, , drop = FALSE]

  sample_labels <- get_display_labels(
    sample_info = sample_info,
    label_column = SAMPLE_LABEL_COLUMN
  )
  plot_column_labels <- paste0(
    strrep(" ", COLUMN_NAME_LEFT_PADDING_SPACES),
    sample_labels
  )

  direction <- ifelse(top_deg[[LOGFC_COLUMN]] > 0, "Up", "Down")
  direction <- factor(direction, levels = c("Up", "Down"))
  row_order <- order(
    direction,
    top_deg[[P_VALUE_COLUMN]],
    -abs(top_deg[[LOGFC_COLUMN]])
  )
  top_deg <- top_deg[row_order, , drop = FALSE]
  gene_ids <- top_deg[[GENE_ID_COLUMN]]
  direction <- direction[row_order]

  expr_z <- prepare_expression_zscore(
    tpm_matrix = tpm_matrix,
    gene_ids = gene_ids,
    sample_ids = sample_info$Sample_ID
  )

  row_labels <- make_row_labels(top_deg)
  rownames(expr_z) <- row_labels
  colnames(expr_z) <- sample_info$Sample_ID

  heatmap_colors <- circlize::colorRamp2(
    c(Z_SCORE_MIN, 0, Z_SCORE_MAX),
    c(HEATMAP_COLOR_LOW, HEATMAP_COLOR_MID, HEATMAP_COLOR_HIGH)
  )

  annotation_colors <- get_annotation_colors(
    sample_info = sample_info,
    group_column = "Display_Group"
  )
  top_annotation <- make_top_annotation(
    sample_info = sample_info,
    group_column = "Display_Group",
    annotation_colors = annotation_colors
  )
  direction_annotation <- make_direction_annotation(direction)
  legend_list <- make_heatmap_legends(
    annotation_colors = annotation_colors,
    heatmap_colors = heatmap_colors
  )

  column_label_parts <- strsplit(plot_column_labels, "\n", fixed = TRUE)
  max_column_label_chars <- max(nchar(unlist(column_label_parts)))
  max_column_label_lines <- max(lengths(column_label_parts))
  pdf_size <- get_pdf_size(
    n_genes = nrow(expr_z),
    n_samples = ncol(expr_z),
    max_column_label_chars = max_column_label_chars,
    max_column_label_lines = max_column_label_lines
  )

  column_split <- if (COLUMN_SPLIT_BY_GROUP) sample_info$Display_Group else NULL
  row_split <- if (ROW_SPLIT_BY_DIRECTION) direction else NULL

  ht <- ComplexHeatmap::Heatmap(
    expr_z,
    name = "Z-Score",
    col = heatmap_colors,
    top_annotation = top_annotation,
    left_annotation = direction_annotation,
    cluster_rows = CLUSTER_ROWS,
    cluster_columns = CLUSTER_COLUMNS,
    clustering_method_rows = CLUSTERING_METHOD,
    clustering_method_columns = CLUSTERING_METHOD,
    row_split = row_split,
    column_split = column_split,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    row_title = NULL,
    column_title_gp = grid::gpar(
      fontsize = COLUMN_TITLE_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    show_row_names = SHOW_ROW_NAMES,
    show_column_names = SHOW_COLUMN_NAMES,
    row_names_side = "right",
    row_names_gp = grid::gpar(
      fontsize = ROW_NAME_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    column_labels = plot_column_labels,
    column_names_gp = grid::gpar(
      fontsize = COLUMN_NAME_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    column_names_rot = COLUMN_NAMES_ROT,
    row_names_max_width = grid::unit(RIGHT_ROW_NAME_MAX_WIDTH_INCH, "in"),
    rect_gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
    border = TRUE,
    border_gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
    row_dend_width = grid::unit(ROW_DEND_WIDTH_MM, "mm"),
    column_dend_height = grid::unit(COLUMN_DEND_HEIGHT_MM, "mm"),
    row_dend_gp = grid::gpar(col = "black", lwd = DENDROGRAM_LINE_WIDTH),
    column_dend_gp = grid::gpar(col = "black", lwd = DENDROGRAM_LINE_WIDTH),
    width = grid::unit(ncol(expr_z) * CELL_WIDTH_MM, "mm"),
    height = grid::unit(nrow(expr_z) * CELL_HEIGHT_MM, "mm"),
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(
      title = "Z-Score",
      border = "black",
      title_gp = grid::gpar(
        fontsize = LEGEND_FONT_SIZE,
        fontface = TEXT_FONT_FACE,
        fontfamily = TEXT_FONT_FAMILY,
        col = TEXT_COLOR
      ),
      labels_gp = grid::gpar(
        fontsize = LEGEND_FONT_SIZE,
        fontface = TEXT_FONT_FACE,
        fontfamily = TEXT_FONT_FAMILY,
        col = TEXT_COLOR
      )
    )
  )

  output_dir <- file.path(PLOT_ROOT, sanitize_file_name(analysis_name))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(output_dir, pattern = "[.](pdf|png)$", full.names = TRUE))

  pdf_file <- file.path(output_dir, "gene_heatmap.pdf")
  output_files <- save_grid_pdf_png(
    pdf_file = pdf_file,
    width = pdf_size$width,
    height = pdf_size$height,
    draw_fun = function() {
      ComplexHeatmap::ht_opt(
        HEATMAP_LEGEND_PADDING = grid::unit(HEATMAP_TO_LEGEND_GAP_MM, "mm")
      )
      ComplexHeatmap::draw(
        ht,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        show_heatmap_legend = TRUE,
        show_annotation_legend = FALSE,
        heatmap_legend_list = legend_list,
        merge_legends = TRUE,
        legend_grouping = "original",
        align_heatmap_legend = "heatmap_top"
      )
    }
  )

  data.frame(
    Analysis_Name = analysis_name,
    Samples = nrow(sample_info),
    Significant_Genes = nrow(significant_result),
    Genes_Plotted = nrow(expr_z),
    Up = sum(direction == "Up"),
    Down = sum(direction == "Down"),
    PDF_Width = round(pdf_size$width, 2),
    PDF_Height = round(pdf_size$height, 2),
    PDF_File = output_files$pdf_file,
    PNG_File = output_files$png_file,
    stringsAsFactors = FALSE
  )
}


# 3. 读取数据并同步配置 --------------------------------------------------------

threshold_config <- sync_de_thresholds_from_script(
  script_file = DE_SCRIPT_FILE,
  p_value_column = P_VALUE_COLUMN,
  p_value_cutoff = 0.05,
  logfc_cutoff = 0.5,
  sync = SYNC_P_VALUE_COLUMN_FROM_01_SCRIPT
)
P_VALUE_COLUMN <- threshold_config$p_value_column

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot(TPM_ASSAY_NAME %in% names(SummarizedExperiment::assays(se)))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot(SAMPLE_LABEL_COLUMN %in% colnames(clinical_data))
stopifnot(CELL_LINE_COLUMN %in% colnames(clinical_data))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

tpm_matrix <- as.matrix(SummarizedExperiment::assay(se, TPM_ASSAY_NAME))
stopifnot(is.numeric(tpm_matrix))

missing_samples <- setdiff(colnames(tpm_matrix), clinical_data$Sample_ID)
stopifnot(length(missing_samples) == 0)

sample_info_all <- clinical_data[
  match(colnames(tpm_matrix), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
rownames(sample_info_all) <- sample_info_all$Sample_ID
stopifnot(all(sample_info_all$Sample_ID == colnames(tpm_matrix)))

analysis_designs <- get_analysis_designs(clinical_data)
file_info <- get_significant_file_info(TABLE_ROOT)
selected_analyses <- get_selected_analysis_names(
  file_info = data.frame(
    Analysis_Name = file_info$Analysis_Name,
    All_Genes_File = file_info$Significant_File,
    stringsAsFactors = FALSE
  ),
  analyses_to_plot = ANALYSES_TO_PLOT
)

missing_designs <- setdiff(selected_analyses, analysis_designs$Analysis_Name)
if (length(missing_designs) > 0) {
  stop("No clinical analysis design was found for: ", paste(missing_designs, collapse = ", "))
}


# 4. 绘制每套差异分析设计的Top差异基因热图 -----------------------------------

if (CLEAN_GENE_HEATMAP_ROOT && dir.exists(PLOT_ROOT)) {
  unlink(PLOT_ROOT, recursive = TRUE)
}
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)

cat("\nRunning top DEG gene heatmap generation...\n")
cat("P value column: ", P_VALUE_COLUMN, "\n", sep = "")
cat("TPM assay: ", TPM_ASSAY_NAME, "\n", sep = "")
cat("Top DEG count: ", TOP_GENE_COUNT, "\n", sep = "")

run_one_top_deg_heatmap <- function(i) {
  analysis_name <- selected_analyses[i]
  draw_top_deg_gene_heatmap(
    analysis_name = analysis_name,
    file_info = file_info,
    analysis_designs = analysis_designs,
    tpm_matrix = tpm_matrix,
    sample_info_all = sample_info_all
  )
}

parallel_strategy <- setup_parallel_strategy(
  total_tasks = length(selected_analyses),
  inner_label = "Gene heatmap inner workers",
  nested_label = "Nested workers"
)

summary_list <- run_indexed_tasks_with_progress(
  total_tasks = length(selected_analyses),
  workers = parallel_strategy$task_workers,
  task_function = run_one_top_deg_heatmap
)
stop_on_parallel_errors(summary_list, task_ids = selected_analyses, label = "top DEG heatmaps")
names(summary_list) <- selected_analyses


# 5. 终端快速汇总 --------------------------------------------------------------

summary_table <- do.call(rbind, summary_list)
rownames(summary_table) <- NULL

cat("\nTop DEG gene heatmap summary:\n")
print(
  summary_table[
    ,
    c(
      "Analysis_Name", "Samples", "Significant_Genes", "Genes_Plotted",
      "Up", "Down", "PDF_Width", "PDF_Height"
    )
  ],
  row.names = FALSE
)

cat("\nTop DEG gene heatmap PDF files:\n")
cat(paste(summary_table$PDF_File, collapse = "\n"), "\n", sep = "")
cat("\nTop DEG gene heatmap PNG files:\n")
cat(paste(summary_table$PNG_File, collapse = "\n"), "\n", sep = "")
cat("\nTop DEG gene heatmap generation finished.\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
