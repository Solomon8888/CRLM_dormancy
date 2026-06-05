# GSE114012批量GSEA绘图
#
# 读取06号脚本生成的GSEA表格与qs2缓存对象，
# 为每套analysis x geneset绘制GseaVis气泡图和单通路GSEA+表达热图。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

SE_RDS_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_clinical_edit.csv"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
GSEA_FUNCTION_FILE <- "scripts/functions/gsea_common_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")

# 测试时可临时设置GSEA_OUTPUT_ROOT，与06号脚本保持一致。
OUTPUT_ROOT <- Sys.getenv("GSEA_OUTPUT_ROOT", unset = RESULT_ROOT)
TABLE_OUTPUT_ROOT <- file.path(OUTPUT_ROOT, "tables")
PLOT_ROOT <- file.path(OUTPUT_ROOT, "plots", "GSEA")

# 需要绘制哪些差异分析设计。设为"all"时自动绘制全部DEG/all_genes.csv对应结果。
ANALYSES_TO_RUN <- "all"
# ANALYSES_TO_RUN <- c("DLD1", "DLD1_HCT15_SW48")

# 物种和排序配置必须与06号GSEA运算脚本保持一致，用于定位qs2缓存对象。
SPECIES <- "human"
GENE_ID_TYPE <- "ENTREZ"  # 可选："ENTREZ", "SYMBOL", "ENSEMBL"
RANK_METRIC_COLUMN <- "t"

GSEA_PARAMS <- list(
  exponent = 1,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  nPerm = 1000,
  method = "multilevel",
  adaptive = FALSE,
  minPerm = 101,
  maxPerm = 1e5,
  pvalThreshold = 0.1
)

GSEA_GENESETS_TO_RUN <- "all"
# GSEA_GENESETS_TO_RUN <- c("H", "C5:GO:BP", "C6")

# 重跑07号脚本时清空旧GSEA图片，避免新旧图片混在一起。
CLEAN_GSEA_PLOT_OUTPUT_DIR <- TRUE

# GseaVis::dotplotGsea气泡图配置。
SIMPLIFY_PATHWAY_PREFIX_IN_PLOT <- TRUE
REPLACE_UNDERSCORE_WITH_SPACE_IN_PLOT <- TRUE

GSEAVIS_DOTPLOT_PARAMS <- list(
  topn = 20,
  pval = NULL,
  pajust = 0.05,
  order.by = "GeneRatio",
  str.width = 45,
  base_size = 10,
  scales = "free_x",
  add.seg = FALSE,
  line.col = "grey80",
  line.size = 1.5,
  line.type = "solid"
)

GSEAVIS_POINT_SIZE_RANGE <- c(4.2, 9.8)

DOTPLOT_BODY_BASE_SIZE <- 4.8
DOTPLOT_LABEL_LINE_HEIGHT <- 0.34
DOTPLOT_TERM_GAP_HEIGHT <- 0.14
DOTPLOT_BODY_MIN_SIZE <- 5.2
DOTPLOT_BODY_MAX_SIZE <- 18.0
DOTPLOT_LABEL_BASE_WIDTH <- 1.8
DOTPLOT_LABEL_WIDTH_PER_CHARACTER <- 0.045
DOTPLOT_LABEL_MIN_WIDTH <- 2.4
DOTPLOT_LABEL_MAX_WIDTH <- 6.2
DOTPLOT_LEGEND_WIDTH <- 1.4
DOTPLOT_VERTICAL_PADDING <- 0.45

# 单通路GSEA图配置。默认绘制每套GSEA结果的Top20通路。
DRAW_SINGLE_PATHWAY_GSEA <- TRUE
SINGLE_PATHWAY_TOP_N <- GSEAVIS_DOTPLOT_PARAMS$topn
SINGLE_PATHWAY_KEYWORDS <- character(0)
SINGLE_PATHWAY_MAX_KEYWORD_TERMS <- 20
SINGLE_PATHWAY_PVALUE_COLUMN <- "p.adjust"
SINGLE_PATHWAY_PVALUE_CUTOFF <- 0.05

TPM_ASSAY_NAME <- "tpm"
SINGLE_PATHWAY_SAMPLE_LABEL_COLUMN <- "Title"

GSEAVIS_SINGLE_PATHWAY_PARAMS <- list(
  subPlot = 2,
  lineSize = 1.0,
  rmSegment = FALSE,
  termWidth = 45,
  segCol = "black",
  curveCol = c("#2166AC", "#D73027", "#762A83"),
  htCol = c("#2166AC", "#D73027"),
  rankCol = c("#2166AC", "white", "#D73027"),
  rankSeq = 5000,
  htHeight = 0.30,
  force = 20,
  max.overlaps = 50,
  geneSize = 4,
  newGsea = TRUE,
  addPoint = TRUE,
  newCurveCol = c("#2166AC", "white", "#D73027"),
  newHtCol = c("#2166AC", "white", "#D73027"),
  rmHt = FALSE,
  addPval = TRUE,
  pvalX = 0.96,
  pvalY = 0.92,
  pvalSize = 4,
  pCol = "black",
  pHjust = 1,
  rmPrefix = TRUE,
  nesDigit = 2,
  pDigit = 3,
  markTopgene = FALSE,
  topGeneN = 5,
  legend.position = "right",
  add.geneExpHt = TRUE,
  scale.exp = TRUE,
  exp.col = c("#2166AC", "white", "#D73027"),
  ht.legend = TRUE,
  ght.relHight = 0.42,
  ght.geneText.size = 7,
  ght.facet = FALSE,
  ght.facet.scale = "free"
)

SINGLE_PATHWAY_PLOT_BASE_SIZE <- 7.2
SINGLE_PATHWAY_PLOT_SIZE_PER_GENE <- 0.090
SINGLE_PATHWAY_PLOT_SIZE_PER_SAMPLE <- 0.115
SINGLE_PATHWAY_PLOT_SIZE_PER_TITLE_LINE <- 0.32
SINGLE_PATHWAY_PLOT_MIN_SIZE <- 7.2
SINGLE_PATHWAY_PLOT_MAX_SIZE <- 34.0
SINGLE_PATHWAY_GENE_TEXT_MIN_SIZE <- 2.2
SINGLE_PATHWAY_GENE_TEXT_MAX_SIZE <- 7.0
SINGLE_PATHWAY_GENE_TEXT_WIDTH_FACTOR <- 0.62
SINGLE_PATHWAY_SAMPLE_TEXT_MIN_SIZE <- 2.6
SINGLE_PATHWAY_SAMPLE_TEXT_MAX_SIZE <- 6.6
SINGLE_PATHWAY_SAMPLE_TEXT_HEIGHT_FACTOR <- 0.74

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 1. 加载包、公共函数和内部声明 ------------------------------------------------

required_packages <- c(
  "SummarizedExperiment",
  "clusterProfiler",
  "msigdbr",
  "GseaVis",
  "ggplot2",
  "qs2",
  "parallel"
)

is_package_available <- function(package_name) {
  suppressWarnings(
    suppressPackageStartupMessages(
      requireNamespace(package_name, quietly = TRUE)
    )
  )
}

missing_packages <- required_packages[
  !vapply(required_packages, is_package_available, logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Please install required R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

source(FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(REPORT_TABLE_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(GSEA_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

READABLE_GENE_SYMBOLS <- TRUE
USE_QS2_CACHE <- TRUE

# 07号默认复用06号生成的qs2缓存对象。若缓存不存在，脚本会按相同配置补算对象。
REFRESH_QS2_CACHE <- Sys.getenv("GSEA_REFRESH_QS2_CACHE", unset = "0") == "1"
QS2_CACHE_DIR <- file.path("temporary", DATA_TYPE, DATASET_ID, "GSEA_qs2_cache")
MSIGDB_REFERENCE_DIR <- file.path("data", "reference", "msigdb")
MSIGDB_REFERENCE_MAX_AGE_DAYS <- 7

PARALLEL_WORKERS <- get_available_worker_count()
configure_parallel_runtime(
  task_workers = PARALLEL_WORKERS,
  inner_workers = PARALLEL_WORKERS,
  qs2_threads = PARALLEL_WORKERS
)


# 2. 准备输入文件、基因集和表达矩阵 ------------------------------------------

stopifnot(GENE_ID_TYPE %in% names(GENE_ID_COLUMN_BY_TYPE))
stopifnot(GENE_ID_TYPE %in% names(MSIGDB_GENE_COLUMN_BY_TYPE))

MSIGDB_GENESET_CATALOG <- build_msigdb_geneset_catalog()
RUNTIME_GENESETS_TO_RUN <- get_runtime_genesets_to_run()
GSEA_GENESET_CONFIG <- select_msigdb_genesets(
  catalog = MSIGDB_GENESET_CATALOG,
  genesets_to_run = RUNTIME_GENESETS_TO_RUN
)

file_info <- get_deg_file_info(TABLE_ROOT)
selected_analyses <- get_selected_analysis_names(
  file_info = file_info,
  analyses_to_plot = ANALYSES_TO_RUN
)

if (CLEAN_GSEA_PLOT_OUTPUT_DIR && dir.exists(PLOT_ROOT)) {
  unlink(PLOT_ROOT, recursive = TRUE)
}
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))
stopifnot(TPM_ASSAY_NAME %in% names(SummarizedExperiment::assays(se)))

analysis_designs <- get_analysis_designs(clinical_data)
expression_matrix_all <- as.matrix(SummarizedExperiment::assay(se, TPM_ASSAY_NAME))
expression_matrix_all <- log2(expression_matrix_all + 1)
gene_annotation_all <- as.data.frame(
  SummarizedExperiment::rowData(se),
  stringsAsFactors = FALSE
)

sample_info_all <- clinical_data[
  match(colnames(expression_matrix_all), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
rownames(sample_info_all) <- sample_info_all$Sample_ID
stopifnot(all(sample_info_all$Sample_ID == colnames(expression_matrix_all)))

cat("\nGSEA plot runtime configuration:\n")
cat("Selected analyses: ", length(selected_analyses), "\n", sep = "")
cat(
  "Selected MSigDB gene set categories: ",
  length(GSEA_GENESET_CONFIG),
  if (length(RUNTIME_GENESETS_TO_RUN) == 1 &&
      tolower(RUNTIME_GENESETS_TO_RUN) == "all") {
    " (all available categories from msigdbr)"
  } else {
    ""
  },
  "\n",
  sep = ""
)
cat("Available workers: ", PARALLEL_WORKERS, "\n", sep = "")
cat("Output root:  ", OUTPUT_ROOT, "\n", sep = "")
cat("Refresh qs2 cache: ", REFRESH_QS2_CACHE, "\n", sep = "")

cat("\nLoading MSigDB gene sets...\n")
geneset_cache <- lapply(names(GSEA_GENESET_CONFIG), function(geneset_name) {
  config <- GSEA_GENESET_CONFIG[[geneset_name]]
  terms <- load_msigdb_terms(geneset_name, config)

  list(
    config = config,
    term2gene = terms$term2gene,
    term2name = terms$term2name,
    cache_source = terms$Cache_Source
  )
})
names(geneset_cache) <- names(GSEA_GENESET_CONFIG)

cat("\nPreparing analysis-level plotting inputs...\n")
analysis_cache <- lapply(selected_analyses, function(analysis_name) {
  deg_index <- match(analysis_name, file_info$Analysis_Name)
  deg_file <- file_info$All_Genes_File[deg_index]
  deg_result <- read_deg_result(file_info, analysis_name)
  gene_list <- prepare_gene_list(
    deg_table = deg_result,
    gene_id_type = GENE_ID_TYPE,
    rank_metric_column = RANK_METRIC_COLUMN
  )
  expression_data <- prepare_single_pathway_expression_table(
    analysis_name = analysis_name,
    expression_matrix_all = expression_matrix_all,
    gene_annotation_all = gene_annotation_all,
    analysis_designs = analysis_designs,
    sample_info_all = sample_info_all
  )

  list(
    deg_file = deg_file,
    gene_list = gene_list,
    expression_data = expression_data
  )
})
names(analysis_cache) <- selected_analyses

task_table <- expand.grid(
  Analysis_Name = selected_analyses,
  GeneSet_Name = names(geneset_cache),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_tasks <- nrow(task_table)

compute_summary_file <- resolve_report_csv_file(
  file.path(OUTPUT_ROOT, "tables", "GSEA_summary", "summary.csv")
)
compute_summary_table <- if (file.exists(compute_summary_file)) {
  read.csv(compute_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame()
}

parallel_strategy <- setup_parallel_strategy(
  total_tasks = total_tasks,
  max_workers = PARALLEL_WORKERS,
  inner_label = "GSEA object nproc per task",
  nested_label = "Single-pathway workers"
)
GSEA_TASK_WORKERS <- parallel_strategy$task_workers
GSEA_INNER_NPROC <- parallel_strategy$inner_workers
SINGLE_PATHWAY_PLOT_WORKERS <- parallel_strategy$nested_workers


# 3. 批量绘制GSEA结果 ----------------------------------------------------------

run_gsea_plot_task <- function(task_id) {
  analysis_name <- task_table$Analysis_Name[task_id]
  geneset_name <- task_table$GeneSet_Name[task_id]
  analysis_input <- analysis_cache[[analysis_name]]
  cache <- geneset_cache[[geneset_name]]
  output_name <- cache$config$output_name
  output_dir_name <- sanitize_file_name(output_name)

  table_output_dir <- file.path(
    TABLE_OUTPUT_ROOT,
    analysis_name,
    "GSEA",
    output_dir_name
  )
  plot_output_dir <- file.path(
    PLOT_ROOT,
    sanitize_file_name(analysis_name),
    output_dir_name
  )
  dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

  csv_file <- file.path(table_output_dir, "gsea_result.csv")
  csv_file <- resolve_report_csv_file(csv_file)
  if (!file.exists(csv_file)) {
    stop("GSEA result table was not found. Please run 06_gsea_analysis.R first: ", csv_file)
  }

  result_table <- read.csv(
    csv_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  gsea_run <- load_or_run_gsea(
    analysis_name = analysis_name,
    deg_file = analysis_input$deg_file,
    geneset_name = geneset_name,
    config = cache$config,
    gene_list = analysis_input$gene_list,
    term2gene = cache$term2gene,
    term2name = cache$term2name
  )

  gsea_result <- gsea_run$result
  summary_source <- gsea_run$source
  if (nrow(compute_summary_table) > 0 &&
      all(c("Analysis_Name", "GeneSet_Name", "Source") %in% colnames(compute_summary_table))) {
    source_index <- which(
      compute_summary_table$Analysis_Name == analysis_name &
        compute_summary_table$GeneSet_Name == geneset_name
    )
    if (length(source_index) > 0) {
      summary_source <- compute_summary_table$Source[source_index[1]]
    }
  }

  if (is_gsea_result_object(gsea_result)) {
    plot_gsea_result <- prepare_gsea_result_for_plot(gsea_result)
    plot_result_table <- as.data.frame(plot_gsea_result)
    dotplot_result <- make_gsea_dotplot(
      gsea_result = plot_gsea_result,
      result_table = plot_result_table,
      analysis_name = analysis_name,
      geneset_name = geneset_name
    )
  } else {
    dotplot_result <- list(
      plot = make_empty_gsea_plot(),
      shown_terms = 0L,
      plot_labels = character(0)
    )
  }

  plot_size <- get_gsea_dotplot_size(
    shown_terms = dotplot_result$shown_terms,
    plot_labels = dotplot_result$plot_labels
  )
  pdf_file <- file.path(plot_output_dir, "dotplot.pdf")
  plot_files <- with_gsea_warnings_suppressed(
    save_ggplot_pdf_png(
      plot = dotplot_result$plot,
      pdf_file = pdf_file,
      width = plot_size$width,
      height = plot_size$height
    )
  )

  single_pathway_count <- if (DRAW_SINGLE_PATHWAY_GSEA &&
                              is_gsea_result_object(gsea_result)) {
    with_gsea_warnings_suppressed(
      save_single_pathway_gsea_plots(
        gsea_result = gsea_result,
        result_table = result_table,
        analysis_name = analysis_name,
        geneset_name = geneset_name,
        plot_output_dir = plot_output_dir,
        expression_data = analysis_input$expression_data
      )
    )
  } else {
    0L
  }

  data.frame(
    Analysis_Name = analysis_name,
    GeneSet_Name = geneset_name,
    Source = summary_source,
    Ranked_Genes = length(analysis_input$gene_list),
    GSEA_Terms = nrow(result_table),
    Positive_NES = count_nes_direction(result_table, "positive"),
    Negative_NES = count_nes_direction(result_table, "negative"),
    Single_Pathway_Plots = single_pathway_count,
    CSV_File = csv_file,
    PDF_File = plot_files$pdf_file,
    PNG_File = plot_files$png_file,
    stringsAsFactors = FALSE
  )
}

cat("\nRunning batch GSEA plotting...\n")
task_ids <- seq_len(total_tasks)
summary_records <- run_parallel_tasks_with_progress(
  task_ids = task_ids,
  task_function = run_gsea_plot_task,
  workers = GSEA_TASK_WORKERS
)
stop_on_parallel_errors(summary_records, task_ids = task_ids, label = "GSEA plotting tasks")


# 4. 终端快速汇总 --------------------------------------------------------------

summary_table <- do.call(rbind, summary_records)
rownames(summary_table) <- NULL

summary_output_dir <- file.path(OUTPUT_ROOT, "tables", "GSEA_summary")
dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
unlink(
  list.files(
    summary_output_dir,
    pattern = "[.](csv|md|tex)$",
    full.names = TRUE
  ),
  force = TRUE
)
summary_csv_file <- write_csv_with_report_previews(
  summary_table,
  file.path(summary_output_dir, "summary.csv"),
  n_rows = 21
)

cat("\nGSEA plot summary:\n")
print(
  summary_table[
    ,
    c(
      "Analysis_Name", "GeneSet_Name", "Source", "Ranked_Genes",
      "GSEA_Terms", "Positive_NES", "Negative_NES",
      "Single_Pathway_Plots"
    )
  ],
  row.names = FALSE
)

cat("\nOutput summary:\n")
cat("GSEA plot directory:  ", file.path(PLOT_ROOT, "<analysis_name>"), "\n", sep = "")
cat("GSEA summary table:   ", summary_csv_file, "\n", sep = "")
cat("PDF/PNG dotplots: ", nrow(summary_table), " each\n", sep = "")
cat("Single-pathway GSEA PDF/PNG plots: ", sum(summary_table$Single_Pathway_Plots), " each\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")

cat("\nBatch GSEA plotting finished.\n")
