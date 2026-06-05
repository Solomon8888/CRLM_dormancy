# GSE114012批量GSEA运算
#
# 自动读取01号脚本输出的每套DEG/all_genes.csv，
# 使用clusterProfiler::GSEA和msigdbr基因集批量运行GSEA分析。
# 本脚本只负责GSEA运算和CSV/MD/TEX表格保存；GSEA图片由07号脚本生成。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
GSEA_FUNCTION_FILE <- "scripts/functions/gsea_common_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")

# 测试时可临时设置GSEA_OUTPUT_ROOT，把结果输出到temporary而不覆盖正式results。
OUTPUT_ROOT <- Sys.getenv("GSEA_OUTPUT_ROOT", unset = RESULT_ROOT)
TABLE_OUTPUT_ROOT <- file.path(OUTPUT_ROOT, "tables")
PLOT_ROOT <- file.path(OUTPUT_ROOT, "plots", "GSEA")

# 需要运行哪些差异分析设计。设为"all"时自动运行全部DEG/all_genes.csv。
ANALYSES_TO_RUN <- "all"
# ANALYSES_TO_RUN <- c("DLD1", "DLD1_HCT15_SW48")

# 物种配置。常用可选项："human"、"Mus musculus"、"Rattus norvegicus"。
SPECIES <- "human"

# GSEA排序基因列表配置。
GENE_ID_TYPE <- "ENTREZ"  # 可选："ENTREZ", "SYMBOL", "ENSEMBL"
RANK_METRIC_COLUMN <- "t"

# GSEA显著性阈值配置。
# GSEA_SIGNIFICANCE_COLUMN用于记录和绘图筛选时采用的p值类型；
# 目前与clusterProfiler::GSEA的pAdjustMethod="BH"配套，统一采用p.adjust阈值。
GSEA_SIGNIFICANCE_COLUMN <- "pvalue"
GSEA_SIGNIFICANCE_CUTOFF <- 0.05

# clusterProfiler::GSEA官方运算参数。
GSEA_PARAMS <- list(
  exponent = 1,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = GSEA_SIGNIFICANCE_CUTOFF,
  pAdjustMethod = "BH",
  verbose = TRUE,
  nPerm = 1000,
  method = "multilevel",
  adaptive = FALSE,
  minPerm = 101,
  maxPerm = 1e5,
  pvalThreshold = 0.1
)

# 当前批量运行哪些MSigDB基因集。
# 键名格式来自msigdbr::msigdbr_collections()中的gs_collection与gs_subcollection。
# 例如Hallmark为"H"；BioCarta为"C2:CP:BIOCARTA"；
# GO Biological Process为"C5:GO:BP"；TFT Legacy为"C3:TFT:TFT_LEGACY"。
# 设为"all"时，会自动运行msigdbr当前数据库里全部可用基因集类别。
GSEA_GENESETS_TO_RUN <- c(
  "H",
  "C2:CP:BIOCARTA",
  "C2:CP:KEGG_MEDICUS",
  "C2:CP:KEGG_LEGACY",
  "C2:CP:REACTOME",
  "C2:CP:WIKIPATHWAYS",
  "C3:TFT:TFT_LEGACY",
  "C3:TFT:GTRD",
  "C5:GO:BP",
  "C5:GO:CC",
  "C5:GO:MF",
  "C5:HPO",
  "C6",
  "C7:IMMUNESIGDB"
)
# GSEA_GENESETS_TO_RUN <- "all"

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 1. 加载包、公共函数和内部声明 ------------------------------------------------

required_packages <- c(
  "clusterProfiler",
  "msigdbr",
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

# 以下为脚本内部固定声明，通常不需要在日常分析中修改。
CLEAN_GSEA_OUTPUT_DIR <- TRUE
READABLE_GENE_SYMBOLS <- TRUE
USE_QS2_CACHE <- TRUE

# 默认重新计算核心GSEA缓存，适合每次调整参数后重新运行。
# 若确认参数完全不变且只想复用缓存，可在终端临时设置：
# GSEA_REFRESH_QS2_CACHE=0 Rscript scripts/GSE114012/06_gsea_analysis.R
REFRESH_QS2_CACHE <- Sys.getenv("GSEA_REFRESH_QS2_CACHE", unset = "1") == "1"
QS2_CACHE_DIR <- file.path("temporary", DATA_TYPE, DATASET_ID, "GSEA_qs2_cache")
MSIGDB_REFERENCE_DIR <- file.path("data", "reference", "msigdb")
MSIGDB_REFERENCE_MAX_AGE_DAYS <- 7

PARALLEL_WORKERS <- get_available_worker_count()
configure_parallel_runtime(
  task_workers = PARALLEL_WORKERS,
  inner_workers = PARALLEL_WORKERS,
  qs2_threads = PARALLEL_WORKERS
)


# 2. 准备输入文件和基因集 ------------------------------------------------------

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

clean_previous_gsea_outputs(selected_analyses)

cat("\nGSEA runtime configuration:\n")
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
cat("Previous GSEA outputs were cleaned before this run.\n")

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

geneset_summary <- do.call(
  rbind,
  lapply(names(geneset_cache), function(geneset_name) {
    cache <- geneset_cache[[geneset_name]]
    data.frame(
      GeneSet_Name = geneset_name,
      Output_Name = cache$config$output_name,
      Terms = length(unique(cache$term2gene$term)),
      Term_Gene_Links = nrow(cache$term2gene),
      Source = cache$cache_source,
      stringsAsFactors = FALSE
    )
  })
)
print(geneset_summary, row.names = FALSE)

cat("\nPreparing analysis-level inputs...\n")
analysis_cache <- lapply(selected_analyses, function(analysis_name) {
  deg_index <- match(analysis_name, file_info$Analysis_Name)
  deg_file <- file_info$All_Genes_File[deg_index]
  deg_result <- read_deg_result(file_info, analysis_name)
  gene_list <- prepare_gene_list(
    deg_table = deg_result,
    gene_id_type = GENE_ID_TYPE,
    rank_metric_column = RANK_METRIC_COLUMN
  )

  list(
    deg_file = deg_file,
    gene_list = gene_list
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

parallel_strategy <- setup_parallel_strategy(
  total_tasks = total_tasks,
  max_workers = PARALLEL_WORKERS,
  inner_label = "GSEA nproc per task",
  nested_label = "Nested workers"
)
GSEA_TASK_WORKERS <- parallel_strategy$task_workers
GSEA_INNER_NPROC <- parallel_strategy$inner_workers


# 3. 批量运行GSEA并保存表格 ----------------------------------------------------

run_gsea_compute_task <- function(task_id) {
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
  dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

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
  csv_file <- file.path(table_output_dir, "gsea_result.csv")
  result_table <- write_gsea_result_tables(gsea_result, csv_file)
  csv_file <- resolve_report_csv_file(csv_file)

  data.frame(
    Analysis_Name = analysis_name,
    GeneSet_Name = geneset_name,
    Source = gsea_run$source,
    Ranked_Genes = length(analysis_input$gene_list),
    GSEA_Terms = nrow(result_table),
    Positive_NES = count_nes_direction(result_table, "positive"),
    Negative_NES = count_nes_direction(result_table, "negative"),
    Single_Pathway_Plots = 0L,
    CSV_File = csv_file,
    PDF_File = "",
    PNG_File = "",
    stringsAsFactors = FALSE
  )
}

cat("\nRunning batch GSEA analyses...\n")
task_ids <- seq_len(total_tasks)
summary_records <- run_parallel_tasks_with_progress(
  task_ids = task_ids,
  task_function = run_gsea_compute_task,
  workers = GSEA_TASK_WORKERS
)
stop_on_parallel_errors(summary_records, task_ids = task_ids, label = "GSEA tasks")


# 4. 终端快速汇总 --------------------------------------------------------------

summary_table <- do.call(rbind, summary_records)
rownames(summary_table) <- NULL

summary_output_dir <- file.path(OUTPUT_ROOT, "tables", "GSEA_summary")
dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
summary_csv_file <- write_csv_with_report_previews(
  summary_table,
  file.path(summary_output_dir, "summary.csv"),
  n_rows = 21
)

cat("\nGSEA compute summary:\n")
print(
  summary_table[
    ,
    c(
      "Analysis_Name", "GeneSet_Name", "Source", "Ranked_Genes",
      "GSEA_Terms", "Positive_NES", "Negative_NES"
    )
  ],
  row.names = FALSE
)

cat("\nOutput summary:\n")
cat("GSEA table directory: ", file.path(TABLE_OUTPUT_ROOT, "<analysis_name>", "GSEA"), "\n", sep = "")
cat("GSEA summary table:   ", summary_csv_file, "\n", sep = "")
cat("CSV/MD/TEX result sets: ", nrow(summary_table), " each\n", sep = "")
cat("GSEA plots are generated by scripts/GSE114012/07_gsea_plot.R\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")

cat("\nBatch GSEA analysis finished.\n")
