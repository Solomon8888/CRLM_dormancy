# GSE50760 differential-expression GSEA.
#
# Reads all_genes.csv files from script 02 and runs MSigDB GSEA using the limma
# t statistic as the ranking metric.


# 0. Small env parsers ---------------------------------------------------------

parse_env_vector <- function(variable_name, default_value) {
  value <- trimws(Sys.getenv(variable_name, unset = ""))
  if (value == "") {
    return(default_value)
  }

  unique(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

parse_env_logical <- function(variable_name, default_value) {
  value <- tolower(trimws(Sys.getenv(variable_name, unset = "")))
  if (value == "") {
    return(default_value)
  }

  value %in% c("1", "true", "t", "yes", "y")
}

parse_env_integer <- function(variable_name, default_value) {
  value <- trimws(Sys.getenv(variable_name, unset = ""))
  if (value == "") {
    return(default_value)
  }

  as.integer(value)
}


# 1. Config -------------------------------------------------------------------

DATASET_ID <- "GSE50760"
DATA_TYPE <- "ngs"

DEG_SCRIPT <- "scripts/GSE50760/02_limma_differential_expression.R"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
GSEA_FUNCTION_FILE <- "scripts/functions/gsea_common_functions.R"
BATCH_VIS_FUNCTION_FILE <- "scripts/functions/gse_batch_visualization_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")
OUTPUT_ROOT <- Sys.getenv("GSE50760_GSEA_OUTPUT_ROOT", unset = RESULT_ROOT)
TABLE_OUTPUT_ROOT <- file.path(OUTPUT_ROOT, "tables")
PLOT_ROOT <- file.path(OUTPUT_ROOT, "plots", "GSEA")

ANALYSES_TO_RUN <- parse_env_vector("GSE50760_DE_GSEA_ANALYSES", "all")
if (length(ANALYSES_TO_RUN) == 1 && tolower(ANALYSES_TO_RUN) == "all") {
  ANALYSES_TO_RUN <- "all"
}

SPECIES <- Sys.getenv("GSE50760_DE_GSEA_SPECIES", unset = "human")
GENE_ID_TYPE <- Sys.getenv("GSE50760_DE_GSEA_GENE_ID_TYPE", unset = "ENTREZ")
RANK_METRIC_COLUMN <- Sys.getenv("GSE50760_DE_GSEA_RANK_COLUMN", unset = "t")

GSEA_SIGNIFICANCE_COLUMN <- Sys.getenv("GSE50760_DE_GSEA_P_COLUMN", unset = "pvalue")
GSEA_SIGNIFICANCE_CUTOFF <- as.numeric(Sys.getenv("GSE50760_DE_GSEA_P_CUTOFF", unset = "0.05"))
GSEA_PARAMS <- list(
  exponent = 1,
  minGSSize = parse_env_integer("GSE50760_DE_GSEA_MIN_GS_SIZE", 5L),
  maxGSSize = parse_env_integer("GSE50760_DE_GSEA_MAX_GS_SIZE", 500L),
  pvalueCutoff = GSEA_SIGNIFICANCE_CUTOFF,
  pAdjustMethod = "BH",
  verbose = parse_env_logical("GSE50760_DE_GSEA_VERBOSE", FALSE),
  nPerm = 1000,
  method = "multilevel",
  adaptive = FALSE,
  minPerm = 101,
  maxPerm = 1e5,
  pvalThreshold = 0.1
)

GSEA_GENESETS_TO_RUN <- parse_env_vector(
  "GSE50760_DE_GSEA_GENESETS",
  c(
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
)

READABLE_GENE_SYMBOLS <- parse_env_logical("GSE50760_DE_GSEA_READABLE_SYMBOLS", TRUE)
USE_QS2_CACHE <- parse_env_logical("GSE50760_DE_GSEA_USE_QS2_CACHE", TRUE)
REFRESH_QS2_CACHE <- parse_env_logical("GSE50760_DE_GSEA_REFRESH_QS2_CACHE", TRUE)
CLEAN_GSEA_OUTPUT_DIR <- parse_env_logical("GSE50760_DE_GSEA_CLEAN", TRUE)

QS2_CACHE_DIR <- file.path("temporary", DATA_TYPE, DATASET_ID, "GSEA_qs2_cache")
MSIGDB_REFERENCE_DIR <- file.path("data", "reference", "msigdb")
MSIGDB_REFERENCE_MAX_AGE_DAYS <- parse_env_integer("GSE50760_DE_GSEA_MSIGDB_MAX_AGE_DAYS", 7L)

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 2. Packages and helpers ------------------------------------------------------

required_packages <- c("clusterProfiler", "msigdbr", "qs2", "parallel")

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
source(BATCH_VIS_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

PARALLEL_WORKERS <- get_available_worker_count()
configure_parallel_runtime(
  task_workers = PARALLEL_WORKERS,
  inner_workers = PARALLEL_WORKERS,
  qs2_threads = PARALLEL_WORKERS
)


# 3. Prepare inputs and gene sets ---------------------------------------------

deg_file_candidates <- c(
  file.path(list.dirs(TABLE_ROOT, recursive = FALSE, full.names = TRUE), "DEG", "all_genes.csv"),
  file.path(list.dirs(TABLE_ROOT, recursive = FALSE, full.names = TRUE), "DEG", "csv", "all_genes.csv")
)
deg_file_candidates <- deg_file_candidates[file.exists(deg_file_candidates)]

if (length(deg_file_candidates) == 0) {
  source(DEG_SCRIPT)
}

stopifnot(GENE_ID_TYPE %in% names(GENE_ID_COLUMN_BY_TYPE))
stopifnot(GENE_ID_TYPE %in% names(MSIGDB_GENE_COLUMN_BY_TYPE))

MSIGDB_GENESET_CATALOG <- build_msigdb_geneset_catalog()
RUNTIME_GENESETS_TO_RUN <- get_runtime_genesets_to_run()
GSEA_GENESET_CONFIG <- select_msigdb_genesets(
  catalog = MSIGDB_GENESET_CATALOG,
  genesets_to_run = RUNTIME_GENESETS_TO_RUN
)

file_info <- collect_batch_deg_file_info(TABLE_ROOT)
file_info <- file_info[file_info$Plot_Category == "Main_DE", , drop = FALSE]
file_info <- file_info[, c("Analysis_Name", "All_Genes_File"), drop = FALSE]
selected_analyses <- get_selected_analysis_names(
  file_info = file_info,
  analyses_to_plot = ANALYSES_TO_RUN
)

clean_previous_gsea_outputs(selected_analyses)

cat("\nGSE50760 DE-GSEA runtime configuration:\n")
cat("Selected analyses: ", length(selected_analyses), "\n", sep = "")
cat("Selected MSigDB gene set categories: ", length(GSEA_GENESET_CONFIG), "\n", sep = "")
cat("Available workers: ", PARALLEL_WORKERS, "\n", sep = "")
cat("Output root: ", OUTPUT_ROOT, "\n", sep = "")
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

cat("\nPreparing analysis-level GSEA inputs...\n")
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

parallel_strategy <- setup_parallel_strategy(
  total_tasks = nrow(task_table),
  max_workers = PARALLEL_WORKERS,
  inner_label = "GSEA nproc per task",
  nested_label = "Nested workers"
)
GSEA_INNER_NPROC <- parallel_strategy$inner_workers


# 4. Run GSEA tasks ------------------------------------------------------------

run_gsea_compute_task <- function(task_id) {
  analysis_name <- task_table$Analysis_Name[task_id]
  geneset_name <- task_table$GeneSet_Name[task_id]
  analysis_input <- analysis_cache[[analysis_name]]
  cache <- geneset_cache[[geneset_name]]
  output_name <- cache$config$output_name
  output_dir_name <- sanitize_file_name(output_name)

  table_output_dir <- file.path(TABLE_OUTPUT_ROOT, analysis_name, "GSEA", output_dir_name)
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
    CSV_File = csv_file,
    stringsAsFactors = FALSE
  )
}

cat("\nRunning GSE50760 batch DE-GSEA analyses...\n")
task_ids <- seq_len(nrow(task_table))
summary_records <- run_parallel_tasks_with_progress(
  task_ids = task_ids,
  task_function = run_gsea_compute_task,
  workers = parallel_strategy$task_workers
)
stop_on_parallel_errors(summary_records, task_ids = task_ids, label = "GSE50760 DE-GSEA tasks")


# 5. Save summary --------------------------------------------------------------

summary_table <- do.call(rbind, summary_records)
rownames(summary_table) <- NULL

summary_output_dir <- file.path(OUTPUT_ROOT, "tables", "GSEA_summary")
dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
summary_csv_file <- write_csv_with_report_previews(
  summary_table,
  file.path(summary_output_dir, "summary.csv"),
  n_rows = 21
)

cat("\nGSE50760 DE-GSEA compute summary:\n")
print(
  summary_table[
    ,
    c("Analysis_Name", "GeneSet_Name", "Source", "Ranked_Genes", "GSEA_Terms", "Positive_NES", "Negative_NES")
  ],
  row.names = FALSE
)

cat("\nOutput summary:\n")
cat("GSEA table directory: ", file.path(TABLE_OUTPUT_ROOT, "<analysis_name>", "GSEA"), "\n", sep = "")
cat("GSEA summary table:   ", summary_csv_file, "\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")

cat("\nGSE50760 batch DE-GSEA analysis finished.\n")
