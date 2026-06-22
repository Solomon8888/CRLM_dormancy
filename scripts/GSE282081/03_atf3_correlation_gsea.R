# GSE282081 ATF3相关性功能分析
#
# 将ATF3表达作为连续变量，在不同样本范围内计算全基因Spearman相关性，
# 再使用相关系数作为排序统计量运行MSigDB GSEA，用于观察ATF3高表达相关
# 的功能通路。当前GSE282081本地样本是SW480体外肝细胞/J2共培养模型，
# 因此这里分析的是liver-niche-like模型中的ATF3功能关联。


# 0. 小型环境变量解析函数 ------------------------------------------------------

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


# 1. 可修改配置 ----------------------------------------------------------------

DATASET_ID <- "GSE282081"
DATA_TYPE <- "ngs"
TARGET_GENE <- Sys.getenv("GSE282081_ATF3_TARGET_GENE", unset = "ATF3")

SE_RDS_FILE <- "data/ngs/GSE282081/data_prepare/GSE282081_se_raw.rds"
CLINICAL_EDIT_FILE <- "data/ngs/GSE282081/data_prepare/GSE282081_clinical_edit.csv"
SAMPLE_DESIGN_SCRIPT <- "scripts/GSE282081/00_prepare_sample_design.R"

LIMMA_FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
GSEA_FUNCTION_FILE <- "scripts/functions/gsea_common_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables", "ATF3_function")

TARGET_EXPRESSION_ASSAY <- Sys.getenv("GSE282081_ATF3_ASSAY", unset = "tpm")
TARGET_EXPRESSION_LOG2 <- parse_env_logical("GSE282081_ATF3_LOG2", TRUE)
GENE_BIOTYPE_FILTER <- "coding"
CORRELATION_METHOD <- "spearman"
MIN_SAMPLES_FOR_CORRELATION <- parse_env_integer("GSE282081_ATF3_MIN_N", 6L)
CORRELATION_P_CUTOFF <- as.numeric(Sys.getenv("GSE282081_ATF3_COR_P_CUTOFF", unset = "0.05"))

CORRELATION_SCOPES_TO_RUN <- parse_env_vector(
  "GSE282081_ATF3_CORRELATION_SCOPES",
  c(
    "All_samples",
    "Hepatocyte_coculture",
    "No_hepatocyte_coculture",
    "2D_culture",
    "3D_organoid"
  )
)

RUN_GSEA <- parse_env_logical("GSE282081_ATF3_RUN_GSEA", TRUE)

SPECIES <- Sys.getenv("GSE282081_ATF3_SPECIES", unset = "human")
GENE_ID_TYPE <- Sys.getenv("GSE282081_ATF3_GENE_ID_TYPE", unset = "ENTREZ")
RANK_METRIC_COLUMN <- "Correlation"

GSEA_SIGNIFICANCE_COLUMN <- Sys.getenv("GSE282081_ATF3_GSEA_P_COLUMN", unset = "pvalue")
GSEA_SIGNIFICANCE_CUTOFF <- as.numeric(Sys.getenv("GSE282081_ATF3_GSEA_P_CUTOFF", unset = "0.05"))
GSEA_PARAMS <- list(
  exponent = 1,
  minGSSize = parse_env_integer("GSE282081_ATF3_GSEA_MIN_GS_SIZE", 5L),
  maxGSSize = parse_env_integer("GSE282081_ATF3_GSEA_MAX_GS_SIZE", 500L),
  pvalueCutoff = GSEA_SIGNIFICANCE_CUTOFF,
  pAdjustMethod = "BH",
  verbose = parse_env_logical("GSE282081_ATF3_GSEA_VERBOSE", FALSE),
  nPerm = 1000,
  method = "multilevel",
  adaptive = FALSE,
  minPerm = 101,
  maxPerm = 1e5,
  pvalThreshold = 0.1
)

GSEA_GENESETS_TO_RUN <- parse_env_vector(
  "GSE282081_ATF3_GSEA_GENESETS",
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

READABLE_GENE_SYMBOLS <- parse_env_logical("GSE282081_ATF3_READABLE_SYMBOLS", TRUE)
USE_QS2_CACHE <- parse_env_logical("GSE282081_ATF3_USE_QS2_CACHE", TRUE)
REFRESH_QS2_CACHE <- parse_env_logical("GSE282081_ATF3_REFRESH_QS2_CACHE", TRUE)
CLEAN_GSEA_OUTPUT_DIR <- parse_env_logical("GSE282081_ATF3_CLEAN_GSEA", TRUE)

QS2_CACHE_DIR <- file.path("temporary", DATA_TYPE, DATASET_ID, "ATF3_GSEA_qs2_cache")
MSIGDB_REFERENCE_DIR <- file.path("data", "reference", "msigdb")
MSIGDB_REFERENCE_MAX_AGE_DAYS <- parse_env_integer("GSE282081_ATF3_MSIGDB_MAX_AGE_DAYS", 7L)

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 2. 加载包和函数 --------------------------------------------------------------

required_packages <- c(
  "SummarizedExperiment",
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

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

source(LIMMA_FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(REPORT_TABLE_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(GSEA_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

PARALLEL_WORKERS <- get_available_worker_count()
configure_parallel_runtime(
  task_workers = PARALLEL_WORKERS,
  inner_workers = PARALLEL_WORKERS,
  qs2_threads = PARALLEL_WORKERS
)


# 3. 读取表达和样本信息 --------------------------------------------------------

if (!file.exists(CLINICAL_EDIT_FILE)) {
  source(SAMPLE_DESIGN_SCRIPT)
}

dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_EDIT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot(TARGET_EXPRESSION_ASSAY %in% names(assays(se)))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot(all(colnames(se) %in% clinical_data$Sample_ID))

clinical_data <- clinical_data[
  match(colnames(se), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
stopifnot(all(clinical_data$Sample_ID == colnames(se)))

exprSet_all <- as.matrix(assay(se, TARGET_EXPRESSION_ASSAY))
if (TARGET_EXPRESSION_LOG2) {
  exprSet_all <- log2(exprSet_all + 1)
}

feature_id <- rownames(exprSet_all)
if (is.null(feature_id)) {
  feature_id <- paste0("Feature_", seq_len(nrow(exprSet_all)))
  rownames(exprSet_all) <- feature_id
}

gene_annotation <- data.frame(
  Feature_ID = feature_id,
  as.data.frame(rowData(se), stringsAsFactors = FALSE),
  check.names = FALSE
)
rownames(gene_annotation) <- rownames(exprSet_all)

target_index_all <- which(toupper(trimws(as.character(gene_annotation$Symbol))) == toupper(TARGET_GENE))
stopifnot(length(target_index_all) >= 1)
target_index_all <- target_index_all[1]
target_feature_id <- rownames(exprSet_all)[target_index_all]
target_expression_all <- as.numeric(exprSet_all[target_index_all, ])

gene_filter <- filter_genes_by_biotype(
  exprSet = exprSet_all,
  gene_annotation = gene_annotation,
  biotype_filter = GENE_BIOTYPE_FILTER
)

exprSet <- gene_filter$exprSet
gene_annotation <- gene_filter$gene_annotation


# 4. 样本范围和相关性计算 ------------------------------------------------------

scope_definitions <- list(
  All_samples = function(sample_info) rep(TRUE, nrow(sample_info)),
  Hepatocyte_coculture = function(sample_info) {
    sample_info$hepatocyte_coculture == "Hepatocyte_coculture"
  },
  No_hepatocyte_coculture = function(sample_info) {
    sample_info$hepatocyte_coculture == "No_hepatocyte_coculture"
  },
  `2D_culture` = function(sample_info) {
    sample_info$culture_model == "2D_culture"
  },
  `3D_organoid` = function(sample_info) {
    sample_info$culture_model == "3D_organoid"
  }
)

missing_scopes <- setdiff(CORRELATION_SCOPES_TO_RUN, names(scope_definitions))
if (length(missing_scopes) > 0) {
  stop(
    "Undefined correlation scopes: ",
    paste(missing_scopes, collapse = ", ")
  )
}

compute_correlation_p_value <- function(rho, n) {
  p_value <- rep(NA_real_, length(rho))
  valid <- is.finite(rho) & n > 2
  rho_for_p <- pmin(pmax(rho[valid], -1 + 1e-15), 1 - 1e-15)
  t_stat <- rho_for_p * sqrt((n - 2) / pmax(1 - rho_for_p^2, .Machine$double.eps))
  p_value[valid] <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
  p_value
}

compute_scope_correlation <- function(scope_name) {
  keep_sample <- scope_definitions[[scope_name]](clinical_data)
  keep_sample[is.na(keep_sample)] <- FALSE
  sample_ids <- clinical_data$Sample_ID[keep_sample]
  sample_count <- length(sample_ids)

  output_dir <- file.path(TABLE_ROOT, scope_name, "correlation")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (sample_count < MIN_SAMPLES_FOR_CORRELATION) {
    summary_table <- data.frame(
      Dataset = DATASET_ID,
      Analysis_Name = scope_name,
      Target_Gene = TARGET_GENE,
      Status = "Skipped: too few samples",
      Samples_Used = sample_count,
      Ranked_Genes = 0L,
      Positive_Correlation_P05 = 0L,
      Negative_Correlation_P05 = 0L,
      Strong_Positive_AbsR_0_5 = 0L,
      Strong_Negative_AbsR_0_5 = 0L,
      Correlation_Method = CORRELATION_METHOD,
      stringsAsFactors = FALSE
    )

    summary_file <- write_csv_with_report_previews(
      summary_table,
      file.path(output_dir, "summary.csv")
    )

    return(list(
      analysis_name = scope_name,
      status = summary_table$Status,
      correlation_file = NA_character_,
      summary_file = summary_file,
      sample_count = sample_count
    ))
  }

  target_expression <- target_expression_all[keep_sample]
  expr_scope <- exprSet[, keep_sample, drop = FALSE]

  gene_sd <- apply(expr_scope, 1, sd, na.rm = TRUE)
  valid_gene <- is.finite(gene_sd) & gene_sd > 0
  target_sd <- sd(target_expression, na.rm = TRUE)

  rho <- rep(NA_real_, nrow(expr_scope))
  if (is.finite(target_sd) && target_sd > 0 && any(valid_gene)) {
    rho[valid_gene] <- as.numeric(cor(
      t(expr_scope[valid_gene, , drop = FALSE]),
      target_expression,
      method = CORRELATION_METHOD,
      use = "pairwise.complete.obs"
    ))
  }

  p_value <- compute_correlation_p_value(rho, sample_count)
  adj_p_value <- p.adjust(p_value, method = "BH")

  result_table <- data.frame(
    Dataset = DATASET_ID,
    Analysis_Name = scope_name,
    Target_Gene = TARGET_GENE,
    Feature_ID = gene_annotation$Feature_ID,
    GeneID = gene_annotation$GeneID,
    Symbol = gene_annotation$Symbol,
    Ensembl = gene_annotation$Ensembl,
    Entrez = gene_annotation$Entrez,
    Biotype = gene_annotation$Biotype,
    Expression_Mean_Log2_TPM_Plus_1 = rowMeans(expr_scope, na.rm = TRUE),
    Expression_SD_Log2_TPM_Plus_1 = gene_sd,
    Correlation = rho,
    Abs_Correlation = abs(rho),
    P.Value = p_value,
    adj.P.Val = adj_p_value,
    Samples_Used = sample_count,
    Correlation_Method = CORRELATION_METHOD,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  result_table <- result_table[
    order(-result_table$Abs_Correlation, result_table$P.Value, result_table$Symbol),
    ,
    drop = FALSE
  ]

  significant_table <- result_table[
    is.finite(result_table$P.Value) &
      result_table$P.Value < CORRELATION_P_CUTOFF,
    ,
    drop = FALSE
  ]

  summary_table <- data.frame(
    Dataset = DATASET_ID,
    Analysis_Name = scope_name,
    Target_Gene = TARGET_GENE,
    Status = "OK",
    Samples_Used = sample_count,
    Ranked_Genes = sum(is.finite(result_table$Correlation)),
    Positive_Correlation_P05 = sum(
      result_table$Correlation > 0 & result_table$P.Value < CORRELATION_P_CUTOFF,
      na.rm = TRUE
    ),
    Negative_Correlation_P05 = sum(
      result_table$Correlation < 0 & result_table$P.Value < CORRELATION_P_CUTOFF,
      na.rm = TRUE
    ),
    Strong_Positive_AbsR_0_5 = sum(result_table$Correlation >= 0.5, na.rm = TRUE),
    Strong_Negative_AbsR_0_5 = sum(result_table$Correlation <= -0.5, na.rm = TRUE),
    Target_Feature_ID = target_feature_id,
    Correlation_Method = CORRELATION_METHOD,
    Expression_Assay = TARGET_EXPRESSION_ASSAY,
    Expression_Log2_Transformed = TARGET_EXPRESSION_LOG2,
    stringsAsFactors = FALSE
  )

  correlation_file <- write_csv_with_report_previews(
    result_table,
    file.path(output_dir, "all_genes.csv")
  )
  significant_file <- write_csv_with_report_previews(
    significant_table,
    file.path(output_dir, "significant_correlated_genes.csv")
  )
  summary_file <- write_csv_with_report_previews(
    summary_table,
    file.path(output_dir, "summary.csv")
  )

  list(
    analysis_name = scope_name,
    status = "OK",
    correlation_file = correlation_file,
    significant_file = significant_file,
    summary_file = summary_file,
    sample_count = sample_count,
    summary = summary_table
  )
}

cat("\nRunning GSE282081 ATF3 correlation analyses...\n")
correlation_results <- lapply(CORRELATION_SCOPES_TO_RUN, compute_scope_correlation)
names(correlation_results) <- CORRELATION_SCOPES_TO_RUN

correlation_summary <- do.call(
  rbind,
  lapply(correlation_results, function(x) {
    if (!is.null(x$summary)) {
      return(x$summary)
    }

    read.csv(x$summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  })
)
rownames(correlation_summary) <- NULL

correlation_summary_file <- write_csv_with_report_previews(
  correlation_summary,
  file.path(TABLE_ROOT, "correlation_summary.csv")
)


# 5. ATF3相关性GSEA -----------------------------------------------------------

if (RUN_GSEA) {
  stopifnot(GENE_ID_TYPE %in% names(GENE_ID_COLUMN_BY_TYPE))
  stopifnot(GENE_ID_TYPE %in% names(MSIGDB_GENE_COLUMN_BY_TYPE))

  gsea_ready_results <- correlation_results[
    vapply(correlation_results, function(x) identical(x$status, "OK"), logical(1))
  ]

  if (length(gsea_ready_results) == 0) {
    stop("No correlation result is available for GSEA.")
  }

  MSIGDB_GENESET_CATALOG <- build_msigdb_geneset_catalog()
  RUNTIME_GENESETS_TO_RUN <- get_runtime_genesets_to_run()
  GSEA_GENESET_CONFIG <- select_msigdb_genesets(
    catalog = MSIGDB_GENESET_CATALOG,
    genesets_to_run = RUNTIME_GENESETS_TO_RUN
  )

  if (CLEAN_GSEA_OUTPUT_DIR) {
    unlink(
      file.path(TABLE_ROOT, names(gsea_ready_results), "GSEA"),
      recursive = TRUE,
      force = TRUE
    )
    unlink(
      file.path(TABLE_ROOT, "GSEA_summary"),
      recursive = TRUE,
      force = TRUE
    )
  }

  cat("\nLoading MSigDB gene sets for ATF3 correlation GSEA...\n")
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

  analysis_cache <- lapply(gsea_ready_results, function(x) {
    correlation_table <- read.csv(
      resolve_report_csv_file(x$correlation_file),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    gene_list <- prepare_gene_list(
      deg_table = correlation_table,
      gene_id_type = GENE_ID_TYPE,
      rank_metric_column = RANK_METRIC_COLUMN
    )

    list(
      correlation_file = x$correlation_file,
      gene_list = gene_list
    )
  })

  task_table <- expand.grid(
    Analysis_Name = names(analysis_cache),
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

  run_gsea_task <- function(task_id) {
    analysis_name <- task_table$Analysis_Name[task_id]
    geneset_name <- task_table$GeneSet_Name[task_id]
    analysis_input <- analysis_cache[[analysis_name]]
    cache <- geneset_cache[[geneset_name]]
    output_name <- cache$config$output_name
    output_dir_name <- sanitize_file_name(output_name)

    table_output_dir <- file.path(TABLE_ROOT, analysis_name, "GSEA", output_dir_name)
    dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

    gsea_run <- load_or_run_gsea(
      analysis_name = paste0("ATF3_", analysis_name),
      deg_file = analysis_input$correlation_file,
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

  cat("\nRunning ATF3 correlation GSEA tasks...\n")
  task_ids <- seq_len(nrow(task_table))
  gsea_summary_records <- run_parallel_tasks_with_progress(
    task_ids = task_ids,
    task_function = run_gsea_task,
    workers = parallel_strategy$task_workers
  )
  stop_on_parallel_errors(gsea_summary_records, task_ids = task_ids, label = "ATF3 correlation GSEA tasks")

  gsea_summary <- do.call(rbind, gsea_summary_records)
  rownames(gsea_summary) <- NULL

  gsea_summary_dir <- file.path(TABLE_ROOT, "GSEA_summary")
  dir.create(gsea_summary_dir, recursive = TRUE, showWarnings = FALSE)
  gsea_summary_file <- write_csv_with_report_previews(
    gsea_summary,
    file.path(gsea_summary_dir, "summary.csv")
  )

  cat("\nATF3 correlation GSEA summary:\n")
  print(
    gsea_summary[
      ,
      c(
        "Analysis_Name", "GeneSet_Name", "Source", "Ranked_Genes",
        "GSEA_Terms", "Positive_NES", "Negative_NES"
      )
    ],
    row.names = FALSE
  )
} else {
  gsea_summary_file <- NA_character_
}


# 6. 终端汇总 ------------------------------------------------------------------

cat("\nGSE282081 ATF3 correlation analysis finished.\n")
cat("Correlation summary: ", correlation_summary_file, "\n", sep = "")
if (RUN_GSEA) {
  cat("GSEA summary:        ", gsea_summary_file, "\n", sep = "")
}
cat("\nCorrelation summary:\n")
print(correlation_summary, row.names = FALSE)
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
