# GSE50760 ATF3 within-tissue functional analysis.
#
# First subsets one tissue group at a time, then evaluates ATF3 function inside
# that tissue context:
# 1. ATF3 median high/low differential expression and GSEA;
# 2. ATF3 continuous-expression Spearman correlation and GSEA.
# The liver metastasis scope directly addresses ATF3-associated functions in
# liver metastasis tumor tissue.


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
TARGET_GENE <- Sys.getenv("GSE50760_ATF3_TARGET_GENE", unset = "ATF3")

SE_RDS_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_se_raw.rds"
CLINICAL_EDIT_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_clinical_edit.csv"
SAMPLE_DESIGN_SCRIPT <- "scripts/GSE50760/00_prepare_sample_design.R"

LIMMA_FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
GSEA_FUNCTION_FILE <- "scripts/functions/gsea_common_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables", "ATF3_function")

TARGET_EXPRESSION_ASSAY <- Sys.getenv("GSE50760_ATF3_ASSAY", unset = "tpm")
TARGET_EXPRESSION_LOG2 <- parse_env_logical("GSE50760_ATF3_LOG2", TRUE)
COUNT_ASSAY_NAME <- Sys.getenv("GSE50760_ATF3_COUNT_ASSAY", unset = "counts")
GENE_BIOTYPE_FILTER <- "coding"
CORRELATION_METHOD <- "spearman"
MIN_SAMPLES_FOR_CORRELATION <- parse_env_integer("GSE50760_ATF3_MIN_N", 8L)
CORRELATION_P_CUTOFF <- as.numeric(Sys.getenv("GSE50760_ATF3_COR_P_CUTOFF", unset = "0.05"))

RUN_HIGH_LOW_DE <- parse_env_logical("GSE50760_ATF3_RUN_HIGH_LOW_DE", TRUE)
MIN_SAMPLES_PER_ATF3_GROUP <- parse_env_integer("GSE50760_ATF3_MIN_HIGH_LOW_GROUP_N", 5L)
HIGH_LOW_P_VALUE_COLUMN <- Sys.getenv("GSE50760_ATF3_HIGH_LOW_P_COLUMN", unset = "P.Value")
HIGH_LOW_P_VALUE_CUTOFF <- as.numeric(Sys.getenv("GSE50760_ATF3_HIGH_LOW_P_CUTOFF", unset = "0.05"))
HIGH_LOW_LOGFC_CUTOFF <- as.numeric(Sys.getenv("GSE50760_ATF3_HIGH_LOW_LOGFC_CUTOFF", unset = "0.5"))

CORRELATION_SCOPES_TO_RUN <- parse_env_vector(
  "GSE50760_ATF3_CORRELATION_SCOPES",
  c(
    "Liver_metastasis",
    "Primary_tumor",
    "Normal_colon"
  )
)

RUN_GSEA <- parse_env_logical("GSE50760_ATF3_RUN_GSEA", TRUE)

SPECIES <- Sys.getenv("GSE50760_ATF3_SPECIES", unset = "human")
GENE_ID_TYPE <- Sys.getenv("GSE50760_ATF3_GENE_ID_TYPE", unset = "ENTREZ")
RANK_METRIC_COLUMN <- "Correlation"

GSEA_SIGNIFICANCE_COLUMN <- Sys.getenv("GSE50760_ATF3_GSEA_P_COLUMN", unset = "pvalue")
GSEA_SIGNIFICANCE_CUTOFF <- as.numeric(Sys.getenv("GSE50760_ATF3_GSEA_P_CUTOFF", unset = "0.05"))
GSEA_PARAMS <- list(
  exponent = 1,
  minGSSize = parse_env_integer("GSE50760_ATF3_GSEA_MIN_GS_SIZE", 5L),
  maxGSSize = parse_env_integer("GSE50760_ATF3_GSEA_MAX_GS_SIZE", 500L),
  pvalueCutoff = GSEA_SIGNIFICANCE_CUTOFF,
  pAdjustMethod = "BH",
  verbose = parse_env_logical("GSE50760_ATF3_GSEA_VERBOSE", FALSE),
  nPerm = 1000,
  method = "multilevel",
  adaptive = FALSE,
  minPerm = 101,
  maxPerm = 1e5,
  pvalThreshold = 0.1
)

GSEA_GENESETS_TO_RUN <- parse_env_vector(
  "GSE50760_ATF3_GSEA_GENESETS",
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

READABLE_GENE_SYMBOLS <- parse_env_logical("GSE50760_ATF3_READABLE_SYMBOLS", TRUE)
USE_QS2_CACHE <- parse_env_logical("GSE50760_ATF3_USE_QS2_CACHE", TRUE)
REFRESH_QS2_CACHE <- parse_env_logical("GSE50760_ATF3_REFRESH_QS2_CACHE", TRUE)
CLEAN_GSEA_OUTPUT_DIR <- parse_env_logical("GSE50760_ATF3_CLEAN_GSEA", TRUE)

QS2_CACHE_DIR <- file.path("temporary", DATA_TYPE, DATASET_ID, "ATF3_GSEA_qs2_cache")
MSIGDB_REFERENCE_DIR <- file.path("data", "reference", "msigdb")
MSIGDB_REFERENCE_MAX_AGE_DAYS <- parse_env_integer("GSE50760_ATF3_MSIGDB_MAX_AGE_DAYS", 7L)

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 2. Packages and helpers ------------------------------------------------------

required_packages <- c("SummarizedExperiment", "limma", "edgeR", "clusterProfiler", "msigdbr", "qs2", "parallel")

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
  library(limma)
  library(edgeR)
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

compute_correlation_p_value <- function(rho, n) {
  p_value <- rep(NA_real_, length(rho))
  valid <- is.finite(rho) & n > 2
  rho_for_p <- pmin(pmax(rho[valid], -1 + 1e-15), 1 - 1e-15)
  t_stat <- rho_for_p * sqrt((n - 2) / pmax(1 - rho_for_p^2, .Machine$double.eps))
  p_value[valid] <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
  p_value
}

prepare_deg_output_table <- function(dat) {
  drop_columns <- c("Feature_ID", "Biotype", "Length")
  keep_columns <- setdiff(colnames(dat), drop_columns)
  dat <- dat[, keep_columns, drop = FALSE]

  if ("Entrez" %in% colnames(dat)) {
    if (is.numeric(dat$Entrez)) {
      dat$Entrez <- ifelse(
        is.na(dat$Entrez),
        "",
        format(dat$Entrez, scientific = FALSE, trim = TRUE)
      )
    } else {
      dat$Entrez <- as.character(dat$Entrez)
    }
  }

  dat
}

make_empty_high_low_summary <- function(scope_name, status, sample_count) {
  metadata <- get_scope_metadata(scope_name)
  data.frame(
    Dataset = DATASET_ID,
    Analysis_Name = scope_name,
    Scope_Group_Column = unname(metadata["Scope_Group_Column"]),
    Scope_Group_Value = unname(metadata["Scope_Group_Value"]),
    Target_Gene = TARGET_GENE,
    Contrast = paste0("ATF3_high_vs_ATF3_low_within_", scope_name),
    Status = status,
    Samples_Used = sample_count,
    ATF3_High_N = NA_integer_,
    ATF3_Low_N = NA_integer_,
    ATF3_Median_Cutoff = NA_real_,
    Up = 0L,
    Down = 0L,
    Total_Significant_Genes = 0L,
    P_Value_Column = HIGH_LOW_P_VALUE_COLUMN,
    P_Value_Cutoff = HIGH_LOW_P_VALUE_CUTOFF,
    LogFC_Cutoff = HIGH_LOW_LOGFC_CUTOFF,
    Genes_Filtered_By_EdgeR = NA_integer_,
    stringsAsFactors = FALSE
  )
}

run_scope_high_low_de <- function(scope_name) {
  keep_sample <- scope_definitions[[scope_name]](clinical_data)
  keep_sample[is.na(keep_sample)] <- FALSE
  sample_info <- clinical_data[keep_sample, , drop = FALSE]
  sample_count <- nrow(sample_info)
  metadata <- get_scope_metadata(scope_name)

  output_dir <- file.path(TABLE_ROOT, scope_name, "ATF3_high_low_DE", "DEG")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (sample_count < MIN_SAMPLES_PER_ATF3_GROUP * 2) {
    summary_table <- make_empty_high_low_summary(
      scope_name,
      "Skipped: too few samples for ATF3 high/low split",
      sample_count
    )
    summary_file <- write_csv_with_report_previews(summary_table, file.path(output_dir, "summary.csv"))
    return(list(
      analysis_name = scope_name,
      status = summary_table$Status,
      deg_file = NA_character_,
      summary_file = summary_file,
      summary = summary_table
    ))
  }

  atf3_expression <- target_expression_all[keep_sample]
  if (length(unique(atf3_expression[is.finite(atf3_expression)])) < 2) {
    summary_table <- make_empty_high_low_summary(
      scope_name,
      "Skipped: ATF3 expression has fewer than two unique values",
      sample_count
    )
    summary_file <- write_csv_with_report_previews(summary_table, file.path(output_dir, "summary.csv"))
    return(list(
      analysis_name = scope_name,
      status = summary_table$Status,
      deg_file = NA_character_,
      summary_file = summary_file,
      summary = summary_table
    ))
  }

  atf3_cutoff <- median(atf3_expression, na.rm = TRUE)
  atf3_group <- ifelse(atf3_expression >= atf3_cutoff, "ATF3_high", "ATF3_low")
  atf3_group <- factor(atf3_group, levels = c("ATF3_low", "ATF3_high"))
  group_counts <- table(atf3_group)

  if (any(group_counts < MIN_SAMPLES_PER_ATF3_GROUP)) {
    summary_table <- make_empty_high_low_summary(
      scope_name,
      "Skipped: too few samples in one ATF3 high/low group",
      sample_count
    )
    summary_table$ATF3_High_N <- unname(group_counts["ATF3_high"])
    summary_table$ATF3_Low_N <- unname(group_counts["ATF3_low"])
    summary_table$ATF3_Median_Cutoff <- atf3_cutoff
    summary_file <- write_csv_with_report_previews(summary_table, file.path(output_dir, "summary.csv"))
    return(list(
      analysis_name = scope_name,
      status = summary_table$Status,
      deg_file = NA_character_,
      summary_file = summary_file,
      summary = summary_table
    ))
  }

  sample_info$ATF3_Group <- atf3_group
  counts_scope <- counts[, sample_info$Sample_ID, drop = FALSE]
  design <- model.matrix(~ 0 + atf3_group)
  colnames(design) <- make.names(levels(atf3_group))
  rownames(design) <- sample_info$Sample_ID

  contrast_name <- "ATF3_high_vs_ATF3_low"
  contrast_matrix <- makeContrasts(
    contrasts = "ATF3_high - ATF3_low",
    levels = design
  )
  colnames(contrast_matrix) <- contrast_name

  analysis_data <- prepare_ngs_data(
    counts = counts_scope,
    gene_annotation = gene_annotation,
    group_list = atf3_group,
    design = design
  )

  fit <- lmFit(analysis_data$data, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  diff_results <- topTable(
    fit2,
    coef = contrast_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P",
    genelist = analysis_data$data$genes
  )

  stopifnot(HIGH_LOW_P_VALUE_COLUMN %in% colnames(diff_results))

  significant_index <- abs(diff_results$logFC) > HIGH_LOW_LOGFC_CUTOFF &
    diff_results[[HIGH_LOW_P_VALUE_COLUMN]] < HIGH_LOW_P_VALUE_CUTOFF
  up_index <- significant_index & diff_results$logFC > HIGH_LOW_LOGFC_CUTOFF
  down_index <- significant_index & diff_results$logFC < -HIGH_LOW_LOGFC_CUTOFF
  significant_results <- diff_results[significant_index, , drop = FALSE]

  all_results_file <- write_csv_with_report_previews(
    prepare_deg_output_table(diff_results),
    file.path(output_dir, "all_genes.csv")
  )
  significant_results_file <- write_csv_with_report_previews(
    prepare_deg_output_table(significant_results),
    file.path(output_dir, "significant_genes.csv")
  )

  summary_table <- data.frame(
    Dataset = DATASET_ID,
    Analysis_Name = scope_name,
    Scope_Group_Column = unname(metadata["Scope_Group_Column"]),
    Scope_Group_Value = unname(metadata["Scope_Group_Value"]),
    Target_Gene = TARGET_GENE,
    Contrast = paste0(contrast_name, "_within_", scope_name),
    Status = "OK",
    Samples_Used = sample_count,
    ATF3_High_N = unname(group_counts["ATF3_high"]),
    ATF3_Low_N = unname(group_counts["ATF3_low"]),
    ATF3_Median_Cutoff = atf3_cutoff,
    Up = sum(up_index),
    Down = sum(down_index),
    Total_Significant_Genes = sum(significant_index),
    P_Value_Column = HIGH_LOW_P_VALUE_COLUMN,
    P_Value_Cutoff = HIGH_LOW_P_VALUE_CUTOFF,
    LogFC_Cutoff = HIGH_LOW_LOGFC_CUTOFF,
    Genes_Filtered_By_EdgeR = analysis_data$filtered_genes,
    stringsAsFactors = FALSE
  )

  summary_file <- write_csv_with_report_previews(summary_table, file.path(output_dir, "summary.csv"))

  stopifnot(file.exists(all_results_file))
  stopifnot(file.exists(significant_results_file))
  stopifnot(file.exists(summary_file))

  list(
    analysis_name = scope_name,
    status = "OK",
    deg_file = all_results_file,
    significant_file = significant_results_file,
    summary_file = summary_file,
    summary = summary_table
  )
}

# 3. Read expression and sample information -----------------------------------

if (!file.exists(CLINICAL_EDIT_FILE)) {
  source(SAMPLE_DESIGN_SCRIPT)
}

dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(TABLE_ROOT, c("All_samples", "Tumor_samples")), recursive = TRUE, force = TRUE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_EDIT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot(TARGET_EXPRESSION_ASSAY %in% names(assays(se)))
stopifnot(COUNT_ASSAY_NAME %in% names(assays(se)))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot("tissue_class" %in% colnames(clinical_data))
stopifnot(all(colnames(se) %in% clinical_data$Sample_ID))

clinical_data <- clinical_data[
  match(colnames(se), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
stopifnot(all(clinical_data$Sample_ID == colnames(se)))

exprSet_all <- as.matrix(assay(se, TARGET_EXPRESSION_ASSAY))
counts_all <- as.matrix(assay(se, COUNT_ASSAY_NAME))
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
counts <- counts_all[rownames(gene_annotation), , drop = FALSE]


# 4. Correlation scopes --------------------------------------------------------

scope_definitions <- list(
  Liver_metastasis = function(sample_info) sample_info$tissue_class == "Liver_metastasis",
  Primary_tumor = function(sample_info) sample_info$tissue_class == "Primary_tumor",
  Normal_colon = function(sample_info) sample_info$tissue_class == "Normal_colon"
)

missing_scopes <- setdiff(CORRELATION_SCOPES_TO_RUN, names(scope_definitions))
if (length(missing_scopes) > 0) {
  stop("Undefined correlation scopes: ", paste(missing_scopes, collapse = ", "))
}

scope_metadata <- list(
  Liver_metastasis = c(Scope_Group_Column = "tissue_class", Scope_Group_Value = "Liver_metastasis"),
  Primary_tumor = c(Scope_Group_Column = "tissue_class", Scope_Group_Value = "Primary_tumor"),
  Normal_colon = c(Scope_Group_Column = "tissue_class", Scope_Group_Value = "Normal_colon")
)

clean_stale_scope_output_dirs <- function() {
  existing_dirs <- list.dirs(TABLE_ROOT, recursive = FALSE, full.names = FALSE)
  keep_dirs <- c(
    CORRELATION_SCOPES_TO_RUN,
    "GSEA_summary",
    "ATF3_high_low_DE_GSEA_summary"
  )
  stale_dirs <- setdiff(existing_dirs, keep_dirs)
  if (length(stale_dirs) > 0) {
    unlink(file.path(TABLE_ROOT, stale_dirs), recursive = TRUE, force = TRUE)
  }

  invisible(stale_dirs)
}

clean_stale_scope_output_dirs()

get_scope_metadata <- function(scope_name) {
  metadata <- scope_metadata[[scope_name]]
  if (is.null(metadata)) {
    return(c(Scope_Group_Column = "custom", Scope_Group_Value = scope_name))
  }

  metadata
}

make_scope_sample_manifest <- function() {
  display_columns <- intersect(
    c("Sample_ID", "Title", "Patient_ID", "tissue", "tissue_class", "ajcc_stage"),
    colnames(clinical_data)
  )

  do.call(
    rbind,
    lapply(CORRELATION_SCOPES_TO_RUN, function(scope_name) {
      keep_sample <- scope_definitions[[scope_name]](clinical_data)
      keep_sample[is.na(keep_sample)] <- FALSE
      metadata <- get_scope_metadata(scope_name)
      sample_table <- clinical_data[keep_sample, display_columns, drop = FALSE]

      data.frame(
        Dataset = DATASET_ID,
        Analysis_Name = scope_name,
        Scope_Group_Column = unname(metadata["Scope_Group_Column"]),
        Scope_Group_Value = unname(metadata["Scope_Group_Value"]),
        sample_table,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
  )
}

scope_sample_manifest <- make_scope_sample_manifest()
scope_sample_manifest_file <- write_csv_with_report_previews(
  scope_sample_manifest,
  file.path(TABLE_ROOT, "scope_sample_manifest.csv")
)

if (RUN_HIGH_LOW_DE) {
  cat("\nRunning GSE50760 within-tissue ATF3 high/low DE analyses...\n")
  high_low_results <- lapply(CORRELATION_SCOPES_TO_RUN, run_scope_high_low_de)
  names(high_low_results) <- CORRELATION_SCOPES_TO_RUN

  high_low_summary <- do.call(
    rbind,
    lapply(high_low_results, function(x) x$summary)
  )
  rownames(high_low_summary) <- NULL

  high_low_summary_file <- write_csv_with_report_previews(
    high_low_summary,
    file.path(TABLE_ROOT, "ATF3_high_low_DE_summary.csv")
  )
} else {
  high_low_results <- list()
  high_low_summary <- data.frame()
  high_low_summary_file <- NA_character_
}

compute_scope_correlation <- function(scope_name) {
  keep_sample <- scope_definitions[[scope_name]](clinical_data)
  keep_sample[is.na(keep_sample)] <- FALSE
  sample_ids <- clinical_data$Sample_ID[keep_sample]
  sample_count <- length(sample_ids)
  metadata <- get_scope_metadata(scope_name)

  output_dir <- file.path(TABLE_ROOT, scope_name, "correlation")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (sample_count < MIN_SAMPLES_FOR_CORRELATION) {
    summary_table <- data.frame(
      Dataset = DATASET_ID,
      Analysis_Name = scope_name,
      Scope_Group_Column = unname(metadata["Scope_Group_Column"]),
      Scope_Group_Value = unname(metadata["Scope_Group_Value"]),
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

    summary_file <- write_csv_with_report_previews(summary_table, file.path(output_dir, "summary.csv"))

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
    Scope_Group_Column = unname(metadata["Scope_Group_Column"]),
    Scope_Group_Value = unname(metadata["Scope_Group_Value"]),
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
    is.finite(result_table$P.Value) & result_table$P.Value < CORRELATION_P_CUTOFF,
    ,
    drop = FALSE
  ]

  summary_table <- data.frame(
    Dataset = DATASET_ID,
    Analysis_Name = scope_name,
    Scope_Group_Column = unname(metadata["Scope_Group_Column"]),
    Scope_Group_Value = unname(metadata["Scope_Group_Value"]),
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

cat("\nRunning GSE50760 ATF3 correlation analyses...\n")
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


# 5. Correlation GSEA ----------------------------------------------------------

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
    unlink(file.path(TABLE_ROOT, names(gsea_ready_results), "GSEA"), recursive = TRUE, force = TRUE)
    unlink(file.path(TABLE_ROOT, "GSEA_summary"), recursive = TRUE, force = TRUE)
  }

  cat("\nLoading MSigDB gene sets for GSE50760 ATF3 correlation GSEA...\n")
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
      Scope_Group_Column = unname(get_scope_metadata(analysis_name)["Scope_Group_Column"]),
      Scope_Group_Value = unname(get_scope_metadata(analysis_name)["Scope_Group_Value"]),
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

  cat("\nRunning GSE50760 ATF3 correlation GSEA tasks...\n")
  task_ids <- seq_len(nrow(task_table))
  gsea_summary_records <- run_parallel_tasks_with_progress(
    task_ids = task_ids,
    task_function = run_gsea_task,
    workers = parallel_strategy$task_workers
  )
  stop_on_parallel_errors(gsea_summary_records, task_ids = task_ids, label = "GSE50760 ATF3 correlation GSEA tasks")

  gsea_summary <- do.call(rbind, gsea_summary_records)
  rownames(gsea_summary) <- NULL

  gsea_summary_dir <- file.path(TABLE_ROOT, "GSEA_summary")
  dir.create(gsea_summary_dir, recursive = TRUE, showWarnings = FALSE)
  gsea_summary_file <- write_csv_with_report_previews(
    gsea_summary,
    file.path(gsea_summary_dir, "summary.csv")
  )

  cat("\nGSE50760 ATF3 correlation GSEA summary:\n")
  print(
    gsea_summary[
      ,
      c("Analysis_Name", "GeneSet_Name", "Source", "Ranked_Genes", "GSEA_Terms", "Positive_NES", "Negative_NES")
    ],
    row.names = FALSE
  )

  if (RUN_HIGH_LOW_DE && length(high_low_results) > 0) {
    high_low_ready_results <- high_low_results[
      vapply(high_low_results, function(x) identical(x$status, "OK"), logical(1))
    ]

    if (length(high_low_ready_results) > 0) {
      if (CLEAN_GSEA_OUTPUT_DIR) {
        unlink(
          file.path(TABLE_ROOT, names(high_low_ready_results), "ATF3_high_low_DE", "GSEA"),
          recursive = TRUE,
          force = TRUE
        )
        unlink(
          file.path(TABLE_ROOT, "ATF3_high_low_DE_GSEA_summary"),
          recursive = TRUE,
          force = TRUE
        )
      }

      old_rank_metric_column <- RANK_METRIC_COLUMN
      RANK_METRIC_COLUMN <- "t"

      high_low_analysis_cache <- lapply(high_low_ready_results, function(x) {
        deg_table <- read.csv(
          resolve_report_csv_file(x$deg_file),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        gene_list <- prepare_gene_list(
          deg_table = deg_table,
          gene_id_type = GENE_ID_TYPE,
          rank_metric_column = "t"
        )

        list(
          deg_file = x$deg_file,
          gene_list = gene_list
        )
      })

      high_low_task_table <- expand.grid(
        Analysis_Name = names(high_low_analysis_cache),
        GeneSet_Name = names(geneset_cache),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      high_low_parallel_strategy <- setup_parallel_strategy(
        total_tasks = nrow(high_low_task_table),
        max_workers = PARALLEL_WORKERS,
        inner_label = "ATF3 high/low GSEA nproc per task",
        nested_label = "Nested workers"
      )
      GSEA_INNER_NPROC <- high_low_parallel_strategy$inner_workers

      run_high_low_gsea_task <- function(task_id) {
        analysis_name <- high_low_task_table$Analysis_Name[task_id]
        geneset_name <- high_low_task_table$GeneSet_Name[task_id]
        analysis_input <- high_low_analysis_cache[[analysis_name]]
        cache <- geneset_cache[[geneset_name]]
        output_name <- cache$config$output_name
        output_dir_name <- sanitize_file_name(output_name)

        table_output_dir <- file.path(
          TABLE_ROOT,
          analysis_name,
          "ATF3_high_low_DE",
          "GSEA",
          output_dir_name
        )
        dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

        gsea_run <- load_or_run_gsea(
          analysis_name = paste0("ATF3_high_low_", analysis_name),
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
          Scope_Group_Column = unname(get_scope_metadata(analysis_name)["Scope_Group_Column"]),
          Scope_Group_Value = unname(get_scope_metadata(analysis_name)["Scope_Group_Value"]),
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

      cat("\nRunning GSE50760 ATF3 high/low DE GSEA tasks...\n")
      high_low_task_ids <- seq_len(nrow(high_low_task_table))
      high_low_gsea_summary_records <- run_parallel_tasks_with_progress(
        task_ids = high_low_task_ids,
        task_function = run_high_low_gsea_task,
        workers = high_low_parallel_strategy$task_workers
      )
      stop_on_parallel_errors(
        high_low_gsea_summary_records,
        task_ids = high_low_task_ids,
        label = "GSE50760 ATF3 high/low DE GSEA tasks"
      )

      RANK_METRIC_COLUMN <- old_rank_metric_column

      high_low_gsea_summary <- do.call(rbind, high_low_gsea_summary_records)
      rownames(high_low_gsea_summary) <- NULL

      high_low_gsea_summary_dir <- file.path(TABLE_ROOT, "ATF3_high_low_DE_GSEA_summary")
      dir.create(high_low_gsea_summary_dir, recursive = TRUE, showWarnings = FALSE)
      high_low_gsea_summary_file <- write_csv_with_report_previews(
        high_low_gsea_summary,
        file.path(high_low_gsea_summary_dir, "summary.csv")
      )

      cat("\nGSE50760 ATF3 high/low DE GSEA summary:\n")
      print(
        high_low_gsea_summary[
          ,
          c("Analysis_Name", "GeneSet_Name", "Source", "Ranked_Genes", "GSEA_Terms", "Positive_NES", "Negative_NES")
        ],
        row.names = FALSE
      )
    } else {
      high_low_gsea_summary_file <- NA_character_
    }
  } else {
    high_low_gsea_summary_file <- NA_character_
  }
} else {
  gsea_summary_file <- NA_character_
  high_low_gsea_summary_file <- NA_character_
}


# 6. Console summary -----------------------------------------------------------

cat("\nGSE50760 ATF3 correlation analysis finished.\n")
cat("Scope sample manifest: ", scope_sample_manifest_file, "\n", sep = "")
if (RUN_HIGH_LOW_DE) {
  cat("ATF3 high/low DE summary: ", high_low_summary_file, "\n", sep = "")
}
cat("Correlation summary: ", correlation_summary_file, "\n", sep = "")
if (RUN_GSEA) {
  cat("Correlation GSEA summary: ", gsea_summary_file, "\n", sep = "")
  cat("High/low DE GSEA summary: ", high_low_gsea_summary_file, "\n", sep = "")
}
if (RUN_HIGH_LOW_DE) {
  cat("\nATF3 high/low DE summary:\n")
  print(high_low_summary, row.names = FALSE)
}
cat("\nCorrelation summary:\n")
print(correlation_summary, row.names = FALSE)
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
