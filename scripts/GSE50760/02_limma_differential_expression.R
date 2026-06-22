# GSE50760 paired differential expression analysis.
#
# Runs edgeR + voom + limma for each analysis_ column created by script 00.
# When complete patient pairs are present, the model includes Patient_ID as a
# fixed effect, which is appropriate for the matched trio design.


# 0. Config -------------------------------------------------------------------

DATASET_ID <- "GSE50760"
DATA_TYPE <- "ngs"
TARGET_GENE <- "ATF3"

SE_RDS_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_clinical_edit.csv"
SAMPLE_DESIGN_SCRIPT <- "scripts/GSE50760/00_prepare_sample_design.R"

FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"

GENE_BIOTYPE_FILTER <- "coding"
P_VALUE_COLUMN <- "P.Value"
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

OUTPUT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID, "tables")
CLEAN_DEG_OUTPUT_DIR <- TRUE
USE_PAIRED_DESIGN_WHEN_AVAILABLE <- TRUE

OUTPUT_DROP_COLUMNS <- c("Feature_ID", "Biotype", "Length")

options(width = 200)

SUMMARY_DISPLAY_COLUMNS <- c(
  "Analysis_Name", "Contrast", "Samples_Used", "Patients_Used", "Paired_Design_Used",
  "Up", "Down", "Total_Significant_Genes",
  "Target_Gene_LogFC", "Target_Gene_P_Value", "Target_Gene_Status"
)


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
})

source(FUNCTION_FILE)
source(REPORT_TABLE_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

prepare_deg_output_table <- function(dat) {
  keep_columns <- setdiff(colnames(dat), OUTPUT_DROP_COLUMNS)
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

get_first_value <- function(dat, column_name, default = NA_character_) {
  if (!column_name %in% colnames(dat)) {
    return(default)
  }

  value <- dat[[column_name]][1]
  if (is.na(value)) {
    return(default)
  }

  as.character(value)
}

get_first_numeric <- function(dat, column_name) {
  if (!column_name %in% colnames(dat)) {
    return(NA_real_)
  }

  as.numeric(dat[[column_name]][1])
}

extract_target_gene_result <- function(
    diff_results,
    analysis_name,
    contrast_name,
    samples_used,
    patients_used,
    experiment_group,
    control_group,
    paired_design_used) {
  target_index <- integer(0)

  if ("Symbol" %in% colnames(diff_results)) {
    target_index <- which(toupper(trimws(as.character(diff_results$Symbol))) == toupper(TARGET_GENE))
  }

  if (length(target_index) == 0 && "GeneID" %in% colnames(diff_results)) {
    target_index <- which(toupper(trimws(as.character(diff_results$GeneID))) == toupper(TARGET_GENE))
  }

  if (length(target_index) == 0) {
    return(data.frame(
      Dataset = DATASET_ID,
      Analysis_Name = analysis_name,
      Contrast = contrast_name,
      Target_Gene = TARGET_GENE,
      Status = "Not found after filtering",
      Samples_Used = samples_used,
      Patients_Used = patients_used,
      Paired_Design_Used = paired_design_used,
      Experiment_Group = experiment_group,
      Control_Group = control_group,
      Feature_ID = NA_character_,
      GeneID = NA_character_,
      Symbol = NA_character_,
      Ensembl = NA_character_,
      Entrez = NA_character_,
      logFC = NA_real_,
      AveExpr = NA_real_,
      t = NA_real_,
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      B = NA_real_,
      Rank_By_PValue = NA_integer_,
      Significant_By_Config = FALSE,
      Direction = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  target_index <- target_index[1]
  target_row <- diff_results[target_index, , drop = FALSE]
  logfc <- get_first_numeric(target_row, "logFC")
  p_value <- get_first_numeric(target_row, P_VALUE_COLUMN)
  significant <- is.finite(logfc) &&
    is.finite(p_value) &&
    abs(logfc) > LOGFC_CUTOFF &&
    p_value < P_VALUE_CUTOFF

  data.frame(
    Dataset = DATASET_ID,
    Analysis_Name = analysis_name,
    Contrast = contrast_name,
    Target_Gene = TARGET_GENE,
    Status = "OK",
    Samples_Used = samples_used,
    Patients_Used = patients_used,
    Paired_Design_Used = paired_design_used,
    Experiment_Group = experiment_group,
    Control_Group = control_group,
    Feature_ID = get_first_value(target_row, "Feature_ID"),
    GeneID = get_first_value(target_row, "GeneID"),
    Symbol = get_first_value(target_row, "Symbol"),
    Ensembl = get_first_value(target_row, "Ensembl"),
    Entrez = get_first_value(target_row, "Entrez"),
    logFC = logfc,
    AveExpr = get_first_numeric(target_row, "AveExpr"),
    t = get_first_numeric(target_row, "t"),
    P.Value = get_first_numeric(target_row, "P.Value"),
    adj.P.Val = get_first_numeric(target_row, "adj.P.Val"),
    B = get_first_numeric(target_row, "B"),
    Rank_By_PValue = target_index,
    Significant_By_Config = significant,
    Direction = ifelse(
      is.na(logfc),
      NA_character_,
      ifelse(logfc > 0, "Higher_in_experiment_group", "Lower_in_experiment_group")
    ),
    stringsAsFactors = FALSE
  )
}

can_use_paired_design <- function(sample_info, group_list) {
  if (!USE_PAIRED_DESIGN_WHEN_AVAILABLE || !"Patient_ID" %in% colnames(sample_info)) {
    return(FALSE)
  }

  pair_table <- table(sample_info$Patient_ID, group_list)
  all(rowSums(pair_table) == 2) && all(pair_table == 1)
}

fit_limma_model <- function(exprSet, gene_annotation, sample_info, group_list, experiment_group, control_group) {
  paired_design_used <- can_use_paired_design(sample_info, group_list)

  if (paired_design_used) {
    sample_info$Patient_ID <- factor(sample_info$Patient_ID)
    group_list <- factor(as.character(group_list), levels = c(control_group, experiment_group))

    design <- model.matrix(~ sample_info$Patient_ID + group_list)
    coefficient_name <- paste0("group_list", experiment_group)

    if (!coefficient_name %in% colnames(design)) {
      stop("Cannot find paired model coefficient: ", coefficient_name)
    }

    analysis_data <- prepare_ngs_data(
      counts = exprSet,
      gene_annotation = gene_annotation,
      group_list = group_list,
      design = design
    )

    fit <- lmFit(analysis_data$data, design)
    fit <- eBayes(fit)
    contrast_name <- paste0(experiment_group, "_vs_", control_group)

    diff_results <- topTable(
      fit,
      coef = coefficient_name,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P",
      genelist = analysis_data$data$genes
    )

    return(list(
      diff_results = diff_results,
      analysis_data = analysis_data,
      contrast_name = contrast_name,
      paired_design_used = TRUE
    ))
  }

  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- make.names(levels(group_list))
  rownames(design) <- colnames(exprSet)

  contrast_name <- paste0(experiment_group, "_vs_", control_group)
  contrast_formula <- paste0(make.names(experiment_group), " - ", make.names(control_group))
  contrast.matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  colnames(contrast.matrix) <- contrast_name

  analysis_data <- prepare_ngs_data(
    counts = exprSet,
    gene_annotation = gene_annotation,
    group_list = group_list,
    design = design
  )

  fit <- lmFit(analysis_data$data, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  diff_results <- topTable(
    fit2,
    coef = contrast_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P",
    genelist = analysis_data$data$genes
  )

  list(
    diff_results = diff_results,
    analysis_data = analysis_data,
    contrast_name = contrast_name,
    paired_design_used = FALSE
  )
}


# 2. Read and check inputs -----------------------------------------------------

if (!file.exists(CLINICAL_FILE)) {
  source(SAMPLE_DESIGN_SCRIPT)
}

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot("Patient_ID" %in% colnames(clinical_data))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

analysis_designs <- get_analysis_designs(clinical_data)
exprSet_all <- get_assay_matrix(se, DATA_TYPE)
stopifnot(is.numeric(exprSet_all))

missing_samples <- setdiff(colnames(exprSet_all), clinical_data$Sample_ID)
stopifnot(length(missing_samples) == 0)

sample_info_all <- clinical_data[
  match(colnames(exprSet_all), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
rownames(sample_info_all) <- sample_info_all$Sample_ID
stopifnot(all(sample_info_all$Sample_ID == colnames(exprSet_all)))

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

gene_filter <- filter_genes_by_biotype(
  exprSet = exprSet_all,
  gene_annotation = gene_annotation,
  biotype_filter = GENE_BIOTYPE_FILTER
)

exprSet_all <- gene_filter$exprSet
gene_annotation <- gene_filter$gene_annotation


# 3. Run limma for each design -------------------------------------------------

run_one_limma_analysis <- function(i) {
  analysis_name <- analysis_designs$Analysis_Name[i]
  analysis_column <- analysis_designs$Column_Name[i]
  analysis_column_index <- analysis_designs$Column_Index[i]
  experiment_group <- analysis_designs$Experiment_Group[i]

  design_samples <- prepare_design_samples(
    sample_info = sample_info_all,
    group_column_index = analysis_column_index,
    experiment_group = experiment_group
  )

  sample_info <- design_samples$sample_info
  group_list <- design_samples$group_list
  control_group <- design_samples$control_group
  exprSet <- exprSet_all[, sample_info$Sample_ID, drop = FALSE]

  fit_result <- fit_limma_model(
    exprSet = exprSet,
    gene_annotation = gene_annotation,
    sample_info = sample_info,
    group_list = group_list,
    experiment_group = experiment_group,
    control_group = control_group
  )

  diff_results <- fit_result$diff_results
  analysis_data <- fit_result$analysis_data
  contrast_name <- fit_result$contrast_name
  paired_design_used <- fit_result$paired_design_used

  stopifnot(P_VALUE_COLUMN %in% colnames(diff_results))

  significant_index <- abs(diff_results$logFC) > LOGFC_CUTOFF &
    diff_results[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF
  up_index <- significant_index & diff_results$logFC > LOGFC_CUTOFF
  down_index <- significant_index & diff_results$logFC < -LOGFC_CUTOFF
  significant_results <- diff_results[significant_index, , drop = FALSE]

  target_result <- extract_target_gene_result(
    diff_results = diff_results,
    analysis_name = analysis_name,
    contrast_name = contrast_name,
    samples_used = nrow(sample_info),
    patients_used = length(unique(sample_info$Patient_ID)),
    experiment_group = experiment_group,
    control_group = control_group,
    paired_design_used = paired_design_used
  )

  diff_results_output <- prepare_deg_output_table(diff_results)
  significant_results_output <- prepare_deg_output_table(significant_results)

  analysis_output_dir <- file.path(OUTPUT_ROOT, analysis_name, "DEG")
  dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

  if (CLEAN_DEG_OUTPUT_DIR) {
    unlink(list.files(analysis_output_dir, pattern = "[.](csv|tex|md)$", full.names = TRUE))
    unlink(file.path(analysis_output_dir, c("csv", "md", "tex")), recursive = TRUE, force = TRUE)
  }

  all_results_file <- file.path(analysis_output_dir, "all_genes.csv")
  significant_results_file <- file.path(analysis_output_dir, "significant_genes.csv")
  summary_file <- file.path(analysis_output_dir, "summary.csv")

  summary_table <- data.frame(
    Dataset = DATASET_ID,
    Data_Type = DATA_TYPE,
    Analysis_Name = analysis_name,
    Analysis_Column = analysis_column,
    Contrast = contrast_name,
    Samples_Used = nrow(sample_info),
    Patients_Used = length(unique(sample_info$Patient_ID)),
    Experiment_Group = experiment_group,
    Control_Group = control_group,
    Paired_Design_Used = paired_design_used,
    Up = sum(up_index),
    Down = sum(down_index),
    Total_Significant_Genes = sum(significant_index),
    P_Value_Column = P_VALUE_COLUMN,
    P_Value_Cutoff = P_VALUE_CUTOFF,
    LogFC_Cutoff = LOGFC_CUTOFF,
    Genes_Filtered_By_EdgeR = analysis_data$filtered_genes,
    Gene_Biotype_Filter = GENE_BIOTYPE_FILTER,
    Target_Gene = TARGET_GENE,
    Target_Gene_LogFC = target_result$logFC,
    Target_Gene_P_Value = target_result[[P_VALUE_COLUMN]],
    Target_Gene_adj_P_Val = target_result$adj.P.Val,
    Target_Gene_Significant = target_result$Significant_By_Config,
    Target_Gene_Status = target_result$Status,
    stringsAsFactors = FALSE
  )

  all_results_file <- write_csv_with_report_previews(diff_results_output, all_results_file)
  significant_results_file <- write_csv_with_report_previews(significant_results_output, significant_results_file)
  summary_file <- write_csv_with_report_previews(summary_table, summary_file)

  saved_significant_results <- read.csv(significant_results_file, stringsAsFactors = FALSE, check.names = FALSE)
  saved_summary <- read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)

  stopifnot(nrow(saved_summary) == 1)
  stopifnot(nrow(saved_significant_results) == summary_table$Total_Significant_Genes)
  stopifnot(saved_summary$Total_Significant_Genes == saved_summary$Up + saved_summary$Down)

  list(
    summary = summary_table,
    target_result = target_result,
    output_check = data.frame(
      Analysis_Name = analysis_name,
      All_Genes_File_Exists = file.exists(all_results_file),
      Significant_Genes_File_Exists = file.exists(significant_results_file),
      Summary_File_Exists = file.exists(summary_file),
      Significant_Genes_Rows = nrow(saved_significant_results),
      Summary_Total_Significant_Genes = saved_summary$Total_Significant_Genes,
      stringsAsFactors = FALSE
    )
  )
}

cat("\nRunning GSE50760 paired limma differential expression analyses...\n")
parallel_strategy <- setup_parallel_strategy(
  total_tasks = nrow(analysis_designs),
  inner_label = "limma inner workers",
  nested_label = "Nested workers"
)

analysis_task_results <- run_indexed_tasks_with_progress(
  total_tasks = nrow(analysis_designs),
  workers = parallel_strategy$task_workers,
  task_function = run_one_limma_analysis
)
stop_on_parallel_errors(
  analysis_task_results,
  task_ids = analysis_designs$Analysis_Name,
  label = "GSE50760 limma analyses"
)

summary_table <- do.call(rbind, lapply(analysis_task_results, `[[`, "summary"))
target_table <- do.call(rbind, lapply(analysis_task_results, `[[`, "target_result"))
output_check <- do.call(rbind, lapply(analysis_task_results, `[[`, "output_check"))

stopifnot(all(output_check$All_Genes_File_Exists))
stopifnot(all(output_check$Significant_Genes_File_Exists))
stopifnot(all(output_check$Summary_File_Exists))
stopifnot(
  all(output_check$Significant_Genes_Rows ==
      output_check$Summary_Total_Significant_Genes)
)

summary_output_dir <- file.path(OUTPUT_ROOT, "DEG_summary")
target_output_dir <- file.path(OUTPUT_ROOT, "ATF3_DEG_summary")
dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(target_output_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- write_csv_with_report_previews(
  summary_table,
  file.path(summary_output_dir, "summary.csv")
)
target_file <- write_csv_with_report_previews(
  target_table,
  file.path(target_output_dir, "ATF3_by_contrast.csv")
)

cat("\nGSE50760 limma analyses finished.\n")
cat("DEG summary:   ", summary_file, "\n", sep = "")
cat("ATF3 summary:  ", target_file, "\n\n", sep = "")
print(summary_table[, SUMMARY_DISPLAY_COLUMNS], row.names = FALSE)
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
