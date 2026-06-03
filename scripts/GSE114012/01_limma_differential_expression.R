# GSE114012 differential expression analysis with limma
#
# The clinical configuration file can contain multiple analysis_XXX columns.
# Each analysis_XXX column defines one differential expression design, where
# XXX is the experiment group. Blank cells indicate samples excluded from that
# specific design.


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"  # 可选："microarray" 或 "ngs"

SE_RDS_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_clinical_edit.csv"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"

# 显著差异筛选阈值
# P_VALUE_COLUMN可选："P.Value" 或 "adj.P.Val"
P_VALUE_COLUMN <- "P.Value"
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.2

OUTPUT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID, "tables")


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
})

if (DATA_TYPE == "ngs") {
  suppressPackageStartupMessages(library(edgeR))
}

source(FUNCTION_FILE)


# 2. 读取数据 -----------------------------------------------------------------

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

class(se)
dim(se)
names(assays(se))
head(clinical_data)


# 3. 检查样本和分析配置 --------------------------------------------------------

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot(DATA_TYPE %in% c("microarray", "ngs"))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

analysis_designs <- get_analysis_designs(clinical_data)

analysis_designs

exprSet_all <- get_assay_matrix(se, DATA_TYPE)

dim(exprSet_all)
head(exprSet_all[, 1:min(4, ncol(exprSet_all))])
summary(as.vector(exprSet_all))

stopifnot(is.numeric(exprSet_all))

missing_samples <- setdiff(colnames(exprSet_all), clinical_data$Sample_ID)
stopifnot(length(missing_samples) == 0)

sample_info_all <- clinical_data[
  match(colnames(exprSet_all), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
rownames(sample_info_all) <- sample_info_all$Sample_ID

head(sample_info_all)
stopifnot(all(sample_info_all$Sample_ID == colnames(exprSet_all)))


# 4. 准备基因注释 --------------------------------------------------------------

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

head(gene_annotation)


# 5. 逐个分析设计运行limma -----------------------------------------------------

for (i in seq_len(nrow(analysis_designs))) {
  analysis_base_name <- analysis_designs$Analysis_Base_Name[i]
  analysis_name <- analysis_designs$Analysis_Name[i]
  analysis_column <- analysis_designs$Column_Name[i]
  analysis_column_index <- analysis_designs$Column_Index[i]
  experiment_group <- analysis_designs$Experiment_Group[i]

  cat("\n\n")
  cat("Analysis:", analysis_name, "\n")
  cat("Design column:", analysis_column, "\n")
  cat("Experiment group:", experiment_group, "\n")

  design_samples <- prepare_design_samples(
    sample_info = sample_info_all,
    group_column_index = analysis_column_index,
    experiment_group = experiment_group
  )

  sample_info <- design_samples$sample_info
  group_list <- design_samples$group_list
  control_group <- design_samples$control_group

  print(table(group_list))
  print(head(sample_info))

  exprSet <- exprSet_all[, sample_info$Sample_ID, drop = FALSE]

  print(dim(exprSet))
  print(head(exprSet[, 1:min(4, ncol(exprSet))]))

  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- make.names(levels(group_list))
  rownames(design) <- colnames(exprSet)

  print(design)

  contrast_name <- paste0(experiment_group, "_vs_", control_group)
  contrast_formula <- paste0(
    make.names(experiment_group),
    " - ",
    make.names(control_group)
  )

  contrast.matrix <- makeContrasts(
    contrasts = contrast_formula,
    levels = design
  )
  colnames(contrast.matrix) <- contrast_name

  print(contrast.matrix)

  if (DATA_TYPE == "ngs") {
    analysis_data <- prepare_ngs_data(
      counts = exprSet,
      gene_annotation = gene_annotation,
      group_list = group_list,
      design = design
    )
    limma_input <- analysis_data$data
    genes_for_output <- limma_input$genes
  }

  if (DATA_TYPE == "microarray") {
    print(get_distribution_diagnostics(limma_input))
    print(analysis_data$log2_transformed)
    print(analysis_data$normalized_between_arrays)
  }

  fit <- lmFit(limma_input, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  diff_results <- topTable(
    fit2,
    coef = contrast_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P",
    genelist = genes_for_output
  )

  print(dim(diff_results))
  print(head(diff_results))
  print(colnames(diff_results))

  stopifnot(P_VALUE_COLUMN %in% colnames(diff_results))

  up_index <- diff_results$logFC > LOGFC_CUTOFF &
    diff_results[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF

  down_index <- diff_results$logFC < -LOGFC_CUTOFF &
    diff_results[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF

  significant_results <- diff_results[up_index | down_index, , drop = FALSE]

  de_summary <- data.frame(
    Up = sum(up_index),
    Down = sum(down_index),
    Not_Significant = nrow(diff_results) - sum(up_index) - sum(down_index)
  )

  print(de_summary)
  print(dim(significant_results))
  print(head(significant_results))

  analysis_output_dir <- file.path(OUTPUT_ROOT, analysis_name)
  dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

  safe_contrast_name <- paste0(
    sanitize_file_name(experiment_group),
    "_vs_",
    sanitize_file_name(control_group)
  )

  all_results_file <- file.path(
    analysis_output_dir,
    paste0(
      DATASET_ID, "_", analysis_name, "_limma_",
      safe_contrast_name, "_all_genes.csv"
    )
  )

  significant_results_file <- file.path(
    analysis_output_dir,
    paste0(
      DATASET_ID, "_", analysis_name, "_limma_",
      safe_contrast_name, "_significant_genes.csv"
    )
  )

  summary_file <- file.path(
    analysis_output_dir,
    paste0(
      DATASET_ID, "_", analysis_name, "_limma_",
      safe_contrast_name, "_summary.csv"
    )
  )

  summary_table <- data.frame(
    Dataset = DATASET_ID,
    Data_Type = DATA_TYPE,
    Analysis_Base_Name = analysis_base_name,
    Analysis_Name = analysis_name,
    Analysis_Column = analysis_column,
    Analysis_Column_Index = analysis_column_index,
    Contrast = contrast_name,
    Control_Group = control_group,
    Experiment_Group = experiment_group,
    Samples_Used = nrow(sample_info),
    Total_Genes_Input = nrow(exprSet_all),
    Total_Genes_Analyzed = nrow(diff_results),
    Up = de_summary$Up,
    Down = de_summary$Down,
    Not_Significant = de_summary$Not_Significant,
    P_Value_Column = P_VALUE_COLUMN,
    P_Value_Cutoff = P_VALUE_CUTOFF,
    LogFC_Cutoff = LOGFC_CUTOFF,
    Microarray_Log2_Transformed = analysis_data$log2_transformed,
    Microarray_Normalized_Between_Arrays = analysis_data$normalized_between_arrays,
    Microarray_Final_Median_Spread = analysis_data$median_spread,
    Microarray_Final_IQR_Spread = analysis_data$iqr_spread,
    NGS_Filtered_Genes = analysis_data$filtered_genes,
    stringsAsFactors = FALSE
  )

  print(summary_table)

  write.csv(diff_results, all_results_file, row.names = FALSE)
  write.csv(significant_results, significant_results_file, row.names = FALSE)
  write.csv(summary_table, summary_file, row.names = FALSE)
}
