# GSE114012差异表达分析
#
# 临床配置文件中可以包含多个analysis_XXX列。
# 每个analysis_XXX列定义一套差异分析设计；空白样本不参与该设计。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"  # 可选："microarray" 或 "ngs"

SE_RDS_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_clinical_edit.csv"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"

# 基因类型过滤
# 可选："coding", "non_coding", "all"
# coding口径与数据准备阶段一致：protein_coding/protein-coding + IG/TR gene。
GENE_BIOTYPE_FILTER <- "coding"

# 显著差异筛选阈值
# P_VALUE_COLUMN可选："P.Value" 或 "adj.P.Val"
P_VALUE_COLUMN <- "P.Value"
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

OUTPUT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID, "tables")

# 重跑时清理当前<analysis_name>/DEG目录内旧表格，避免旧长文件名残留。
CLEAN_DEG_OUTPUT_DIR <- TRUE

# 输出结果中不再保留的辅助注释列。
OUTPUT_DROP_COLUMNS <- c("Feature_ID", "Biotype", "Length")

options(width = 200)

SUMMARY_DISPLAY_COLUMNS <- c(
  "Analysis_Name", "Contrast", "Samples_Used",
  "Up", "Down", "Total_Significant_Genes",
  "P_Value_Column", "P_Value_Cutoff", "LogFC_Cutoff"
)


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
})

if (DATA_TYPE == "ngs") {
  suppressPackageStartupMessages(library(edgeR))
}

source(FUNCTION_FILE)
source(REPORT_TABLE_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()


# 1.1 输出表格整理函数 --------------------------------------------------------

prepare_deg_output_table <- function(dat) {
  # Feature_ID/Biotype/Length主要用于内部注释与过滤，最终DEG表不再保存。
  keep_columns <- setdiff(colnames(dat), OUTPUT_DROP_COLUMNS)
  dat <- dat[, keep_columns, drop = FALSE]

  # Entrez编号作为ID展示，统一按字符写出，避免CSV中出现科学计数法。
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


# 2. 读取和检查数据 ------------------------------------------------------------

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot(DATA_TYPE %in% c("microarray", "ngs"))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

analysis_designs <- get_analysis_designs(clinical_data)
stopifnot(all(analysis_designs$Analysis_Order == seq_len(nrow(analysis_designs))))
stopifnot(all(diff(analysis_designs$Column_Index) > 0))

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


# 3. 准备基因注释 --------------------------------------------------------------

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


# 4. 根据基因类型筛选分析基因 --------------------------------------------------

gene_filter <- filter_genes_by_biotype(
  exprSet = exprSet_all,
  gene_annotation = gene_annotation,
  biotype_filter = GENE_BIOTYPE_FILTER
)

exprSet_all <- gene_filter$exprSet
gene_annotation <- gene_filter$gene_annotation


# 5. 逐个分析设计运行limma -----------------------------------------------------

run_one_limma_analysis <- function(i) {
  analysis_order <- analysis_designs$Analysis_Order[i]
  analysis_base_name <- analysis_designs$Analysis_Base_Name[i]
  analysis_duplicate_order <- analysis_designs$Duplicate_Order[i]
  analysis_name <- analysis_designs$Analysis_Name[i]
  analysis_column <- analysis_designs$Column_Name[i]
  analysis_column_index <- analysis_designs$Column_Index[i]
  experiment_group <- analysis_designs$Experiment_Group[i]

  # 5.1 提取当前设计的样本，并剔除该设计列中的空白样本
  design_samples <- prepare_design_samples(
    sample_info = sample_info_all,
    group_column_index = analysis_column_index,
    experiment_group = experiment_group
  )

  sample_info <- design_samples$sample_info
  group_list <- design_samples$group_list
  control_group <- design_samples$control_group

  exprSet <- exprSet_all[, sample_info$Sample_ID, drop = FALSE]

  # 5.2 建立limma设计矩阵与比较矩阵
  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- make.names(levels(group_list))
  rownames(design) <- colnames(exprSet)

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

  # 5.3 根据数据类型准备limma输入
  if (DATA_TYPE == "microarray") {
    analysis_data <- prepare_microarray_data(exprSet)
    limma_input <- analysis_data$data
    genes_for_output <- gene_annotation
  }

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

  # 5.4 拟合线性模型并提取limma标准结果
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

  stopifnot(P_VALUE_COLUMN %in% colnames(diff_results))

  # 5.5 按配置阈值筛选显著差异基因
  significant_index <- abs(diff_results$logFC) > LOGFC_CUTOFF &
    diff_results[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF

  up_index <- significant_index &
    diff_results$logFC > LOGFC_CUTOFF

  down_index <- significant_index &
    diff_results$logFC < -LOGFC_CUTOFF

  significant_results <- diff_results[significant_index, , drop = FALSE]

  diff_results_output <- prepare_deg_output_table(diff_results)
  significant_results_output <- prepare_deg_output_table(significant_results)

  de_summary <- data.frame(
    Up = sum(up_index),
    Down = sum(down_index),
    Total_Significant_Genes = sum(significant_index)
  )

  stopifnot(nrow(significant_results) == de_summary$Total_Significant_Genes)
  stopifnot(de_summary$Total_Significant_Genes == de_summary$Up + de_summary$Down)

  # 5.6 保存当前设计的结果文件
  # 数据集编号、分析名和DEG类型已经由目录表达，文件名只保留结果类型。
  analysis_output_dir <- file.path(OUTPUT_ROOT, analysis_name, "DEG")
  dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

  if (CLEAN_DEG_OUTPUT_DIR) {
    unlink(list.files(
      analysis_output_dir,
      pattern = "[.](csv|tex|md)$",
      full.names = TRUE
    ))
  }

  all_results_file <- file.path(analysis_output_dir, "all_genes.csv")

  significant_results_file <- file.path(analysis_output_dir, "significant_genes.csv")

  summary_file <- file.path(analysis_output_dir, "summary.csv")

  summary_table <- data.frame(
    Dataset = DATASET_ID,
    Data_Type = DATA_TYPE,
    Analysis_Name = analysis_name,
    Contrast = contrast_name,
    Samples_Used = nrow(sample_info),
    Up = de_summary$Up,
    Down = de_summary$Down,
    Total_Significant_Genes = de_summary$Total_Significant_Genes,
    P_Value_Column = P_VALUE_COLUMN,
    P_Value_Cutoff = P_VALUE_CUTOFF,
    LogFC_Cutoff = LOGFC_CUTOFF,
    Microarray_Log2_Transformed = analysis_data$log2_transformed,
    stringsAsFactors = FALSE
  )

  write_csv_with_report_previews(diff_results_output, all_results_file)
  write_csv_with_report_previews(significant_results_output, significant_results_file)
  write_csv_with_report_previews(summary_table, summary_file)

  stopifnot(file.exists(all_results_file))
  stopifnot(file.exists(significant_results_file))
  stopifnot(file.exists(summary_file))

  saved_significant_results <- read.csv(
    significant_results_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  saved_summary <- read.csv(
    summary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stopifnot(nrow(saved_summary) == 1)
  stopifnot(nrow(saved_significant_results) == summary_table$Total_Significant_Genes)
  stopifnot(!any(OUTPUT_DROP_COLUMNS %in% colnames(saved_significant_results)))
  stopifnot(saved_summary$Up == summary_table$Up)
  stopifnot(saved_summary$Down == summary_table$Down)
  stopifnot(saved_summary$Total_Significant_Genes == summary_table$Total_Significant_Genes)
  stopifnot(saved_summary$Total_Significant_Genes == saved_summary$Up + saved_summary$Down)

  list(
    summary = summary_table,
    output_check = data.frame(
      Analysis_Name = analysis_name,
      All_Genes_File_Exists = file.exists(all_results_file),
      Significant_Genes_File_Exists = file.exists(significant_results_file),
      Summary_File_Exists = file.exists(summary_file),
      Significant_Genes_Rows = nrow(saved_significant_results),
      Summary_Total_Significant_Genes = saved_summary$Total_Significant_Genes,
      Summary_Total_Equals_Up_Down = saved_summary$Total_Significant_Genes ==
        saved_summary$Up + saved_summary$Down,
      stringsAsFactors = FALSE
    )
  )
}

cat("\nRunning limma differential expression analyses...\n")
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
  label = "limma analyses"
)

summary_list <- lapply(analysis_task_results, `[[`, "summary")
output_check_list <- lapply(analysis_task_results, `[[`, "output_check")

all_summary <- do.call(rbind, summary_list)
output_check <- do.call(rbind, output_check_list)

stopifnot(all(output_check$All_Genes_File_Exists))
stopifnot(all(output_check$Significant_Genes_File_Exists))
stopifnot(all(output_check$Summary_File_Exists))
stopifnot(
  all(output_check$Significant_Genes_Rows ==
      output_check$Summary_Total_Significant_Genes)
)
stopifnot(all(output_check$Summary_Total_Equals_Up_Down))

cat("\nAll analyses finished.\n")
cat("Output file check passed.\n\n")
print(all_summary[, SUMMARY_DISPLAY_COLUMNS], row.names = FALSE)
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
