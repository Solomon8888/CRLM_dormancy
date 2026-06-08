# 本地GTEx样本：目标基因中位数分组差异分析 + GSEA
#
# 设计目的：
# 1. 读取本项目已经准备好的GTEx SummarizedExperiment对象；
# 2. 默认提取GTEx COLON全部样本，且仅分析protein-coding基因；
# 3. 以指定基因在全部GTEx样本中的表达中位数为阈值，分为高表达组和低表达组；
# 4. 复用NGS流程中的edgeR + voom + limma差异分析逻辑；
# 5. 复用NGS流程中的clusterProfiler + msigdbr + GseaVis GSEA全套表格和绘图逻辑；
# 6. 结果统一保存到results/quickanalysis/local_tcga/<目标基因>_deg。
#
# 默认示例：
# - 目标基因：ATF3
# - 组织：COLON
# - 样本：GTEx COLON全部样本


# 0. 定位项目根目录并加载共用helper ------------------------------------------

.local_get_current_script_file <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(normalizePath(frames[[i]]$ofile, winslash = "/", mustWork = TRUE))
    }
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }

  NA_character_
}

.local_find_project_root <- function(script_file = NA_character_) {
  start_points <- unique(c(
    getwd(),
    if (!is.na(script_file)) dirname(script_file) else character(0)
  ))

  for (start_point in start_points) {
    current <- normalizePath(start_point, winslash = "/", mustWork = TRUE)
    repeat {
      marker <- file.path(current, "scripts", "functions", "local_tcga_gsea_common_functions.R")
      if (file.exists(marker) && dir.exists(file.path(current, "scripts", "functions"))) {
        return(current)
      }
      parent <- dirname(current)
      if (identical(parent, current)) {
        break
      }
      current <- parent
    }
  }

  stop("Cannot locate project root. Please run this script from the project directory.")
}

SCRIPT_FILE <- .local_get_current_script_file()
PROJECT_ROOT <- .local_find_project_root(SCRIPT_FILE)
source(file.path(PROJECT_ROOT, "scripts", "functions", "local_tcga_gsea_common_functions.R"))


# 1. 分析设计区：日常主要修改这里 ---------------------------------------------

# 指定要分析的目标基因。默认示例为ATF3。
# 临时换基因可使用：
# GTEX_MEDIAN_DE_TARGET_GENES=ATF3,MYC Rscript scripts/quickanalysis/05_local_gtex_median_de_gsea.R
TARGET_GENES <- parse_env_vector("GTEX_MEDIAN_DE_TARGET_GENES", c("ATF3"))

# 指定GTEx组织数据集。默认分析COLON。
# 多组织可写成：GTEX_MEDIAN_DE_TARGET_TISSUES=COLON,LIVER
TARGET_TISSUES <- parse_env_vector("GTEX_MEDIAN_DE_TARGET_TISSUES", c("COLON"))

# 本地SE文件模板。默认读取data/GTEx/<组织>/data_prepare/<组织>_se_raw.rds。
SE_FILE_TEMPLATE <- NULL

# GTEx默认使用全部样本；如需筛选可配置样本注释列、取值或Sample_ID正则。
SAMPLE_FILTER_COLUMN <- Sys.getenv("GTEX_MEDIAN_DE_SAMPLE_FILTER_COLUMN", unset = "")
SAMPLE_FILTER_VALUES <- parse_env_vector("GTEX_MEDIAN_DE_SAMPLE_FILTER_VALUES", character(0))
SAMPLE_BARCODE_PATTERN <- Sys.getenv("GTEX_MEDIAN_DE_SAMPLE_PATTERN", unset = "")

# 目标基因分组使用的表达矩阵。TPM取log2(TPM + 1)后分组；单调变换不改变中位数分组方向。
TARGET_EXPRESSION_ASSAY <- Sys.getenv("GTEX_MEDIAN_DE_TARGET_ASSAY", unset = "tpm")
TARGET_EXPRESSION_LOG2 <- parse_env_logical("GTEX_MEDIAN_DE_TARGET_LOG2", TRUE)

# 差异分析使用raw counts，并仅保留protein-coding基因。
COUNT_ASSAY_NAME <- Sys.getenv("GTEX_MEDIAN_DE_COUNT_ASSAY", unset = "counts")
GENE_BIOTYPE_FILTER <- "coding"

# 中位数分组标签。阈值规则：表达值 >= median为High，表达值 < median为Low。
HIGH_GROUP_LABEL <- Sys.getenv("GTEX_MEDIAN_DE_HIGH_LABEL", unset = "High")
LOW_GROUP_LABEL <- Sys.getenv("GTEX_MEDIAN_DE_LOW_LABEL", unset = "Low")
MIN_SAMPLES_PER_GROUP <- parse_env_integer("GTEX_MEDIAN_DE_MIN_GROUP_N", 5L)

# limma显著差异筛选阈值。GSEA使用all_genes.csv中的t统计量，不受显著阈值限制。
P_VALUE_COLUMN <- Sys.getenv("GTEX_MEDIAN_DE_P_VALUE_COLUMN", unset = "P.Value")
P_VALUE_CUTOFF <- as.numeric(Sys.getenv("GTEX_MEDIAN_DE_P_VALUE_CUTOFF", unset = "0.05"))
LOGFC_CUTOFF <- as.numeric(Sys.getenv("GTEX_MEDIAN_DE_LOGFC_CUTOFF", unset = "0.5"))
OUTPUT_DROP_COLUMNS <- c("Feature_ID", "Biotype", "Length")


# 2. GSEA配置 -----------------------------------------------------------------

SPECIES <- Sys.getenv("GTEX_MEDIAN_DE_SPECIES", unset = "human")
GENE_ID_TYPE <- Sys.getenv("GTEX_MEDIAN_DE_GENE_ID_TYPE", unset = "ENTREZ")
RANK_METRIC_COLUMN <- Sys.getenv("GTEX_MEDIAN_DE_RANK_COLUMN", unset = "t")

GSEA_SIGNIFICANCE_COLUMN <- Sys.getenv("GTEX_MEDIAN_DE_GSEA_P_COLUMN", unset = "pvalue")
GSEA_SIGNIFICANCE_CUTOFF <- as.numeric(Sys.getenv("GTEX_MEDIAN_DE_GSEA_P_CUTOFF", unset = "0.05"))

GSEA_PARAMS <- list(
  exponent = 1,
  minGSSize = parse_env_integer("GTEX_MEDIAN_DE_GSEA_MIN_GS_SIZE", 5L),
  maxGSSize = parse_env_integer("GTEX_MEDIAN_DE_GSEA_MAX_GS_SIZE", 500L),
  pvalueCutoff = GSEA_SIGNIFICANCE_CUTOFF,
  pAdjustMethod = "BH",
  verbose = parse_env_logical("GTEX_MEDIAN_DE_VERBOSE", FALSE),
  nPerm = 1000,
  method = "multilevel",
  adaptive = FALSE,
  minPerm = 101,
  maxPerm = 1e5,
  pvalThreshold = 0.1
)

# 默认沿用NGS流程中的常用MSigDB类别。测试时可设为H快速跑通：
# GTEX_MEDIAN_DE_GSEA_GENESETS=H Rscript scripts/quickanalysis/05_local_gtex_median_de_gsea.R
GSEA_GENESETS_TO_RUN <- parse_env_vector(
  "GTEX_MEDIAN_DE_GSEA_GENESETS",
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

READABLE_GENE_SYMBOLS <- parse_env_logical("GTEX_MEDIAN_DE_READABLE_SYMBOLS", TRUE)
USE_QS2_CACHE <- parse_env_logical("GTEX_MEDIAN_DE_USE_QS2_CACHE", TRUE)
REFRESH_QS2_CACHE <- parse_env_logical("GTEX_MEDIAN_DE_REFRESH_QS2_CACHE", TRUE)
MSIGDB_REFERENCE_MAX_AGE_DAYS <- parse_env_integer("GTEX_MEDIAN_DE_MSIGDB_MAX_AGE_DAYS", 7L)

# GSEA绘图配置：默认输出每类基因集的dotplot。
# 单通路GSEA/表达热图非常耗时，quickanalysis默认关闭；
# 如需完全复刻NGS 07号脚本的TopN单通路图，可设置GTEX_MEDIAN_DE_DRAW_SINGLE_PATHWAY=1。
RUN_GSEA_PLOTS <- parse_env_logical("GTEX_MEDIAN_DE_RUN_GSEA_PLOTS", TRUE)
DRAW_SINGLE_PATHWAY_GSEA <- parse_env_logical("GTEX_MEDIAN_DE_DRAW_SINGLE_PATHWAY", FALSE)
SINGLE_PATHWAY_TOP_N <- parse_env_integer("GTEX_MEDIAN_DE_SINGLE_PATHWAY_TOP_N", 10L)
SINGLE_PATHWAY_KEYWORDS <- parse_env_vector("GTEX_MEDIAN_DE_SINGLE_PATHWAY_KEYWORDS", character(0))
SINGLE_PATHWAY_MAX_KEYWORD_TERMS <- parse_env_integer("GTEX_MEDIAN_DE_SINGLE_PATHWAY_MAX_KEYWORD", 20L)
SINGLE_PATHWAY_PVALUE_COLUMN <- GSEA_SIGNIFICANCE_COLUMN
SINGLE_PATHWAY_PVALUE_CUTOFF <- GSEA_SIGNIFICANCE_CUTOFF

SIMPLIFY_PATHWAY_PREFIX_IN_PLOT <- TRUE
REPLACE_UNDERSCORE_WITH_SPACE_IN_PLOT <- TRUE
GSEAVIS_DOTPLOT_PARAMS <- list(
  topn = 10,
  pval = if (GSEA_SIGNIFICANCE_COLUMN == "pvalue") GSEA_SIGNIFICANCE_CUTOFF else NULL,
  pajust = if (GSEA_SIGNIFICANCE_COLUMN == "p.adjust") GSEA_SIGNIFICANCE_CUTOFF else NULL,
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

TPM_ASSAY_NAME <- TARGET_EXPRESSION_ASSAY
SINGLE_PATHWAY_SAMPLE_LABEL_COLUMN <- "Sample_ID"
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


# 3. 输出、缓存与运行控制 ------------------------------------------------------

DATASET_ID <- "local_tcga"
SCRIPT_OUTPUT_PREFIX <- "gtex"
ANALYSIS_MODULE <- qa_make_gene_analysis_slug(TARGET_GENES, "deg")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "results", "quickanalysis", DATASET_ID, ANALYSIS_MODULE)
TABLE_ROOT <- file.path(OUTPUT_ROOT, "tables")
TABLE_OUTPUT_ROOT <- TABLE_ROOT
DEG_TABLE_ROOT <- file.path(TABLE_ROOT, "DEG")
PLOT_ROOT <- file.path(OUTPUT_ROOT, "plots", "GSEA")
VOLCANO_PLOT_ROOT <- file.path(OUTPUT_ROOT, "plots", "DEG")
TEMP_ROOT <- file.path(PROJECT_ROOT, "temporary", "quickanalysis", DATASET_ID, ANALYSIS_MODULE, "gtex_median_de")
QS2_CACHE_DIR <- file.path(TEMP_ROOT, "GSEA_qs2_cache")
MSIGDB_REFERENCE_DIR <- file.path(PROJECT_ROOT, "data", "reference", "msigdb")

# 默认在运行前清空本脚本前次运行产生的全部结果与中间文件。
# 这只删除当前GTEx差异分析名和gtex_前缀对应的旧文件，以及本脚本TEMP_ROOT；
# 不会删除同一ATF3_deg目录中的TCGA结果、本地SE对象或MSigDB参考缓存。
CLEAR_PREVIOUS_RUN_OUTPUTS <- parse_env_logical("GTEX_MEDIAN_DE_CLEAR_PREVIOUS", TRUE)
MAX_PARALLEL_WORKERS <- parse_env_integer(
  "GTEX_MEDIAN_DE_PARALLEL_WORKERS",
  parse_env_integer(
    "QUICKANALYSIS_PARALLEL_WORKERS",
    parse_env_integer("PARALLEL_RUNTIME_WORKERS", NA_integer_)
  )
)
QUICKANALYSIS_VERBOSE <- parse_env_logical("GTEX_MEDIAN_DE_VERBOSE", FALSE)
DISABLE_FORK_PARALLEL <- parse_env_logical("GTEX_MEDIAN_DE_DISABLE_FORK", interactive())
PARALLEL_BACKEND <- Sys.getenv(
  "GTEX_MEDIAN_DE_PARALLEL_BACKEND",
  unset = Sys.getenv(
    "QUICKANALYSIS_PARALLEL_BACKEND",
    unset = Sys.getenv("PARALLEL_RUNTIME_BACKEND", unset = "auto")
  )
)
options(width = 200)
options(lifecycle_verbosity = "quiet")
options(bitmapType = "cairo")
options(quickanalysis_verbose = QUICKANALYSIS_VERBOSE)
options(parallel_runtime_force_single_line_progress = TRUE)
options(parallel_runtime_quiet_strategy = !QUICKANALYSIS_VERBOSE)
options(parallel_runtime_disable_fork = DISABLE_FORK_PARALLEL)
options(parallel_runtime_backend = PARALLEL_BACKEND)


# 4. 加载包和项目共用函数 ------------------------------------------------------

qa_require_packages(c(
  "SummarizedExperiment",
  "limma",
  "edgeR",
  "clusterProfiler",
  "msigdbr",
  "GseaVis",
  "ggplot2",
  "qs2",
  "parallel",
  "Cairo",
  "RColorBrewer",
  "org.Hs.eg.db"
))

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
})

qa_source_project_functions(PROJECT_ROOT)
SCRIPT_START_TIME <- start_runtime_timer()

PARALLEL_WORKERS <- if (is.na(MAX_PARALLEL_WORKERS)) {
  get_available_worker_count()
} else {
  max(1L, MAX_PARALLEL_WORKERS)
}


# 5. 清理旧结果并准备任务 ------------------------------------------------------

analysis_task_table <- expand.grid(
  Tissue = TARGET_TISSUES,
  Target_Gene = TARGET_GENES,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
analysis_task_table$Analysis_Name <- sanitize_file_name(
  paste(
    analysis_task_table$Tissue,
    "GTEx",
    analysis_task_table$Target_Gene,
    "median_high_vs_low",
    sep = "_"
  )
)
analysis_task_table$Task_ID <- seq_len(nrow(analysis_task_table))

if (CLEAR_PREVIOUS_RUN_OUTPUTS) {
  qa_clean_quickanalysis_local_outputs(
    output_root = OUTPUT_ROOT,
    analysis_names = analysis_task_table$Analysis_Name,
    summary_prefix = SCRIPT_OUTPUT_PREFIX,
    table_categories = c("DEG", "GSEA"),
    plot_categories = c("DEG", "GSEA"),
    temp_root = TEMP_ROOT
  )
}
dir.create(DEG_TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TABLE_ROOT, "GSEA"), recursive = TRUE, showWarnings = FALSE)
dir.create(VOLCANO_PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)

sample_filter_text <- if (nzchar(SAMPLE_FILTER_COLUMN) && length(SAMPLE_FILTER_VALUES) > 0L) {
  paste0(SAMPLE_FILTER_COLUMN, " in ", paste(SAMPLE_FILTER_VALUES, collapse = ", "))
} else if (nzchar(SAMPLE_BARCODE_PATTERN)) {
  paste0("Sample_ID pattern ", SAMPLE_BARCODE_PATTERN)
} else {
  "all GTEx samples in each selected tissue"
}

if (QUICKANALYSIS_VERBOSE) {
  cat("\nLocal GTEx median-expression DE + GSEA configuration:\n")
  cat("Target genes: ", paste(TARGET_GENES, collapse = ", "), "\n", sep = "")
  cat("Target tissues: ", paste(TARGET_TISSUES, collapse = ", "), "\n", sep = "")
  cat("Sample filter: ", sample_filter_text, "\n", sep = "")
  cat("Gene biotype filter: protein-coding\n")
  cat("Grouping assay: ", TARGET_EXPRESSION_ASSAY, if (TARGET_EXPRESSION_LOG2) " log2(x + 1)" else "", "\n", sep = "")
  cat("DE assay: ", COUNT_ASSAY_NAME, "\n", sep = "")
  cat("GSEA rank column: ", RANK_METRIC_COLUMN, "\n", sep = "")
  cat("Result root: ", OUTPUT_ROOT, "\n", sep = "")
  cat("Prepared analyses: ", nrow(analysis_task_table), "\n", sep = "")
}


# 6. 单个目标基因/癌种差异分析 -------------------------------------------------

run_one_median_de_analysis <- function(task_id) {
  tissue <- analysis_task_table$Tissue[task_id]
  target_gene <- analysis_task_table$Target_Gene[task_id]
  analysis_name <- analysis_task_table$Analysis_Name[task_id]

  inputs <- qa_load_local_se_inputs(
    project_root = PROJECT_ROOT,
    data_source = "GTEx",
    dataset_id = tissue,
    se_file_template = SE_FILE_TEMPLATE,
    count_assay_name = COUNT_ASSAY_NAME,
    expression_assay_name = TARGET_EXPRESSION_ASSAY,
    expression_log2_transform = TARGET_EXPRESSION_LOG2,
    gene_biotype_filter = GENE_BIOTYPE_FILTER,
    sample_filter_column = SAMPLE_FILTER_COLUMN,
    sample_filter_values = SAMPLE_FILTER_VALUES,
    sample_filter_regex = SAMPLE_BARCODE_PATTERN,
    sample_filter_label = "GTEx configured samples"
  )

  target_feature <- qa_find_target_feature(
    gene_annotation = inputs$gene_annotation,
    expression_matrix = inputs$expression,
    target_gene = target_gene
  )

  target_expression <- as.numeric(inputs$expression[target_feature$feature_id, ])
  names(target_expression) <- colnames(inputs$expression)
  valid_sample <- is.finite(target_expression)
  if (sum(valid_sample) < MIN_SAMPLES_PER_GROUP * 2L) {
    stop("Too few samples have finite target-gene expression for ", analysis_name, ".")
  }

  sample_ids <- names(target_expression)[valid_sample]
  target_expression <- target_expression[valid_sample]
  median_cutoff <- median(target_expression, na.rm = TRUE)
  group_values <- ifelse(
    target_expression >= median_cutoff,
    HIGH_GROUP_LABEL,
    LOW_GROUP_LABEL
  )

  group_counts <- table(group_values)
  required_groups <- c(HIGH_GROUP_LABEL, LOW_GROUP_LABEL)
  if (!all(required_groups %in% names(group_counts))) {
    stop(
      "Median split did not produce both expression groups for ",
      analysis_name,
      ": ",
      paste(names(group_counts), group_counts, sep = "=", collapse = ", ")
    )
  }

  group_count_vector <- as.integer(group_counts[required_groups])
  names(group_count_vector) <- required_groups
  if (any(group_count_vector < MIN_SAMPLES_PER_GROUP)) {
    stop(
      "Median split produced too few samples for ",
      analysis_name,
      ": ",
      paste(names(group_counts), group_counts, sep = "=", collapse = ", ")
    )
  }

  sample_info <- inputs$sample_info[
    match(sample_ids, inputs$sample_info$Sample_ID),
    ,
    drop = FALSE
  ]
  rownames(sample_info) <- sample_info$Sample_ID
  counts <- inputs$counts[, sample_ids, drop = FALSE]
  expression_matrix <- inputs$expression[, sample_ids, drop = FALSE]

  design_info <- qa_make_analysis_design(
    sample_info = sample_info,
    analysis_name = analysis_name,
    group_values = group_values,
    experiment_group = HIGH_GROUP_LABEL
  )

  design_samples <- prepare_design_samples(
    sample_info = design_info$sample_info,
    group_column_index = design_info$analysis_designs$Column_Index,
    experiment_group = HIGH_GROUP_LABEL
  )
  sample_info_used <- design_samples$sample_info
  group_list <- design_samples$group_list
  control_group <- design_samples$control_group
  counts_used <- counts[, sample_info_used$Sample_ID, drop = FALSE]

  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- make.names(levels(group_list))
  rownames(design) <- colnames(counts_used)

  contrast_name <- paste0(HIGH_GROUP_LABEL, "_vs_", control_group)
  contrast_formula <- paste0(
    make.names(HIGH_GROUP_LABEL),
    " - ",
    make.names(control_group)
  )
  contrast_matrix <- makeContrasts(
    contrasts = contrast_formula,
    levels = design
  )
  colnames(contrast_matrix) <- contrast_name

  analysis_data <- prepare_ngs_data(
    counts = counts_used,
    gene_annotation = inputs$gene_annotation,
    group_list = group_list,
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
  stopifnot(P_VALUE_COLUMN %in% colnames(diff_results))

  regulation <- qa_count_regulation(
    dat = diff_results,
    p_value_column = P_VALUE_COLUMN,
    p_value_cutoff = P_VALUE_CUTOFF,
    effect_column = "logFC",
    effect_cutoff = LOGFC_CUTOFF
  )
  significant_results <- diff_results[regulation$significant, , drop = FALSE]

  analysis_output_dir <- DEG_TABLE_ROOT
  dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

  all_results_file <- write_csv_with_report_previews(
    qa_prepare_ranked_output_table(diff_results, OUTPUT_DROP_COLUMNS),
    file.path(analysis_output_dir, paste0(analysis_name, "_all_genes.csv")),
    n_rows = 21,
    na = "NA"
  )
  significant_results_file <- write_csv_with_report_previews(
    qa_prepare_ranked_output_table(significant_results, OUTPUT_DROP_COLUMNS),
    file.path(analysis_output_dir, paste0(analysis_name, "_significant_genes.csv")),
    n_rows = 21,
    na = "NA"
  )

  sample_group_table <- data.frame(
    Tissue = tissue,
    Target_Gene = target_gene,
    Target_Feature_ID = target_feature$feature_id,
    Sample_ID = sample_info_used$Sample_ID,
    Subject_ID = if ("Subject_ID" %in% colnames(sample_info_used)) sample_info_used$Subject_ID else NA_character_,
    Target_Expression_Assay = TARGET_EXPRESSION_ASSAY,
    Target_Expression_Log2_Transformed = TARGET_EXPRESSION_LOG2,
    Target_Expression_Value = target_expression[sample_info_used$Sample_ID],
    Median_Cutoff = median_cutoff,
    Expression_Group = as.character(group_list),
    stringsAsFactors = FALSE
  )
  sample_group_file <- write_csv_with_report_previews(
    sample_group_table,
    file.path(analysis_output_dir, paste0(analysis_name, "_sample_groups.csv")),
    n_rows = 21,
    na = "NA"
  )

  summary_table <- data.frame(
    Dataset = "GTEx",
    Tissue = tissue,
    Analysis_Name = analysis_name,
    Analysis_Type = "target_gene_median_split_differential_expression",
    Target_Gene = target_gene,
    Target_Feature_ID = target_feature$feature_id,
    Target_Symbol = if ("Symbol" %in% colnames(target_feature$annotation)) target_feature$annotation$Symbol else target_gene,
    Target_Entrez = if ("Entrez" %in% colnames(target_feature$annotation)) as.character(target_feature$annotation$Entrez) else NA_character_,
    Target_Candidate_Features = target_feature$candidate_count,
    Sample_Filter = inputs$sample_filter_source,
    Samples_Used = nrow(sample_info_used),
    High_Group = HIGH_GROUP_LABEL,
    Low_Group = control_group,
    High_Group_N = unname(group_count_vector[HIGH_GROUP_LABEL]),
    Low_Group_N = unname(group_count_vector[LOW_GROUP_LABEL]),
    Median_Cutoff = median_cutoff,
    Expression_Assay = TARGET_EXPRESSION_ASSAY,
    Expression_Log2_Transformed = TARGET_EXPRESSION_LOG2,
    Count_Assay = COUNT_ASSAY_NAME,
    Gene_Biotype_Filter = "protein_coding",
    Protein_Coding_Genes_Input = inputs$selected_gene_count,
    Low_Count_Filtered_Genes = analysis_data$filtered_genes,
    Ranked_Genes = nrow(diff_results),
    Contrast = contrast_name,
    Up = sum(regulation$up),
    Down = sum(regulation$down),
    Total_Significant_Genes = sum(regulation$significant),
    P_Value_Column = P_VALUE_COLUMN,
    P_Value_Cutoff = P_VALUE_CUTOFF,
    LogFC_Cutoff = LOGFC_CUTOFF,
    SE_RDS_File = inputs$se_file,
    All_Genes_File = all_results_file,
    Significant_Genes_File = significant_results_file,
    Sample_Group_File = sample_group_file,
    stringsAsFactors = FALSE
  )
  summary_file <- write_csv_with_report_previews(
    summary_table,
    file.path(analysis_output_dir, paste0(analysis_name, "_summary.csv")),
    n_rows = 21,
    na = "NA"
  )

  stopifnot(file.exists(all_results_file))
  stopifnot(file.exists(significant_results_file))
  stopifnot(file.exists(sample_group_file))
  stopifnot(file.exists(summary_file))

  list(
    summary = summary_table,
    ranked_file_info = data.frame(
      Analysis_Name = analysis_name,
      All_Genes_File = all_results_file,
      SE_RDS_File = inputs$se_file,
      stringsAsFactors = FALSE
    ),
    expression_input = list(
      expression_matrix = expression_matrix,
      gene_annotation = inputs$gene_annotation,
      analysis_designs = design_info$analysis_designs,
      sample_info = design_info$sample_info
    )
  )
}


# 7. 批量运行差异分析 ----------------------------------------------------------

parallel_strategy <- setup_parallel_strategy(
  total_tasks = nrow(analysis_task_table),
  max_workers = PARALLEL_WORKERS,
  inner_label = "limma inner workers",
  nested_label = "GSEA workers"
)

qa_log("\nRunning median-split differential expression analyses...\n")
analysis_results <- run_parallel_tasks_with_progress(
  task_ids = analysis_task_table$Task_ID,
  task_function = run_one_median_de_analysis,
  workers = parallel_strategy$task_workers,
  progress_label = "DE"
)
stop_on_parallel_errors(
  analysis_results,
  task_ids = paste(analysis_task_table$Tissue, analysis_task_table$Target_Gene, sep = "_"),
  label = "median-split DE analyses"
)

analysis_summaries <- do.call(rbind, lapply(analysis_results, `[[`, "summary"))
ranked_file_info <- do.call(rbind, lapply(analysis_results, `[[`, "ranked_file_info"))
expression_input_list <- lapply(analysis_results, `[[`, "expression_input")
names(expression_input_list) <- ranked_file_info$Analysis_Name

run_summary_dir <- DEG_TABLE_ROOT
dir.create(run_summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_csv <- write_csv_with_report_previews(
  analysis_summaries,
  file.path(run_summary_dir, paste0(SCRIPT_OUTPUT_PREFIX, "_median_de_summary.csv")),
  n_rows = 21,
  na = "NA"
)
ranked_input_csv <- write_csv_with_report_previews(
  ranked_file_info,
  file.path(run_summary_dir, paste0(SCRIPT_OUTPUT_PREFIX, "_gsea_ranked_input_files.csv")),
  n_rows = 21,
  na = "NA"
)

if (QUICKANALYSIS_VERBOSE) {
  cat("\nMedian-split DE summary:\n")
  print(
    analysis_summaries[
      ,
      c(
        "Tissue", "Target_Gene", "Samples_Used", "High_Group_N",
        "Low_Group_N", "Ranked_Genes", "Up", "Down",
        "Total_Significant_Genes"
      )
    ],
    row.names = FALSE
  )
}


# 8. 绘制传统火山图 ------------------------------------------------------------

volcano_result <- qa_run_traditional_volcano_plots(
  ranked_file_info = ranked_file_info,
  plot_root = VOLCANO_PLOT_ROOT,
  table_output_root = DEG_TABLE_ROOT,
  parallel_workers = PARALLEL_WORKERS,
  p_value_column = P_VALUE_COLUMN,
  p_value_cutoff = P_VALUE_CUTOFF,
  logfc_cutoff = LOGFC_CUTOFF,
  clean_outputs = TRUE,
  summary_file_prefix = SCRIPT_OUTPUT_PREFIX
)


# 9. 运行GSEA计算和绘图 --------------------------------------------------------

SE_RDS_FILE <- ranked_file_info$SE_RDS_File[1]

gsea_result <- qa_run_gsea_compute_and_plot(
  ranked_file_info = ranked_file_info,
  expression_input_list = expression_input_list,
  table_output_root = TABLE_OUTPUT_ROOT,
  plot_root = PLOT_ROOT,
  parallel_workers = PARALLEL_WORKERS,
  clean_outputs = TRUE,
  summary_file_prefix = SCRIPT_OUTPUT_PREFIX,
  run_plots = RUN_GSEA_PLOTS
)

if (QUICKANALYSIS_VERBOSE) {
  cat("\nOutput summary:\n")
  cat("DE summary table:    ", summary_csv, "\n", sep = "")
  cat("Volcano summary:     ", volcano_result$summary_csv_file, "\n", sep = "")
  cat("GSEA input table:    ", ranked_input_csv, "\n", sep = "")
  cat("GSEA summary table:  ", gsea_result$summary_csv_file, "\n", sep = "")
  cat("DEG table root:      ", DEG_TABLE_ROOT, "\n", sep = "")
  cat("GSEA table root:     ", file.path(TABLE_ROOT, "GSEA"), "\n", sep = "")
  cat("DEG plot root:       ", VOLCANO_PLOT_ROOT, "\n", sep = "")
  cat("GSEA plot root:      ", PLOT_ROOT, "\n", sep = "")
}

cat("\n05 local GTEx DEG + GSEA finished: ", OUTPUT_ROOT, "\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
