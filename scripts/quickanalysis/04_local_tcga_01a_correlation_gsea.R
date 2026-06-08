# 本地TCGA 01A肿瘤样本：目标基因全基因相关性排序 + GSEA
#
# 设计目的：
# 1. 读取本项目已经准备好的TCGA SummarizedExperiment对象；
# 2. 仅提取01A原发肿瘤样本，且仅分析protein-coding基因；
# 3. 以指定基因为连续变量，计算其他全部protein-coding基因与它的相关性；
# 4. 相关性结果按相关系数从高到低排序保存；
# 5. 复用NGS流程中的clusterProfiler + msigdbr + GseaVis GSEA全套表格和绘图逻辑；
# 6. 结果统一保存到results/quickanalysis/local_tcga/<目标基因>_correlation。
#
# 默认示例：
# - 目标基因：ATF3
# - 癌种：COAD，不包含READ
# - 样本：TCGA COAD 01A primary solid tumor
# - 相关性方法：Spearman，可在脚本头部或环境变量切换Pearson


# 0. 定位项目根目录并加载quickanalysis局部helper ------------------------------

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
# TCGA_CORRELATION_TARGET_GENES=ATF3,MYC Rscript scripts/quickanalysis/04_local_tcga_01a_correlation_gsea.R
TARGET_GENES <- parse_env_vector("TCGA_CORRELATION_TARGET_GENES", c("ATF3"))

# 指定TCGA癌种。默认只分析COAD，不包含READ。
TARGET_CANCERS <- parse_env_vector("TCGA_CORRELATION_TARGET_CANCERS", c("COAD"))

# 本地SE文件模板。默认读取data/TCGA/<癌种>/data_prepare/<癌种>_se_raw.rds。
SE_FILE_TEMPLATE <- NULL

# 01A样本筛选。优先使用colData(se)$group_detail；若该列不存在，则回退到Sample_ID正则。
SAMPLE_DETAIL_FILTER <- parse_env_vector(
  "TCGA_CORRELATION_SAMPLE_DETAIL",
  c("tumor_primary_01a")
)
SAMPLE_BARCODE_PATTERN <- Sys.getenv(
  "TCGA_CORRELATION_SAMPLE_PATTERN",
  unset = "-01A$"
)

# 相关性分析使用的表达矩阵。默认使用log2(TPM + 1)。
EXPRESSION_ASSAY <- Sys.getenv("TCGA_CORRELATION_ASSAY", unset = "tpm")
EXPRESSION_LOG2 <- parse_env_logical("TCGA_CORRELATION_LOG2", TRUE)

# 相关性方法可选pearson或spearman。
CORRELATION_METHOD <- tolower(Sys.getenv("TCGA_CORRELATION_METHOD", unset = "spearman"))
stopifnot(CORRELATION_METHOD %in% c("pearson", "spearman"))

# 仅保留protein-coding基因；目标基因本身默认不纳入“其他基因”的相关性排序。
COUNT_ASSAY_NAME <- Sys.getenv("TCGA_CORRELATION_COUNT_ASSAY", unset = "counts")
GENE_BIOTYPE_FILTER <- "coding"
EXCLUDE_TARGET_GENE_FROM_RANKING <- parse_env_logical("TCGA_CORRELATION_EXCLUDE_TARGET", TRUE)
MIN_COMPLETE_SAMPLES <- parse_env_integer("TCGA_CORRELATION_MIN_COMPLETE_N", 10L)

# 相关性显著性筛选阈值。GSEA使用all_genes.csv中的Correlation排序，不受显著阈值限制。
P_VALUE_COLUMN <- Sys.getenv("TCGA_CORRELATION_P_VALUE_COLUMN", unset = "P.Value")
P_VALUE_CUTOFF <- as.numeric(Sys.getenv("TCGA_CORRELATION_P_VALUE_CUTOFF", unset = "0.05"))
CORRELATION_CUTOFF <- as.numeric(Sys.getenv("TCGA_CORRELATION_EFFECT_CUTOFF", unset = "0.3"))
TOP_CORRELATED_N <- parse_env_integer("TCGA_CORRELATION_TOP_N", 100L)
OUTPUT_DROP_COLUMNS <- c("Feature_ID", "Biotype", "Length")

# 单通路GSEA表达热图仍需要样本分组展示顺序；这里仅用于热图排序，不参与相关性计算。
HIGH_GROUP_LABEL <- Sys.getenv("TCGA_CORRELATION_HIGH_LABEL", unset = "High")
LOW_GROUP_LABEL <- Sys.getenv("TCGA_CORRELATION_LOW_LABEL", unset = "Low")


# 2. GSEA配置 -----------------------------------------------------------------

SPECIES <- Sys.getenv("TCGA_CORRELATION_SPECIES", unset = "human")
GENE_ID_TYPE <- Sys.getenv("TCGA_CORRELATION_GENE_ID_TYPE", unset = "ENTREZ")
RANK_METRIC_COLUMN <- Sys.getenv("TCGA_CORRELATION_RANK_COLUMN", unset = "Correlation")

GSEA_SIGNIFICANCE_COLUMN <- Sys.getenv("TCGA_CORRELATION_GSEA_P_COLUMN", unset = "pvalue")
GSEA_SIGNIFICANCE_CUTOFF <- as.numeric(Sys.getenv("TCGA_CORRELATION_GSEA_P_CUTOFF", unset = "0.05"))

GSEA_PARAMS <- list(
  exponent = 1,
  minGSSize = parse_env_integer("TCGA_CORRELATION_GSEA_MIN_GS_SIZE", 5L),
  maxGSSize = parse_env_integer("TCGA_CORRELATION_GSEA_MAX_GS_SIZE", 500L),
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

# 默认沿用NGS流程中的常用MSigDB类别。测试时可设为H快速跑通：
# TCGA_CORRELATION_GSEA_GENESETS=H Rscript scripts/quickanalysis/04_local_tcga_01a_correlation_gsea.R
GSEA_GENESETS_TO_RUN <- parse_env_vector(
  "TCGA_CORRELATION_GSEA_GENESETS",
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

READABLE_GENE_SYMBOLS <- parse_env_logical("TCGA_CORRELATION_READABLE_SYMBOLS", TRUE)
USE_QS2_CACHE <- parse_env_logical("TCGA_CORRELATION_USE_QS2_CACHE", TRUE)
REFRESH_QS2_CACHE <- parse_env_logical("TCGA_CORRELATION_REFRESH_QS2_CACHE", TRUE)
MSIGDB_REFERENCE_MAX_AGE_DAYS <- parse_env_integer("TCGA_CORRELATION_MSIGDB_MAX_AGE_DAYS", 7L)

# GSEA绘图配置：默认输出每类基因集的dotplot。
# 单通路GSEA/表达热图非常耗时，quickanalysis默认关闭；
# 如需完全复刻NGS 07号脚本的TopN单通路图，可设置TCGA_CORRELATION_DRAW_SINGLE_PATHWAY=1。
RUN_GSEA_PLOTS <- parse_env_logical("TCGA_CORRELATION_RUN_GSEA_PLOTS", TRUE)
DRAW_SINGLE_PATHWAY_GSEA <- parse_env_logical("TCGA_CORRELATION_DRAW_SINGLE_PATHWAY", FALSE)
SINGLE_PATHWAY_TOP_N <- parse_env_integer("TCGA_CORRELATION_SINGLE_PATHWAY_TOP_N", 10L)
SINGLE_PATHWAY_KEYWORDS <- parse_env_vector("TCGA_CORRELATION_SINGLE_PATHWAY_KEYWORDS", character(0))
SINGLE_PATHWAY_MAX_KEYWORD_TERMS <- parse_env_integer("TCGA_CORRELATION_SINGLE_PATHWAY_MAX_KEYWORD", 20L)
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

TPM_ASSAY_NAME <- EXPRESSION_ASSAY
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
SCRIPT_OUTPUT_PREFIX <- "tcga"
ANALYSIS_MODULE <- qa_make_gene_analysis_slug(TARGET_GENES, "correlation")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "results", "quickanalysis", DATASET_ID, ANALYSIS_MODULE)
TABLE_ROOT <- file.path(OUTPUT_ROOT, "tables")
TABLE_OUTPUT_ROOT <- TABLE_ROOT
CORRELATION_TABLE_ROOT <- file.path(TABLE_ROOT, "Correlation")
PLOT_ROOT <- file.path(OUTPUT_ROOT, "plots", "GSEA")
TEMP_ROOT <- file.path(PROJECT_ROOT, "temporary", "quickanalysis", DATASET_ID, ANALYSIS_MODULE, "tcga_correlation")
QS2_CACHE_DIR <- file.path(TEMP_ROOT, "GSEA_qs2_cache")
MSIGDB_REFERENCE_DIR <- file.path(PROJECT_ROOT, "data", "reference", "msigdb")

# 默认在运行前清空本脚本前次运行产生的全部结果与中间文件。
# 这只删除当前TCGA相关性分析名和tcga_前缀对应的旧文件，以及本脚本TEMP_ROOT；
# 不会删除同一ATF3_correlation目录中的GTEx结果、本地SE对象或MSigDB参考缓存。
CLEAR_PREVIOUS_RUN_OUTPUTS <- parse_env_logical("TCGA_CORRELATION_CLEAR_PREVIOUS", TRUE)
MAX_PARALLEL_WORKERS <- parse_env_integer("TCGA_CORRELATION_PARALLEL_WORKERS", NA_integer_)
options(width = 200)
options(lifecycle_verbosity = "quiet")
options(bitmapType = "cairo")


# 4. 加载包和项目共用函数 ------------------------------------------------------

qa_require_packages(c(
  "SummarizedExperiment",
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

suppressPackageStartupMessages(library(SummarizedExperiment))

qa_source_project_functions(PROJECT_ROOT)
SCRIPT_START_TIME <- start_runtime_timer()

PARALLEL_WORKERS <- if (is.na(MAX_PARALLEL_WORKERS)) {
  get_available_worker_count()
} else {
  max(1L, MAX_PARALLEL_WORKERS)
}


# 5. 清理旧结果并准备任务 ------------------------------------------------------

analysis_task_table <- expand.grid(
  Cancer = TARGET_CANCERS,
  Target_Gene = TARGET_GENES,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
analysis_task_table$Analysis_Name <- sanitize_file_name(
  paste(
    analysis_task_table$Cancer,
    "01A",
    analysis_task_table$Target_Gene,
    CORRELATION_METHOD,
    "correlation",
    sep = "_"
  )
)
analysis_task_table$Task_ID <- seq_len(nrow(analysis_task_table))

if (CLEAR_PREVIOUS_RUN_OUTPUTS) {
  qa_clean_quickanalysis_local_outputs(
    output_root = OUTPUT_ROOT,
    analysis_names = analysis_task_table$Analysis_Name,
    summary_prefix = SCRIPT_OUTPUT_PREFIX,
    table_categories = c("Correlation", "GSEA"),
    plot_categories = "GSEA",
    temp_root = TEMP_ROOT
  )
}
dir.create(CORRELATION_TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TABLE_ROOT, "GSEA"), recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)

cat("\nLocal TCGA 01A target-gene correlation + GSEA configuration:\n")
cat("Target genes: ", paste(TARGET_GENES, collapse = ", "), "\n", sep = "")
cat("Target cancers: ", paste(TARGET_CANCERS, collapse = ", "), "\n", sep = "")
cat("Sample filter: ", paste(SAMPLE_DETAIL_FILTER, collapse = ", "), "\n", sep = "")
cat("Gene biotype filter: protein-coding\n")
cat("Expression assay: ", EXPRESSION_ASSAY, if (EXPRESSION_LOG2) " log2(x + 1)" else "", "\n", sep = "")
cat("Correlation method: ", CORRELATION_METHOD, "\n", sep = "")
cat("GSEA rank column: ", RANK_METRIC_COLUMN, "\n", sep = "")
cat("Result root: ", OUTPUT_ROOT, "\n", sep = "")
cat("Prepared analyses: ", nrow(analysis_task_table), "\n", sep = "")


# 6. 单个目标基因/癌种相关性分析 ---------------------------------------------

qa_correlation_pvalues <- function(correlation, complete_n) {
  correlation <- as.numeric(correlation)
  complete_n <- as.integer(complete_n)

  p_value <- rep(NA_real_, length(correlation))
  valid <- is.finite(correlation) & complete_n >= 3L
  if (!any(valid)) {
    return(p_value)
  }

  r <- pmax(pmin(correlation[valid], 1 - .Machine$double.eps), -1 + .Machine$double.eps)
  t_stat <- r * sqrt((complete_n[valid] - 2) / pmax(1 - r^2, .Machine$double.eps))
  p_value[valid] <- 2 * stats::pt(-abs(t_stat), df = complete_n[valid] - 2)
  p_value[valid & abs(correlation[valid]) >= 1] <- 0
  p_value
}

run_one_correlation_analysis <- function(task_id) {
  cancer <- analysis_task_table$Cancer[task_id]
  target_gene <- analysis_task_table$Target_Gene[task_id]
  analysis_name <- analysis_task_table$Analysis_Name[task_id]

  inputs <- qa_load_tcga_01a_inputs(
    project_root = PROJECT_ROOT,
    cancer = cancer,
    se_file_template = SE_FILE_TEMPLATE,
    count_assay_name = COUNT_ASSAY_NAME,
    expression_assay_name = EXPRESSION_ASSAY,
    expression_log2_transform = EXPRESSION_LOG2,
    gene_biotype_filter = GENE_BIOTYPE_FILTER,
    sample_detail_filter = SAMPLE_DETAIL_FILTER,
    sample_barcode_pattern = SAMPLE_BARCODE_PATTERN
  )

  target_feature <- qa_find_target_feature(
    gene_annotation = inputs$gene_annotation,
    expression_matrix = inputs$expression,
    target_gene = target_gene
  )

  target_expression <- as.numeric(inputs$expression[target_feature$feature_id, ])
  names(target_expression) <- colnames(inputs$expression)
  valid_target_sample <- is.finite(target_expression)
  if (sum(valid_target_sample) < MIN_COMPLETE_SAMPLES) {
    stop("Too few samples have finite target-gene expression for ", analysis_name, ".")
  }

  sample_ids <- names(target_expression)[valid_target_sample]
  target_expression <- target_expression[valid_target_sample]
  expression_matrix <- inputs$expression[, sample_ids, drop = FALSE]
  sample_info <- inputs$sample_info[
    match(sample_ids, inputs$sample_info$Sample_ID),
    ,
    drop = FALSE
  ]
  rownames(sample_info) <- sample_info$Sample_ID

  target_median <- median(target_expression, na.rm = TRUE)
  heatmap_group <- ifelse(
    target_expression >= target_median,
    HIGH_GROUP_LABEL,
    LOW_GROUP_LABEL
  )
  design_info <- qa_make_analysis_design(
    sample_info = sample_info,
    analysis_name = analysis_name,
    group_values = heatmap_group,
    experiment_group = HIGH_GROUP_LABEL
  )

  complete_n <- rowSums(is.finite(expression_matrix))
  gene_sd <- apply(expression_matrix, 1, stats::sd, na.rm = TRUE)
  keep_gene <- complete_n >= MIN_COMPLETE_SAMPLES &
    is.finite(gene_sd) &
    gene_sd > 0

  if (EXCLUDE_TARGET_GENE_FROM_RANKING) {
    keep_gene[rownames(expression_matrix) == target_feature$feature_id] <- FALSE
  }

  if (sum(keep_gene) < 2L) {
    stop("Too few genes remained for correlation analysis: ", analysis_name)
  }

  expression_for_correlation <- expression_matrix[keep_gene, , drop = FALSE]
  correlation <- suppressWarnings(as.numeric(stats::cor(
    t(expression_for_correlation),
    target_expression,
    method = CORRELATION_METHOD,
    use = "pairwise.complete.obs"
  )))
  p_value <- qa_correlation_pvalues(
    correlation = correlation,
    complete_n = complete_n[keep_gene]
  )
  adj_p_value <- p.adjust(p_value, method = "BH")

  result_table <- data.frame(
    inputs$gene_annotation[keep_gene, , drop = FALSE],
    Target_Gene = target_gene,
    Target_Feature_ID = target_feature$feature_id,
    Correlation_Method = CORRELATION_METHOD,
    Correlation = correlation,
    P.Value = p_value,
    adj.P.Val = adj_p_value,
    Samples_Used = complete_n[keep_gene],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  result_table <- result_table[
    order(-result_table$Correlation, result_table$P.Value, result_table$Symbol),
    ,
    drop = FALSE
  ]

  significant_index <- is.finite(result_table$Correlation) &
    is.finite(result_table[[P_VALUE_COLUMN]]) &
    abs(result_table$Correlation) >= CORRELATION_CUTOFF &
    result_table[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF

  positive_index <- significant_index & result_table$Correlation > 0
  negative_index <- significant_index & result_table$Correlation < 0
  significant_table <- result_table[significant_index, , drop = FALSE]

  analysis_output_dir <- CORRELATION_TABLE_ROOT
  dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

  all_results_file <- write_csv_with_report_previews(
    qa_prepare_ranked_output_table(result_table, OUTPUT_DROP_COLUMNS),
    file.path(analysis_output_dir, paste0(analysis_name, "_all_genes.csv")),
    n_rows = 21,
    na = "NA"
  )
  significant_results_file <- write_csv_with_report_previews(
    qa_prepare_ranked_output_table(significant_table, OUTPUT_DROP_COLUMNS),
    file.path(analysis_output_dir, paste0(analysis_name, "_significant_genes.csv")),
    n_rows = 21,
    na = "NA"
  )
  top_positive_file <- write_csv_with_report_previews(
    qa_prepare_ranked_output_table(
      head(result_table[result_table$Correlation > 0, , drop = FALSE], TOP_CORRELATED_N),
      OUTPUT_DROP_COLUMNS
    ),
    file.path(analysis_output_dir, paste0(analysis_name, "_top_positive_genes.csv")),
    n_rows = 21,
    na = "NA"
  )
  top_negative_file <- write_csv_with_report_previews(
    qa_prepare_ranked_output_table(
      {
        negative_table <- result_table[result_table$Correlation < 0, , drop = FALSE]
        negative_table <- negative_table[
          order(negative_table$Correlation, negative_table$P.Value),
          ,
          drop = FALSE
        ]
        head(negative_table, TOP_CORRELATED_N)
      },
      OUTPUT_DROP_COLUMNS
    ),
    file.path(analysis_output_dir, paste0(analysis_name, "_top_negative_genes.csv")),
    n_rows = 21,
    na = "NA"
  )

  sample_expression_table <- data.frame(
    Cancer = cancer,
    Target_Gene = target_gene,
    Target_Feature_ID = target_feature$feature_id,
    Sample_ID = sample_info$Sample_ID,
    Patient_ID = if ("Patient_ID" %in% colnames(sample_info)) sample_info$Patient_ID else NA_character_,
    Expression_Assay = EXPRESSION_ASSAY,
    Expression_Log2_Transformed = EXPRESSION_LOG2,
    Target_Expression_Value = target_expression[sample_info$Sample_ID],
    Median_Cutoff = target_median,
    Median_Group_For_GSEA_Heatmap = heatmap_group[sample_info$Sample_ID],
    stringsAsFactors = FALSE
  )
  sample_expression_file <- write_csv_with_report_previews(
    sample_expression_table,
    file.path(analysis_output_dir, paste0(analysis_name, "_target_gene_sample_expression.csv")),
    n_rows = 21,
    na = "NA"
  )

  summary_table <- data.frame(
    Dataset = "TCGA",
    Cancer = cancer,
    Analysis_Name = analysis_name,
    Analysis_Type = "target_gene_correlation",
    Target_Gene = target_gene,
    Target_Feature_ID = target_feature$feature_id,
    Target_Symbol = if ("Symbol" %in% colnames(target_feature$annotation)) target_feature$annotation$Symbol else target_gene,
    Target_Entrez = if ("Entrez" %in% colnames(target_feature$annotation)) as.character(target_feature$annotation$Entrez) else NA_character_,
    Target_Candidate_Features = target_feature$candidate_count,
    Sample_Filter = inputs$sample_filter_source,
    Samples_Used = length(target_expression),
    Expression_Assay = EXPRESSION_ASSAY,
    Expression_Log2_Transformed = EXPRESSION_LOG2,
    Correlation_Method = CORRELATION_METHOD,
    P_Value_Method = "large-sample t approximation from correlation coefficient",
    Gene_Biotype_Filter = "protein_coding",
    Protein_Coding_Genes_Input = inputs$selected_gene_count,
    Ranked_Genes = nrow(result_table),
    Positive_Significant_Genes = sum(positive_index),
    Negative_Significant_Genes = sum(negative_index),
    Total_Significant_Genes = sum(significant_index),
    P_Value_Column = P_VALUE_COLUMN,
    P_Value_Cutoff = P_VALUE_CUTOFF,
    Correlation_Cutoff = CORRELATION_CUTOFF,
    Exclude_Target_Gene_From_Ranking = EXCLUDE_TARGET_GENE_FROM_RANKING,
    SE_RDS_File = inputs$se_file,
    All_Genes_File = all_results_file,
    Significant_Genes_File = significant_results_file,
    Top_Positive_Genes_File = top_positive_file,
    Top_Negative_Genes_File = top_negative_file,
    Target_Sample_Expression_File = sample_expression_file,
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
  stopifnot(file.exists(top_positive_file))
  stopifnot(file.exists(top_negative_file))
  stopifnot(file.exists(sample_expression_file))
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


# 7. 批量运行相关性分析 --------------------------------------------------------

parallel_strategy <- setup_parallel_strategy(
  total_tasks = nrow(analysis_task_table),
  max_workers = PARALLEL_WORKERS,
  inner_label = "correlation inner workers",
  nested_label = "GSEA workers"
)

cat("\nRunning target-gene correlation analyses...\n")
analysis_results <- run_parallel_tasks_with_progress(
  task_ids = analysis_task_table$Task_ID,
  task_function = run_one_correlation_analysis,
  workers = parallel_strategy$task_workers,
  progress_label = "Correlation"
)
stop_on_parallel_errors(
  analysis_results,
  task_ids = paste(analysis_task_table$Cancer, analysis_task_table$Target_Gene, sep = "_"),
  label = "target-gene correlation analyses"
)

analysis_summaries <- do.call(rbind, lapply(analysis_results, `[[`, "summary"))
ranked_file_info <- do.call(rbind, lapply(analysis_results, `[[`, "ranked_file_info"))
expression_input_list <- lapply(analysis_results, `[[`, "expression_input")
names(expression_input_list) <- ranked_file_info$Analysis_Name

run_summary_dir <- CORRELATION_TABLE_ROOT
dir.create(run_summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_csv <- write_csv_with_report_previews(
  analysis_summaries,
  file.path(run_summary_dir, paste0(SCRIPT_OUTPUT_PREFIX, "_correlation_summary.csv")),
  n_rows = 21,
  na = "NA"
)
ranked_input_csv <- write_csv_with_report_previews(
  ranked_file_info,
  file.path(run_summary_dir, paste0(SCRIPT_OUTPUT_PREFIX, "_gsea_ranked_input_files.csv")),
  n_rows = 21,
  na = "NA"
)

cat("\nCorrelation summary:\n")
print(
  analysis_summaries[
    ,
    c(
      "Cancer", "Target_Gene", "Samples_Used", "Ranked_Genes",
      "Positive_Significant_Genes", "Negative_Significant_Genes",
      "Total_Significant_Genes"
    )
  ],
  row.names = FALSE
)


# 8. 运行GSEA计算和绘图 --------------------------------------------------------

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

cat("\nOutput summary:\n")
cat("Correlation summary table: ", summary_csv, "\n", sep = "")
cat("GSEA input table:          ", ranked_input_csv, "\n", sep = "")
cat("GSEA summary table:        ", gsea_result$summary_csv_file, "\n", sep = "")
cat("Correlation table root:    ", CORRELATION_TABLE_ROOT, "\n", sep = "")
cat("GSEA table root:           ", file.path(TABLE_ROOT, "GSEA"), "\n", sep = "")
cat("GSEA plot root:            ", PLOT_ROOT, "\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")

cat("\nLocal TCGA target-gene correlation + GSEA analysis finished.\n")
