# GSE114012转录因子富集分析
#
# 输入：
# 1. 01号脚本生成的每套差异分析显著基因列表：
#    results/ngs/GSE114012/tables/<analysis_name>/DEG/significant_genes.csv
# 2. 02号脚本生成的每套交集基因列表：
#    results/ngs/GSE114012/intersect/<intersection_scheme>/gene_list.csv
#
# 注意：
# - ORA/API类方法使用显著基因列表或交集gene_list作为离散输入；
# - VIPER/CollecTRI类activity方法使用对应差异分析的全量all_genes.csv作为连续ranked signature，
#   这是官方方法对连续分子读数/统计签名的要求，不使用交集后很短的deg_results.csv构建signature。
#
# 方法：
# - DoRothEA：
#   官方定位是signed TF-target regulon资源，不是单一统计检验。
#   这里按DoRothEA vignette推荐使用A/B/C置信度TF-target网络，
#   并用decoupleR::run_ora对显著基因列表做TF target ORA。
# - ChEA3：
#   使用rChEA3::queryChEA3调用ChEA3官方API。
#   ChEA3输入为HGNC gene symbol列表；输出保留官网返回的8个library/Integrated结果表，
#   包括Integrated--meanRank、Integrated--topRank、ENCODE/ReMap/Literature ChIP-seq、
#   ARCHS4/GTEx co-expression和Enrichr query来源。
# - VIPER：
#   使用DoRothEA signed regulon，结合Bioconductor viper::viper对DEG统计量签名推断TF活性。
#   输入不是离散基因列表，而是每个基因的t统计量；缺失t时回退到logFC。
# - ENRICHR：
#   使用enrichR::enrichr调用Ma'ayan Lab Enrichr API，默认运行TF相关数据库。
#   输入为显著基因symbol列表；输出保留enrichR官方返回的Term、P.value、Adjusted.P.value等列。
# - TRRUST：
#   使用OmnipathR::trrust_download获取TRRUST v2官方TF-target网络，
#   再用decoupleR::run_ora对显著基因列表做TF target ORA。
# - CollecTRI：
#   使用decoupleR::get_collectri官方接口获取CollecTRI网络；
#   若当前OmnipathR版本解析失败，则通过OmnipathR::static_tables定位官方静态TSV兜底读取。
#   统计部分使用decoupleR::run_viper基于连续DEG签名推断TF活性。
#
# 输出：
# results/ngs/GSE114012/TF/DEG/<analysis_name>/<method>/
# results/ngs/GSE114012/TF/intersect/<intersection_scheme>/<method>/
# 同步保存csv、md、tex三种表格格式。
# 文件名也带方法前缀，例如dorothea_deg_tf_enrichment.csv、
# chea3_deg_Integrated--meanRank.csv、collectri_intersect_tf_activity.csv。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
NETWORK_CACHE_FUNCTION_FILE <- "scripts/functions/network_cache_functions.R"
TF_FUNCTION_FILE <- "scripts/functions/tf_enrichment_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")
INTERSECT_ROOT <- file.path(RESULT_ROOT, "intersect")

# 测试时可临时设置TF_OUTPUT_ROOT，把结果写到temporary而不覆盖正式results。
OUTPUT_ROOT <- Sys.getenv("TF_OUTPUT_ROOT", unset = RESULT_ROOT)
TF_ROOT <- file.path(OUTPUT_ROOT, "TF")

# 远程API结果与网络TF资源缓存目录。
# ChEA3、Enrichr、TRRUST、CollecTRI等涉及网络访问的内容会优先读取本地缓存；
# 缓存缺失或超过7天后才重新获取，避免重复占用远程资源。
TF_REFERENCE_ROOT <- file.path("data", "reference", "TF")
TF_REFERENCE_MAX_AGE_DAYS <- 7
USE_TF_REFERENCE_CACHE <- TRUE

# 需要运行的转录因子富集方法。
# 可选："dorothea"、"chea3"、"viper"、"enrichr"、"trrust"、"collectri"。
METHODS_TO_RUN <- c("dorothea", "chea3", "viper", "enrichr", "trrust", "collectri")

# 需要纳入的输入类型。
# DEG      = 01号脚本输出的显著差异表达基因列表；
# INTERSECT = 02号脚本输出的交集基因列表。
INPUT_TYPES_TO_RUN <- c("DEG", "INTERSECT")

# 需要运行哪些差异分析设计。设为"all"时自动读取全部DEG/significant_genes.csv。
ANALYSES_TO_RUN <- "all"
# ANALYSES_TO_RUN <- c("DLD1", "HCT15", "SW48")

# 需要运行哪些交集方案。设为"all"时自动读取全部intersect/<scheme>/gene_list.csv。
INTERSECTION_SCHEMES_TO_RUN <- "all"
# INTERSECTION_SCHEMES_TO_RUN <- c("DLD1_HCT15", "DLD1_HCT15_SW48")

# 基因Symbol列。DoRothEA、ChEA3、VIPER当前均以gene symbol作为主要输入。
SYMBOL_COLUMN <- "Symbol"

# VIPER使用的差异分析排序统计量。
# t统计量同时包含效应方向和标准误信息，通常比单独logFC更适合构建连续签名。
RANK_METRIC_COLUMN <- "t"
FALLBACK_RANK_COLUMNS <- c("t", "logFC")

# DoRothEA regulon配置。
# A/B/C为官方推荐的高至中等置信度集合；若想提高覆盖率可加入"D"。
DOROTHEA_SPECIES <- "human"
DOROTHEA_CONFIDENCE_LEVELS <- c("A", "B", "C")

# decoupleR::run_ora参数，用于DoRothEA离散基因列表富集。
DOROTHEA_ORA_PARAMS <- list(
  n_background = 20000,
  with_ties = TRUE,
  seed = 20260605,
  minsize = 5
)

# TRRUST v2官方数据库配置。
TRRUST_SPECIES <- "human"
TRRUST_ORA_PARAMS <- list(
  n_background = 20000,
  with_ties = TRUE,
  seed = 20260605,
  minsize = 5
)

# viper::viper参数，用于DoRothEA regulon的VIPER TF活性推断。
# method="none"适合当前每个contrast只有一个ranked DEG signature的情况；
# method="scale"/"rank"/"mad"通常更适合有多列样本表达矩阵时做样本内标准化。
# minsize：每个TF至少需要多少个目标基因命中输入signature；
# eset.filter：输入已经是差异分析统计量签名，默认不再由VIPER过滤表达矩阵；
# nes=TRUE：输出normalized enrichment score，便于跨TF比较。
VIPER_PARAMS <- list(
  method = "none",
  minsize = 5,
  nes = TRUE,
  pleiotropy = FALSE,
  eset.filter = FALSE,
  verbose = FALSE
)

# CollecTRI官方网络配置。split_complexes=FALSE保留复合物名称，符合decoupleR默认建议。
COLLECTRI_SPECIES <- "human"
COLLECTRI_SPLIT_COMPLEXES <- FALSE
COLLECTRI_VIPER_PARAMS <- list(
  minsize = 5,
  verbose = FALSE,
  pleiotropy = TRUE,
  eset.filter = FALSE
)

# ChEA3官方API配置。
# ChEA3使用离散gene symbol列表进行TF over-representation/ranking。
CHEA3_API_URL <- "https://maayanlab.cloud/chea3/api/enrich/"
CHEA3_MIN_GENES <- 5

# ENRICHR官方API配置。默认选择TF相关数据库；可根据课题需要继续增删。
ENRICHR_SITE <- "Enrichr"
ENRICHR_DATABASES <- c(
  "ChEA_2022",
  "TRRUST_Transcription_Factors_2019",
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
  "ENCODE_TF_ChIP-seq_2015",
  "Transcription_Factor_PPIs"
)
ENRICHR_MIN_GENES <- 5
ENRICHR_BACKGROUND <- NULL
ENRICHR_INCLUDE_OVERLAP <- TRUE
ENRICHR_SLEEP_TIME <- 1

# ChEA3/ENRICHR是远程API方法。
# 这两类方法按官方普通调用方式逐个运行，并使用本地缓存减少真实远程请求；
# DoRothEA/VIPER/TRRUST/CollecTRI均为本地统计运算，继续使用parallel函数全速并行。
REMOTE_API_METHODS <- c("chea3", "enrichr")

# 表格预览行数。CSV保存完整结果；md/tex预览前若干行。
TABLE_PREVIEW_ROWS <- 21

# 重跑时清空本次方法对应的旧结果，避免新旧表格混合。
CLEAN_TF_OUTPUT_DIR <- TRUE

# 若缺少官方包，自动安装。若服务器/网络不稳定，可改为FALSE并手动安装。
AUTO_INSTALL_MISSING_PACKAGES <- TRUE

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 1. 安装/加载依赖与公共函数 --------------------------------------------------

install_tf_missing_packages <- function(auto_install = TRUE) {
  bioc_packages <- c("dorothea", "decoupleR", "viper")
  cran_packages <- c("httr", "jsonlite", "enrichR")
  runiverse_packages <- c("rChEA3")
  bioc_or_cran_packages <- c("OmnipathR")

  missing_bioc <- bioc_packages[
    !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  missing_cran <- cran_packages[
    !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  missing_runiverse <- runiverse_packages[
    !vapply(runiverse_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  missing_bioc_or_cran <- bioc_or_cran_packages[
    !vapply(bioc_or_cran_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  all_missing <- c(missing_bioc, missing_cran, missing_runiverse, missing_bioc_or_cran)
  if (length(all_missing) == 0L) {
    return(invisible(TRUE))
  }

  if (!auto_install) {
    stop(
      "Please install required R packages before running this script: ",
      paste(all_missing, collapse = ", ")
    )
  }

  if (length(missing_bioc) > 0L) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }

  if (length(missing_bioc_or_cran) > 0L) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing_bioc_or_cran, ask = FALSE, update = FALSE)
  }

  if (length(missing_cran) > 0L) {
    install.packages(missing_cran)
  }

  if (length(missing_runiverse) > 0L) {
    install.packages(
      missing_runiverse,
      repos = c("https://ckntav.r-universe.dev", "https://cloud.r-project.org")
    )
  }

  invisible(TRUE)
}

install_tf_missing_packages(AUTO_INSTALL_MISSING_PACKAGES)

source(REPORT_TABLE_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(NETWORK_CACHE_FUNCTION_FILE)
source(TF_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

run_quietly <- function(expr) {
  # 缓存命中或常规资源加载时不刷屏；真实错误仍会正常抛出。
  result <- NULL
  invisible(capture.output({
    result <- suppressMessages(expr)
  }))
  result
}


# 2. 准备输入、regulon和输出目录 ----------------------------------------------

runtime_methods <- tolower(get_runtime_vector("TF_TEST_METHODS", METHODS_TO_RUN))
runtime_methods <- intersect(
  runtime_methods,
  c("dorothea", "chea3", "viper", "enrichr", "trrust", "collectri")
)
if (length(runtime_methods) == 0L) {
  stop("No valid TF enrichment method was selected.")
}

runtime_input_types <- toupper(get_runtime_vector("TF_TEST_INPUT_TYPES", INPUT_TYPES_TO_RUN))
runtime_input_types <- intersect(runtime_input_types, c("DEG", "INTERSECT"))
if (length(runtime_input_types) == 0L) {
  stop("No valid TF enrichment input type was selected.")
}

all_inputs <- list()
if ("DEG" %in% runtime_input_types) {
  all_inputs <- c(
    all_inputs,
    get_tf_deg_inputs(TABLE_ROOT, analyses_to_run = ANALYSES_TO_RUN)
  )
}
if ("INTERSECT" %in% runtime_input_types) {
  all_inputs <- c(
    all_inputs,
    get_tf_intersection_inputs(
      INTERSECT_ROOT,
      schemes_to_run = INTERSECTION_SCHEMES_TO_RUN,
      table_root = TABLE_ROOT
    )
  )
}

test_max_inputs <- suppressWarnings(
  as.integer(Sys.getenv("TF_TEST_MAX_INPUTS", unset = "0"))
)
if (!is.na(test_max_inputs) && test_max_inputs > 0L) {
  all_inputs <- all_inputs[seq_len(min(test_max_inputs, length(all_inputs)))]
}

if (length(all_inputs) == 0L) {
  stop("No TF enrichment input files were found.")
}

if (CLEAN_TF_OUTPUT_DIR) {
  # 清理当前新目录结构，同时删除旧版TF/<method>/结构，避免历史结果混淆。
  clean_dirs <- c(
    file.path(TF_ROOT, c("DEG", "intersect", "summary")),
    file.path(TF_ROOT, runtime_methods)
  )
  for (clean_dir in unique(clean_dirs)) {
    if (dir.exists(clean_dir)) {
      unlink(clean_dir, recursive = TRUE, force = TRUE)
    }
  }
}

dorothea_network <- data.frame()
if (any(c("dorothea", "viper") %in% runtime_methods)) {
  dorothea_network <- run_quietly(
    load_dorothea_regulon(
      species = DOROTHEA_SPECIES,
      confidence_levels = DOROTHEA_CONFIDENCE_LEVELS
    )
  )
}

trrust_network <- data.frame()
if ("trrust" %in% runtime_methods) {
  trrust_network <- run_quietly(
    load_trrust_network_cached(
      species = TRRUST_SPECIES,
      reference_root = TF_REFERENCE_ROOT,
      max_age_days = TF_REFERENCE_MAX_AGE_DAYS,
      use_cache = USE_TF_REFERENCE_CACHE
    )
  )
}

collectri_network <- data.frame()
if ("collectri" %in% runtime_methods) {
  collectri_network <- run_quietly(
    load_collectri_network_cached(
      species = COLLECTRI_SPECIES,
      split_complexes = COLLECTRI_SPLIT_COMPLEXES,
      reference_root = TF_REFERENCE_ROOT,
      max_age_days = TF_REFERENCE_MAX_AGE_DAYS,
      use_cache = USE_TF_REFERENCE_CACHE
    )
  )
}

# 3. 构建并行任务 -------------------------------------------------------------

make_tf_tasks <- function(methods, inputs) {
  tasks <- list()
  task_index <- 0L

  for (method_name in methods) {
    for (input_index in seq_along(inputs)) {
      task_index <- task_index + 1L
      tasks[[task_index]] <- list(
        method = method_name,
        input = inputs[[input_index]]
      )
    }
  }

  tasks
}

tf_tasks <- make_tf_tasks(runtime_methods, all_inputs)
remote_task_ids <- which(vapply(
  tf_tasks,
  function(x) x$method %in% REMOTE_API_METHODS,
  logical(1)
))
local_task_ids <- setdiff(seq_along(tf_tasks), remote_task_ids)

cat("\nTF enrichment tasks prepared.\n")
cat("Local tasks:  ", length(local_task_ids), "\n", sep = "")
cat("Remote tasks: ", length(remote_task_ids), "\n", sep = "")

task_strategy <- setup_parallel_strategy(
  total_tasks = length(tf_tasks),
  inner_label = "TF method inner workers",
  nested_label = "Nested workers",
  print_strategy = FALSE
)

run_one_tf_task <- function(task_id) {
  run_tf_enrichment_task(
    task = tf_tasks[[task_id]],
    tf_root = TF_ROOT,
    dorothea_network = dorothea_network,
    trrust_network = trrust_network,
    collectri_network = collectri_network,
    symbol_column = SYMBOL_COLUMN,
    rank_metric_column = RANK_METRIC_COLUMN,
    fallback_rank_columns = FALLBACK_RANK_COLUMNS,
    dorothea_ora_params = DOROTHEA_ORA_PARAMS,
    trrust_ora_params = TRRUST_ORA_PARAMS,
    viper_params = VIPER_PARAMS,
    collectri_viper_params = COLLECTRI_VIPER_PARAMS,
    chea3_min_genes = CHEA3_MIN_GENES,
    chea3_api_url = CHEA3_API_URL,
    enrichr_databases = ENRICHR_DATABASES,
    enrichr_min_genes = ENRICHR_MIN_GENES,
    enrichr_background = ENRICHR_BACKGROUND,
    enrichr_include_overlap = ENRICHR_INCLUDE_OVERLAP,
    enrichr_sleep_time = ENRICHR_SLEEP_TIME,
    enrichr_site = ENRICHR_SITE,
    reference_root = TF_REFERENCE_ROOT,
    reference_max_age_days = TF_REFERENCE_MAX_AGE_DAYS,
    use_reference_cache = USE_TF_REFERENCE_CACHE,
    preview_rows = TABLE_PREVIEW_ROWS
  )
}


# 4. 并行运行TF富集分析 -------------------------------------------------------

summary_records <- list()

if (length(local_task_ids) > 0L) {
  cat("\nRunning local TF enrichment tasks...\n")
  local_results <- run_parallel_tasks_with_progress(
    task_ids = local_task_ids,
    task_function = run_one_tf_task,
    workers = min(task_strategy$task_workers, length(local_task_ids)),
    progress_label = "Local TF"
  )
  stop_on_parallel_errors(
    local_results,
    task_ids = local_task_ids,
    label = "local TF enrichment tasks"
  )
  summary_records <- c(summary_records, local_results)
}

if (length(remote_task_ids) > 0L) {
  cat("\nRunning remote API TF enrichment tasks with local cache...\n")
  remote_results <- run_parallel_tasks_with_progress(
    task_ids = remote_task_ids,
    task_function = run_one_tf_task,
    workers = 1L,
    progress_label = "Remote API"
  )
  stop_on_parallel_errors(
    remote_results,
    task_ids = remote_task_ids,
    label = "remote API TF enrichment tasks"
  )
  summary_records <- c(summary_records, remote_results)
}

summary_table <- do.call(rbind, summary_records)
rownames(summary_table) <- NULL
summary_table <- summary_table[order(
  summary_table$Method,
  summary_table$Input_Type,
  summary_table$Input_Name
), ]

summary_output_dir <- file.path(TF_ROOT, "summary")
dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
summary_csv_file <- write_csv_with_report_previews(
  dat = summary_table,
  csv_file = file.path(summary_output_dir, "summary.csv"),
  n_rows = TABLE_PREVIEW_ROWS
)


# 5. 终端简要汇总 -------------------------------------------------------------

cat("\nTF enrichment finished.\n")
cat("Result root:   ", TF_ROOT, "\n", sep = "")
cat("Summary table: ", summary_csv_file, "\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
