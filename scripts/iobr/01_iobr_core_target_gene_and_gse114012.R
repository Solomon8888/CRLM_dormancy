# IOBR核心快速分析：ATF3本地TCGA/GTEx + GSE114012表型变化
#
# 分析目的：
# 1. 研究指定基因在TCGA COAD 01A肿瘤样本中的TME、免疫/代谢/肿瘤signature功能关联；
# 2. 研究指定基因在GTEx COLON全部样本中的正常组织signature功能关联；
# 3. 复用GSE114012临床配置中的全部analysis_*设计，比较不同方案下IOBR表型变化；
# 4. 使用IOBR的deconvo_tme、calculate_sig_score、batch_cor、batch_wilcoxon/
#    batch_kruskal、sig_box、iobr_pca、cell_bar_plot等核心能力；
# 5. 结果统一保存到results/iobr/core，缓存与中间文件保存到data/iobr和temporary/iobr。


# 0. 项目定位与共用函数 --------------------------------------------------------

.iobr_get_current_script_file <- function() {
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

.iobr_find_project_root <- function(script_file = NA_character_) {
  start_points <- unique(c(
    getwd(),
    if (!is.na(script_file)) dirname(script_file) else character(0)
  ))

  for (start_point in start_points) {
    current <- normalizePath(start_point, winslash = "/", mustWork = TRUE)
    repeat {
      marker <- file.path(current, "scripts", "functions", "iobr_common_functions.R")
      if (file.exists(marker) && dir.exists(file.path(current, "data"))) {
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

SCRIPT_FILE <- .iobr_get_current_script_file()
PROJECT_ROOT <- .iobr_find_project_root(SCRIPT_FILE)
source(file.path(PROJECT_ROOT, "scripts", "functions", "iobr_common_functions.R"))
iobr_source_project_functions(PROJECT_ROOT)


# 1. 分析设计区：日常主要修改这里 --------------------------------------------

# 指定目标基因。默认示例为ATF3；多基因用英文逗号分隔。
TARGET_GENES <- iobr_parse_env_vector("IOBR_TARGET_GENES", c("ATF3"))

# TCGA主线：只分析COAD，样本筛选为01A primary tumor，不纳入READ。
TCGA_CANCERS <- iobr_parse_env_vector("IOBR_TCGA_CANCERS", c("COAD"))
TCGA_SAMPLE_DETAIL <- iobr_parse_env_vector("IOBR_TCGA_SAMPLE_DETAIL", c("tumor_primary_01a"))
TCGA_SAMPLE_PATTERN <- Sys.getenv("IOBR_TCGA_SAMPLE_PATTERN", unset = "-01A$")

# GTEx主线：默认使用COLON全部样本。
GTEX_TISSUES <- iobr_parse_env_vector("IOBR_GTEX_TISSUES", c("COLON"))

# GSE114012主线：默认使用全部样本，并按临床表中的analysis_*列逐一比较表型。
RUN_GSE114012 <- iobr_parse_env_logical("IOBR_RUN_GSE114012", TRUE)
GSE114012_SE_FILE <- file.path(PROJECT_ROOT, "data", "ngs", "GSE114012", "data_prepare", "GSE114012_se_raw.rds")
GSE114012_CLINICAL_FILE <- file.path(PROJECT_ROOT, "data", "ngs", "GSE114012", "data_prepare", "GSE114012_clinical_edit.csv")

# IOBR分析模块。deconvolution包含IOBR当前导出的全部TME估计算法；
# signature scoring默认使用IOBR支持的pca/zscore/ssgsea/integration四套算法。
RUN_DECONVOLUTION <- iobr_parse_env_logical("IOBR_RUN_DECONVOLUTION", TRUE)
RUN_GSE114012_DECONVOLUTION <- iobr_parse_env_logical("IOBR_RUN_GSE114012_DECONVOLUTION", FALSE)
DECONVOLUTION_METHODS <- iobr_parse_env_vector("IOBR_DECONVOLUTION_METHODS", iobr_get_tme_methods())
SIGNATURE_SCORE_METHODS <- iobr_parse_env_vector("IOBR_SIGNATURE_METHODS", iobr_get_signature_methods())

# signature范围。默认使用IOBR内置signature_collection全部signature；
# 测试时可用IOBR_MAX_SIGNATURES=30快速跑通。
SIGNATURE_GROUPS <- iobr_parse_env_vector("IOBR_SIGNATURE_GROUPS", character(0))
SIGNATURE_NAMES <- iobr_parse_env_vector("IOBR_SIGNATURE_NAMES", character(0))
MAX_SIGNATURES <- iobr_parse_env_integer("IOBR_MAX_SIGNATURES", Inf)
MINI_GENE_COUNT <- iobr_parse_env_integer("IOBR_MINI_GENE_COUNT", 3L)

# 表型统计和图片输出。
CORRELATION_METHOD <- tolower(Sys.getenv("IOBR_CORRELATION_METHOD", unset = "spearman"))
TOP_FEATURE_N <- iobr_parse_env_integer("IOBR_TOP_FEATURE_N", 15L)
CIBERSORT_PERM <- iobr_parse_env_integer("IOBR_CIBERSORT_PERM", 100L)
DRAW_ALL_SIGNATURE_HEATMAP <- iobr_parse_env_logical("IOBR_DRAW_ALL_SIGNATURE_HEATMAP", TRUE)
ALL_SIGNATURE_HEATMAP_MAX_FEATURES <- iobr_parse_env_integer("IOBR_ALL_SIGNATURE_HEATMAP_MAX_FEATURES", Inf)


# 2. 路径、缓存、并行与清理 ----------------------------------------------------

OUTPUT_ROOT <- file.path(PROJECT_ROOT, "results", "iobr", "core")
TEMP_ROOT <- file.path(PROJECT_ROOT, "temporary", "iobr", "core")
DATA_ROOT <- file.path(PROJECT_ROOT, "data", "iobr")
CACHE_ROOT <- file.path(DATA_ROOT, "prepared_inputs")
IOBR_CACHE_DIR <- file.path(DATA_ROOT, "IOBR_cache")

CLEAR_PREVIOUS_RUN_OUTPUTS <- iobr_parse_env_logical("IOBR_CLEAR_PREVIOUS", TRUE)
IOBR_VERBOSE <- iobr_parse_env_logical("IOBR_VERBOSE", FALSE)
DISABLE_FORK_PARALLEL <- iobr_parse_env_logical("IOBR_DISABLE_FORK", interactive())
PARALLEL_BACKEND <- Sys.getenv(
  "IOBR_PARALLEL_BACKEND",
  unset = Sys.getenv("PARALLEL_RUNTIME_BACKEND", unset = "auto")
)
MAX_PARALLEL_WORKERS <- iobr_parse_env_integer(
  "IOBR_PARALLEL_WORKERS",
  iobr_parse_env_integer("PARALLEL_RUNTIME_WORKERS", NA_integer_)
)

options(width = 200)
options(bitmapType = "cairo")
options(lifecycle_verbosity = "quiet")
options(parallel_runtime_force_single_line_progress = TRUE)
options(parallel_runtime_quiet_strategy = !IOBR_VERBOSE)
options(parallel_runtime_disable_fork = DISABLE_FORK_PARALLEL)
options(parallel_runtime_backend = PARALLEL_BACKEND)
options(parallel_runtime_worker_packages = c(
  "IOBR", "ggplot2", "dplyr", "SummarizedExperiment", "qs2", "Cairo"
))

iobr_prepare_output_tree(
  output_root = OUTPUT_ROOT,
  temporary_root = TEMP_ROOT,
  clear_previous = CLEAR_PREVIOUS_RUN_OUTPUTS
)
dir.create(DATA_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)

iobr_setup_runtime(
  project_root = PROJECT_ROOT,
  iobr_cache_dir = IOBR_CACHE_DIR,
  auto_install = TRUE,
  parallel_backend = PARALLEL_BACKEND,
  quiet_strategy = !IOBR_VERBOSE
)

suppressPackageStartupMessages({
  library(IOBR)
  library(ggplot2)
  library(dplyr)
  library(SummarizedExperiment)
})

SCRIPT_START_TIME <- start_runtime_timer()
PARALLEL_WORKERS <- if (is.na(MAX_PARALLEL_WORKERS)) {
  get_available_worker_count()
} else {
  max(1L, MAX_PARALLEL_WORKERS)
}


# 3. 写出IOBR能力目录 ---------------------------------------------------------

iobr_write_capability_catalog(OUTPUT_ROOT, scope = "core")
iobr_write_signature_catalog(OUTPUT_ROOT)

run_config <- data.frame(
  Key = c(
    "IOBR_version", "Target_genes", "TCGA_cancers", "GTEx_tissues",
    "Run_GSE114012", "Deconvolution_methods", "Signature_methods",
    "Signature_groups", "Signature_names", "Max_signatures",
    "Selected_signature_count", "Draw_all_signature_heatmap",
    "All_signature_heatmap_max_features", "Correlation_method",
    "Clear_previous_run_outputs", "Parallel_backend", "Parallel_workers",
    "Output_root", "Temporary_root", "Data_root"
  ),
  Value = c(
    as.character(utils::packageVersion("IOBR")),
    paste(TARGET_GENES, collapse = ", "),
    paste(TCGA_CANCERS, collapse = ", "),
    paste(GTEX_TISSUES, collapse = ", "),
    RUN_GSE114012,
    paste(DECONVOLUTION_METHODS, collapse = ", "),
    paste(SIGNATURE_SCORE_METHODS, collapse = ", "),
    paste(SIGNATURE_GROUPS, collapse = ", "),
    paste(SIGNATURE_NAMES, collapse = ", "),
    MAX_SIGNATURES,
    NA_character_,
    DRAW_ALL_SIGNATURE_HEATMAP,
    ALL_SIGNATURE_HEATMAP_MAX_FEATURES,
    CORRELATION_METHOD,
    CLEAR_PREVIOUS_RUN_OUTPUTS,
    PARALLEL_BACKEND,
    PARALLEL_WORKERS,
    OUTPUT_ROOT,
    TEMP_ROOT,
    DATA_ROOT
  ),
  stringsAsFactors = FALSE
)
iobr_write_module_csv(run_config, OUTPUT_ROOT, "run_summary", "000_iobr_core_run_configuration")

signatures <- iobr_select_signatures(
  selected_groups = SIGNATURE_GROUPS,
  selected_signatures = SIGNATURE_NAMES,
  max_signatures = MAX_SIGNATURES
)
if (length(signatures) == 0L) {
  stop("No IOBR signatures selected. Check IOBR_SIGNATURE_GROUPS/IOBR_SIGNATURE_NAMES.")
}

selected_signature_table <- data.frame(
  Signature = names(signatures),
  Gene_Count = vapply(signatures, length, integer(1)),
  stringsAsFactors = FALSE
)
iobr_write_module_csv(selected_signature_table, OUTPUT_ROOT, "run_summary", "001_selected_signatures")

run_config$Value[run_config$Key == "Selected_signature_count"] <- length(signatures)
iobr_write_module_csv(run_config, OUTPUT_ROOT, "run_summary", "000_iobr_core_run_configuration")


# 4. 准备并缓存输入数据 --------------------------------------------------------

input_manifest_list <- list()

for (target_gene in TARGET_GENES) {
  for (cancer in TCGA_CANCERS) {
    se_file <- file.path(PROJECT_ROOT, "data", "TCGA", cancer, "data_prepare", paste0(cancer, "_se_raw.rds"))
    input_manifest_list[[length(input_manifest_list) + 1L]] <- iobr_prepare_and_cache_local_input(
      project_root = PROJECT_ROOT,
      cache_root = CACHE_ROOT,
      dataset_id = paste0("TCGA_", cancer, "_01A"),
      se_file = se_file,
      target_gene = target_gene,
      project_id = cancer,
      assay_name = "tpm",
      deconvolution_assay = "tpm",
      score_assay = "tpm",
      score_log2 = TRUE,
      sample_detail_values = TCGA_SAMPLE_DETAIL,
      sample_detail_column = "group_detail",
      sample_barcode_pattern = TCGA_SAMPLE_PATTERN,
      clinical_file = file.path(PROJECT_ROOT, "data", "TCGA", cancer, "data_prepare", paste0(cancer, "_clinical_raw.csv")),
      tumor = TRUE,
      timer_indication = tolower(cancer)
    )
  }

  for (tissue in GTEX_TISSUES) {
    se_file <- file.path(PROJECT_ROOT, "data", "GTEx", tissue, "data_prepare", paste0(tissue, "_se_raw.rds"))
    input_manifest_list[[length(input_manifest_list) + 1L]] <- iobr_prepare_and_cache_local_input(
      project_root = PROJECT_ROOT,
      cache_root = CACHE_ROOT,
      dataset_id = paste0("GTEx_", tissue),
      se_file = se_file,
      target_gene = target_gene,
      project_id = paste0("GTEx_", tissue),
      assay_name = "tpm",
      deconvolution_assay = "tpm",
      score_assay = "tpm",
      score_log2 = TRUE,
      sample_detail_values = NULL,
      sample_barcode_pattern = NULL,
      clinical_file = file.path(PROJECT_ROOT, "data", "GTEx", tissue, "data_prepare", paste0(tissue, "_clinical_raw.csv")),
      tumor = FALSE,
      timer_indication = NULL
    )
  }

  if (RUN_GSE114012) {
    input_manifest_list[[length(input_manifest_list) + 1L]] <- iobr_prepare_and_cache_local_input(
      project_root = PROJECT_ROOT,
      cache_root = CACHE_ROOT,
      dataset_id = "GSE114012",
      se_file = GSE114012_SE_FILE,
      target_gene = target_gene,
      project_id = "GSE114012",
      assay_name = "tpm",
      deconvolution_assay = "tpm",
      score_assay = "tpm",
      score_log2 = TRUE,
      sample_detail_values = NULL,
      sample_barcode_pattern = NULL,
      clinical_file = GSE114012_CLINICAL_FILE,
      tumor = TRUE,
      timer_indication = NULL
    )
  }
}

input_manifest <- dplyr::bind_rows(input_manifest_list)
iobr_write_module_csv(input_manifest, OUTPUT_ROOT, "run_summary", "002_prepared_input_manifest")


# 5. 构建任务表 ----------------------------------------------------------------

task_list <- list()

if (RUN_DECONVOLUTION) {
  for (i in seq_len(nrow(input_manifest))) {
    dataset_id <- input_manifest$Dataset_ID[i]
    is_gse <- identical(dataset_id, "GSE114012")
    if (is_gse && !RUN_GSE114012_DECONVOLUTION) {
      next
    }

    methods <- DECONVOLUTION_METHODS
    if (is.na(input_manifest$TIMER_Indication[i]) || !nzchar(input_manifest$TIMER_Indication[i])) {
      methods <- setdiff(methods, "timer")
    }

    for (method in methods) {
      task_list[[length(task_list) + 1L]] <- data.frame(
        Task_Type = "deconvolution",
        Dataset_ID = dataset_id,
        Target_Gene = input_manifest$Target_Gene[i],
        Input_File = input_manifest$Cache_File[i],
        Method = method,
        Perm = CIBERSORT_PERM,
        Mini_Gene_Count = MINI_GENE_COUNT,
        Inner_Workers = 1L,
        stringsAsFactors = FALSE
      )
    }
  }
}

for (i in seq_len(nrow(input_manifest))) {
  for (method in SIGNATURE_SCORE_METHODS) {
    task_list[[length(task_list) + 1L]] <- data.frame(
      Task_Type = "signature_score",
      Dataset_ID = input_manifest$Dataset_ID[i],
      Target_Gene = input_manifest$Target_Gene[i],
      Input_File = input_manifest$Cache_File[i],
      Method = method,
      Perm = CIBERSORT_PERM,
      Mini_Gene_Count = MINI_GENE_COUNT,
      Inner_Workers = 1L,
      stringsAsFactors = FALSE
    )
  }
}

core_task_table <- dplyr::bind_rows(task_list)
core_task_table$Task_ID <- seq_len(nrow(core_task_table))
iobr_write_module_csv(core_task_table, OUTPUT_ROOT, "run_summary", "003_iobr_core_task_table")


# 6. 并行运行IOBR任务 ----------------------------------------------------------

run_one_iobr_core_task <- function(task_id) {
  task <- core_task_table[task_id, , drop = FALSE]

  if (identical(task$Task_Type, "deconvolution")) {
    return(iobr_run_deconvolution_task(
      task = task,
      output_root = OUTPUT_ROOT,
      correlation_method = CORRELATION_METHOD,
      top_n = TOP_FEATURE_N,
      draw_all_feature_heatmap = DRAW_ALL_SIGNATURE_HEATMAP,
      all_feature_heatmap_max_features = ALL_SIGNATURE_HEATMAP_MAX_FEATURES
    ))
  }

  if (identical(task$Task_Type, "signature_score")) {
    return(iobr_run_signature_score_task(
      task = task,
      output_root = OUTPUT_ROOT,
      signatures = signatures,
      correlation_method = CORRELATION_METHOD,
      top_n = TOP_FEATURE_N,
      draw_all_feature_heatmap = DRAW_ALL_SIGNATURE_HEATMAP,
      all_feature_heatmap_max_features = ALL_SIGNATURE_HEATMAP_MAX_FEATURES
    ))
  }

  iobr_status_row(
    dataset_id = task$Dataset_ID,
    module = task$Task_Type,
    task = task$Method,
    status = "failed",
    message = paste("Unknown task type:", task$Task_Type)
  )
}

parallel_strategy <- setup_parallel_strategy(
  total_tasks = nrow(core_task_table),
  max_workers = PARALLEL_WORKERS,
  inner_label = "IOBR inner workers",
  nested_label = "signature scoring workers"
)

task_results <- run_parallel_tasks_with_progress(
  task_ids = core_task_table$Task_ID,
  task_function = run_one_iobr_core_task,
  workers = parallel_strategy$task_workers,
  progress_label = "IOBR core"
)

status_table <- iobr_bind_rows(task_results)
if (nrow(status_table) == 0L) {
  status_table <- data.frame(
    Dataset_ID = character(0),
    Module = character(0),
    Task = character(0),
    Status = character(0),
    Message = character(0),
    Rows = integer(0),
    Columns = integer(0),
    Output_File = character(0)
  )
}
iobr_write_module_csv(status_table, OUTPUT_ROOT, "run_summary", "004_iobr_core_task_status")

runtime_seconds <- print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
runtime_table <- data.frame(
  Total_Runtime_Seconds = runtime_seconds,
  Total_Runtime = format_runtime_seconds(runtime_seconds),
  Finished_At = as.character(Sys.time()),
  Output_Root = OUTPUT_ROOT,
  stringsAsFactors = FALSE
)
iobr_write_module_csv(runtime_table, OUTPUT_ROOT, "run_summary", "005_iobr_core_runtime")

failed <- status_table[status_table$Status != "success", , drop = FALSE]
if (nrow(failed) > 0L) {
  iobr_write_module_csv(failed, OUTPUT_ROOT, "run_summary", "006_iobr_core_failed_tasks")
}

cat("\n01 IOBR core analysis finished: ", OUTPUT_ROOT, "\n", sep = "")
