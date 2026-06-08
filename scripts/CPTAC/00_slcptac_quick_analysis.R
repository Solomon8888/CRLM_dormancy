# SLCPTAC快速CPTAC多组学全套分析脚本
#
# 设计目的：
# 1. 使用SLCPTAC包快速完成ATF3在CPTAC-COAD及泛癌中的多组学探索；
# 2. 覆盖SLCPTAC官方说明中的17个分析场景，并把每个场景登记到结果目录中的catalog表；
# 3. 正式结果统一保存到results/ngs/cptac，不让SLCPTAC自动生成的slcptac_output散落在项目根目录；
# 4. CPTAC bulk数据默认放在data/slcptac/bulk_data，并通过SL_BULK_DATA传给SLCPTAC；
# 5. SLCPTAC产生的临时图片、中间文件统一放在temporary/slcptac；
# 6. 复用scripts/functions中的表格保存、并行运行、进度条、耗时统计和缓存函数。
#
# 当前默认设计：
# - 主分析基因：ATF3
# - 主分析癌种：COAD
# - 泛癌分析：SLCPTAC支持的10个CPTAC癌种
# - COAD无Phospho/Methylation覆盖，因此这些任务默认记录为skipped；
#   如需展示磷酸化/甲基化场景，可在头部启用RUN_AUXILIARY_AVAILABLE_SCENARIOS。


# 0. 命令行覆盖配置的小工具 ---------------------------------------------------

parse_env_vector <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  parts <- unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

parse_env_logical <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(default)
  }

  value %in% c("1", "true", "t", "yes", "y")
}

parse_env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed)) {
    return(default)
  }

  parsed
}


# 1. 分析设计区：日常主要修改这里 ---------------------------------------------

# 主分析基因。临时换基因可使用：
# SLCPTAC_TARGET_GENES=ATF3,MYC Rscript scripts/CPTAC/00_slcptac_quick_analysis.R
TARGET_GENES <- parse_env_vector("SLCPTAC_TARGET_GENES", c("ATF3"))

# 主分析癌种。SLCPTAC中COAD仅代表CPTAC-COAD，不包含READ。
TARGET_CANCERS <- parse_env_vector("SLCPTAC_TARGET_CANCERS", c("COAD"))

# 泛癌分析癌种。默认使用SLCPTAC官方支持的全部10个CPTAC癌种。
PAN_CANCERS <- parse_env_vector(
  "SLCPTAC_PAN_CANCERS",
  c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")
)

# COAD没有Phospho/Methylation覆盖。若这里设为TRUE，会额外用AUXILIARY_CANCERS
# 跑一组SLCPTAC官方示例式的磷酸化/甲基化可用场景，方便快速看包的完整能力。
RUN_AUXILIARY_AVAILABLE_SCENARIOS <- parse_env_logical(
  "SLCPTAC_RUN_AUXILIARY_AVAILABLE_SCENARIOS",
  FALSE
)
AUXILIARY_CANCERS <- parse_env_vector("SLCPTAC_AUXILIARY_CANCERS", c("BRCA", "LUAD"))
AUXILIARY_PHOSPHO_GENE <- parse_env_vector("SLCPTAC_AUXILIARY_PHOSPHO_GENE", c("AKT1"))[1]

# 用于“多个连续变量”场景的候选基因。TARGET_GENES会自动放在第一位。
CONTINUOUS_PANEL_GENES <- unique(c(
  TARGET_GENES,
  parse_env_vector("SLCPTAC_CONTINUOUS_PANEL_GENES", c("JUN", "FOS", "DDIT3", "CEBPB"))
))

# 用于突变、共突变、突变驱动富集的候选基因。这里给COAD常见driver作为默认值。
MUTATION_PANEL_GENES <- parse_env_vector(
  "SLCPTAC_MUTATION_PANEL_GENES",
  c("APC", "TP53", "KRAS", "PIK3CA")
)
PRIMARY_MUTATION_GENE <- parse_env_vector("SLCPTAC_PRIMARY_MUTATION_GENE", c("KRAS"))[1]
SECONDARY_MUTATION_GENES <- parse_env_vector(
  "SLCPTAC_SECONDARY_MUTATION_GENES",
  c("TP53", "APC")
)

# 临床变量。SLCPTAC会按内部别名匹配，并把Age/BMI二分化。
CLINICAL_VARIABLES <- parse_env_vector(
  "SLCPTAC_CLINICAL_VARIABLES",
  c("Age", "Tumor_Stage", "Gender", "BMI")
)
PRIMARY_CLINICAL_VARIABLE <- parse_env_vector(
  "SLCPTAC_PRIMARY_CLINICAL_VARIABLE",
  c("Tumor_Stage")
)[1]

# 组学层与统计参数。
CORRELATION_METHOD <- Sys.getenv("SLCPTAC_CORRELATION_METHOD", unset = "pearson")
P_ADJUST_METHOD <- Sys.getenv("SLCPTAC_P_ADJUST_METHOD", unset = "BH")
ALPHA <- as.numeric(Sys.getenv("SLCPTAC_ALPHA", unset = "0.05"))
TOP_N <- parse_env_integer("SLCPTAC_TOP_N", 30L)
SURVIVAL_TYPES <- parse_env_vector("SLCPTAC_SURVIVAL_TYPES", c("OS", "PFS"))
SURVIVAL_CUTOFF_TYPE <- Sys.getenv("SLCPTAC_SURVIVAL_CUTOFF_TYPE", unset = "optimal")
SURVIVAL_MINPROP <- as.numeric(Sys.getenv("SLCPTAC_SURVIVAL_MINPROP", unset = "0.1"))
SURVIVAL_PERCENT <- as.numeric(Sys.getenv("SLCPTAC_SURVIVAL_PERCENT", unset = "0.25"))

# 富集数据库。默认只跑MsigDB Hallmark，保持“快速全套”。
# 如需更完整但更慢的路径数据库，可设置：
# SLCPTAC_ENRICH_DATABASES=MsigDB,GO,KEGG,Reactome
ENRICH_DATABASES <- parse_env_vector("SLCPTAC_ENRICH_DATABASES", c("MsigDB"))
GO_ONTOLOGIES <- parse_env_vector("SLCPTAC_GO_ONTOLOGIES", c("BP"))
GENOME_SCAN_MODALS <- parse_env_vector("SLCPTAC_GENOME_SCAN_MODALS", c("Protein"))
MSIGDB_CATEGORY <- Sys.getenv("SLCPTAC_MSIGDB_CATEGORY", unset = "H")
KEGG_CATEGORY <- Sys.getenv("SLCPTAC_KEGG_CATEGORY", unset = "pathway")
HGD_SOURCE <- Sys.getenv("SLCPTAC_HGDISEASE_SOURCE", unset = "do")
MESH_METHOD <- Sys.getenv("SLCPTAC_MESH_METHOD", unset = "gendoo")
MESH_CATEGORY <- Sys.getenv("SLCPTAC_MESH_CATEGORY", unset = "A")
ENRICHRDB_LIBRARY <- Sys.getenv(
  "SLCPTAC_ENRICHRDB_LIBRARY",
  unset = "Cancer_Cell_Line_Encyclopedia"
)

# 任务选择：
# - "all"：登记并运行本脚本封装的全部SLCPTAC场景；
# - "literature_quick_scan"：偏单基因文章常见的表达、临床、富集和生存；
# - 也可用逗号指定，例如：
#   SLCPTAC_ANALYSES=data_summary,scenario01,scenario13,scenario16 Rscript ...
ANALYSES_TO_RUN <- parse_env_vector("SLCPTAC_ANALYSES", "all")

# 是否额外生成“同一基因的可用组学层两两组合”相关/关联图。
RUN_MODAL_PAIR_GRID <- parse_env_logical("SLCPTAC_RUN_MODAL_PAIR_GRID", TRUE)

# 数据目录。SLCPTAC要求SL_BULK_DATA指向bulk根目录：
#   SL_BULK_DATA/
#     CPTAC_Omics_Split/TP53_cptac.qs
#     CPTAC_Omics_Split/<GENE>_cptac.qs
#     LinkedOmicsKB_PanCancer_Protein_Quantification.qs
#     LinkedOmicsKB_PanCancer_Clin.qs
#     ...
PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
SCRIPT_FILE <- file.path(PROJECT_ROOT, "scripts", "CPTAC", "00_slcptac_quick_analysis.R")
DATASET_ID <- "cptac"
DATA_TYPE <- "ngs"
RESULT_ROOT <- file.path(PROJECT_ROOT, "results", DATA_TYPE, DATASET_ID)
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "SLCPTAC")
TABLE_ROOT <- file.path(RESULT_ROOT, "tables", "SLCPTAC")
PLOT_PDF_DIR <- file.path(PLOT_ROOT, "pdf")
PLOT_PNG_DIR <- file.path(PLOT_ROOT, "png")
DATA_ROOT <- file.path(PROJECT_ROOT, "data", "slcptac")
TEMP_ROOT <- file.path(PROJECT_ROOT, "temporary", "slcptac")
SLCPTAC_BULK_DATA_ROOT <- Sys.getenv(
  "SL_BULK_DATA",
  unset = file.path(DATA_ROOT, "bulk_data")
)
SLCPTAC_BULK_DATA_ROOT <- Sys.getenv(
  "SLCPTAC_BULK_DATA",
  unset = SLCPTAC_BULK_DATA_ROOT
)
SLCPTAC_REFERENCE_CACHE_ROOT <- file.path(DATA_ROOT, "reference_cache")
SLCPTAC_TASK_CACHE_ROOT <- file.path(SLCPTAC_REFERENCE_CACHE_ROOT, "task_manifest")

# 结果与运行控制。
CLEAR_PREVIOUS_RUN_OUTPUTS <- parse_env_logical("SLCPTAC_CLEAR_PREVIOUS_OUTPUTS", TRUE)
CLEAN_TASK_OUTPUT_DIR <- parse_env_logical("SLCPTAC_CLEAN_OUTPUT", TRUE)
SAVE_RAW_DATA_TABLES <- parse_env_logical("SLCPTAC_SAVE_RAW_DATA_TABLES", TRUE)
SAVE_RESULT_OBJECTS <- parse_env_logical("SLCPTAC_SAVE_RESULT_OBJECTS", FALSE)
ADD_SCRIPT_TITLES <- parse_env_logical("SLCPTAC_ADD_SCRIPT_TITLES", TRUE)
STOP_IF_BULK_DATA_MISSING <- parse_env_logical("SLCPTAC_STOP_IF_DATA_MISSING", FALSE)

# 包安装控制。默认不在分析脚本中改动R库；如需自动尝试安装：
# SLCPTAC_AUTO_INSTALL=1 Rscript scripts/CPTAC/00_slcptac_quick_analysis.R
AUTO_INSTALL_SLCPTAC <- parse_env_logical("SLCPTAC_AUTO_INSTALL", FALSE)

# 并行配置。外层任务并行时，SLCPTAC内部GSEA worker会自动限制，避免过度抢核。
MAX_PARALLEL_WORKERS <- parse_env_integer("SLCPTAC_PARALLEL_WORKERS", NA_integer_)
ENRICHMENT_WORKERS <- parse_env_integer("SLCPTAC_ENRICHMENT_WORKERS", NA_integer_)
USE_SLCPTAC_TASK_CACHE <- parse_env_logical("SLCPTAC_USE_TASK_CACHE", TRUE)
SLCPTAC_TASK_CACHE_MAX_AGE_DAYS <- parse_env_integer("SLCPTAC_TASK_CACHE_MAX_AGE_DAYS", 30L)

options(width = 200)
options(bitmapType = "cairo")


# 2. 加载R包和项目共用函数 ----------------------------------------------------

install_slcptac_if_requested <- function() {
  if (!AUTO_INSTALL_SLCPTAC) {
    return(invisible(FALSE))
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop(
      "SLCPTAC is not installed, and package 'remotes' is unavailable. ",
      "Please install remotes first or install SLCPTAC manually."
    )
  }

  install_cran_if_missing <- function(packages, repos = getOption("repos")) {
    packages <- setdiff(packages, rownames(installed.packages()))
    if (length(packages) == 0) {
      return(invisible(TRUE))
    }
    install.packages(packages, repos = repos)
    invisible(TRUE)
  }

  ensure_qs_dependency <- function() {
    if (requireNamespace("qs", quietly = TRUE)) {
      return(invisible(TRUE))
    }

    # qs目前在部分新版R环境中无法从普通CRAN直接安装；按CRAN/r-universe/GitHub
    # 的顺序尝试，失败时给出明确原因，避免后续SLCPTAC安装报错过长。
    install_attempts <- list(
      function() install.packages("qs"),
      function() install.packages(
        "qs",
        repos = c(qsbase = "https://qsbase.r-universe.dev", CRAN = "https://cloud.r-project.org")
      ),
      function() remotes::install_github("qsbase/qs", upgrade = "never")
    )

    for (attempt in install_attempts) {
      try(suppressWarnings(attempt()), silent = TRUE)
      if (requireNamespace("qs", quietly = TRUE)) {
        return(invisible(TRUE))
      }
    }

    stop(
      "Dependency package 'qs' could not be installed. ",
      "SLCPTAC 1.2.0 imports qs, but qs 0.27.3 may fail to compile under R ",
      as.character(getRversion()),
      " because it uses R internal C API macros unavailable in this R build. ",
      "Please run this script in an R 4.4/4.5 environment, or install a local patched qs build first."
    )
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(
    c("fgsea", "ComplexHeatmap", "limma"),
    ask = FALSE,
    update = FALSE
  )

  install_cran_if_missing(
    c("dplyr", "ggplot2", "survival", "patchwork", "jsonlite",
      "Hmisc", "tidyr", "ggnewscale", "ggrepel", "scales", "circlize",
      "igraph", "ggraph", "Cairo", "png", "pdftools")
  )
  ensure_qs_dependency()

  remotes::install_github("GangLiLab/geneset", upgrade = "never")
  remotes::install_github("Zaoqu-Liu/genekitr2", upgrade = "never")
  remotes::install_github("SolvingLab/SLCPTAC", upgrade = "never")
  invisible(TRUE)
}

if (!requireNamespace("SLCPTAC", quietly = TRUE)) {
  install_slcptac_if_requested()
}

required_packages <- c("SLCPTAC", "ggplot2", "grid", "Cairo", "png")
missing_packages <- required_packages[
  !vapply(
    required_packages,
    function(package) {
      suppressPackageStartupMessages(
        suppressWarnings(requireNamespace(package, quietly = TRUE))
      )
    },
    logical(1)
  )
]
if (length(missing_packages) > 0) {
  stop(
    "Please install required R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    "\nOfficial install reference: remotes::install_github('SolvingLab/SLCPTAC')",
    "\nIf qs is the missing package under R >= 4.6, use R 4.4/4.5 or install a patched qs build first."
  )
}

suppressPackageStartupMessages({
  library(SLCPTAC)
})

PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
TABLE_IO_FUNCTION_FILE <- "scripts/functions/result_table_io_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
NETWORK_CACHE_FUNCTION_FILE <- "scripts/functions/network_cache_functions.R"

source(PLOTTING_FUNCTION_FILE)
source(TABLE_IO_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(NETWORK_CACHE_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

PARALLEL_WORKERS <- if (is.na(MAX_PARALLEL_WORKERS)) {
  get_available_worker_count()
} else {
  max(1L, MAX_PARALLEL_WORKERS)
}

Sys.setenv(SL_BULK_DATA = normalizePath(SLCPTAC_BULK_DATA_ROOT, winslash = "/", mustWork = FALSE))
Sys.setenv(SLCPTAC_CACHE_DIR = SLCPTAC_REFERENCE_CACHE_ROOT)
Sys.setenv(R_USER_CACHE_DIR = SLCPTAC_REFERENCE_CACHE_ROOT)


# 3. 通用工具函数 --------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    return(y)
  }
  x
}

safe_name <- function(x, default = "analysis") {
  sanitize_file_name(x, default = default)
}

paste_compact <- function(x, collapse = "_", default = "all") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    return(default)
  }
  paste(x, collapse = collapse)
}

write_table <- function(dat, output_dir, file_stem) {
  write_csv_with_report_previews(
    dat = dat,
    csv_file = file.path(output_dir, paste0(file_stem, ".csv")),
    na = "NA"
  )
}

format_task_timestamp <- function(time) {
  if (is.null(time) || is.na(time)) {
    return("")
  }
  format(time, "%Y-%m-%d %H:%M:%S")
}

capture_task <- function(expr) {
  warnings <- character(0)
  messages <- character(0)

  result <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      },
      message = function(message) {
        messages <<- c(messages, conditionMessage(message))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(error) error
  )

  list(
    result = result,
    warnings = unique(warnings),
    messages = unique(messages)
  )
}

object_dim_text <- function(x) {
  dims <- dim(x)
  if (is.null(dims)) {
    return(as.character(length(x)))
  }
  paste(dims, collapse = " x ")
}

list_to_data_frame <- function(x) {
  data.frame(
    Name = names(x) %||% seq_along(x),
    Class = vapply(x, function(value) paste(class(value), collapse = ";"), character(1)),
    Dimension = vapply(x, object_dim_text, character(1)),
    stringsAsFactors = FALSE
  )
}

normalize_result_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  for (column_name in colnames(dat)) {
    if (is.list(dat[[column_name]])) {
      dat[[column_name]] <- vapply(
        dat[[column_name]],
        function(x) paste(unlist(x), collapse = ","),
        character(1)
      )
    }
  }
  dat
}


# 4. SLCPTAC官方能力与数据覆盖 -------------------------------------------------

VALID_CANCERS <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")
VALID_MODALS <- c("RNAseq", "Protein", "Phospho", "Mutation", "Clinical", "logCNA", "Methylation", "Survival")

MODAL_AVAILABILITY <- list(
  RNAseq = VALID_CANCERS,
  Protein = VALID_CANCERS,
  Phospho = c("BRCA", "CCRCC", "GBM", "HNSCC", "LUAD", "LUSC", "PDAC", "UCEC"),
  Mutation = VALID_CANCERS,
  Clinical = VALID_CANCERS,
  logCNA = VALID_CANCERS,
  Methylation = c("CCRCC", "GBM", "HNSCC", "LUAD", "LUSC", "PDAC", "UCEC"),
  Survival = VALID_CANCERS
)

modal_available_for_cancers <- function(modal, cancers) {
  if (!modal %in% names(MODAL_AVAILABILITY)) {
    return(FALSE)
  }
  all(cancers %in% MODAL_AVAILABILITY[[modal]])
}

available_cancers_for_modal <- function(modal, cancers = PAN_CANCERS) {
  intersect(cancers, MODAL_AVAILABILITY[[modal]] %||% character(0))
}

make_modal_availability_table <- function() {
  rows <- list()
  for (modal in names(MODAL_AVAILABILITY)) {
    rows[[length(rows) + 1L]] <- data.frame(
      Modal = modal,
      Available_Cancers = paste(MODAL_AVAILABILITY[[modal]], collapse = ";"),
      Unavailable_In_PAN_CANCERS = paste(setdiff(PAN_CANCERS, MODAL_AVAILABILITY[[modal]]), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

LITERATURE_QUICK_SCAN_ANALYSES <- c(
  "data_summary",
  "modal_pair_grid",
  "scenario01",
  "scenario04",
  "scenario05",
  "scenario08",
  "scenario09",
  "scenario12",
  "scenario13",
  "scenario16",
  "scenario17"
)

ALL_SCENARIO_ANALYSES <- paste0("scenario", sprintf("%02d", 1:17))
ALL_WRAPPED_ANALYSES <- unique(c("data_summary", "modal_pair_grid", ALL_SCENARIO_ANALYSES))

make_analysis_catalog <- function(selected_analyses = character(0)) {
  rows <- list()
  add_row <- function(
      scenario,
      function_name,
      module,
      variable_design,
      statistics,
      plot_type,
      default_scope,
      notes = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      Script_Analysis = scenario,
      Function = function_name,
      Module = module,
      Variable_Design = variable_design,
      Statistics = statistics,
      Plot_Type = plot_type,
      Default_Scope = default_scope,
      Notes_CN = notes,
      Implemented_In_Script = TRUE,
      Selected_This_Run = scenario %in% selected_analyses,
      stringsAsFactors = FALSE
    )
  }

  add_row("data_summary", "cptac_load_modality", "数据审计", "包版本、数据根目录、组学覆盖、任务所需文件", "文件存在性检查", "CSV表", "本项目配置", "不调用SLCPTAC下载数据；只审计SL_BULK_DATA结构")
  add_row("modal_pair_grid", "cptac_correlation", "组学两两关联", "同一目标基因在可用组学层之间两两比较", "自动选择Pearson/Wilcoxon/Chi-square等", "Scatter/Lollipop/Box/Heatmap", "COAD及泛癌可用层", "覆盖7x7能力矩阵中的常用单基因组合")
  add_row("scenario01", "cptac_correlation", "相关/关联场景1", "1 continuous vs 1 continuous", "Pearson/Spearman/Kendall", "CorPlot或LollipopPlot", "ATF3 RNAseq vs Protein", "")
  add_row("scenario02", "cptac_correlation", "相关/关联场景2", "1 continuous vs multiple continuous", "相关分析+FDR", "LollipopPlot或DotPlot", "ATF3 Protein vs连续变量面板", "")
  add_row("scenario03", "cptac_correlation", "相关/关联场景3", "multiple continuous vs multiple continuous", "相关矩阵+FDR", "DotPlot", "连续变量面板RNAseq vs Protein", "")
  add_row("scenario04", "cptac_correlation", "相关/关联场景4", "1 categorical vs 1 continuous", "Wilcoxon或Kruskal-Wallis", "BoxPlot", "KRAS Mutation vs ATF3 Protein", "")
  add_row("scenario05", "cptac_correlation", "相关/关联场景5", "1 continuous vs multiple categorical", "多组Wilcoxon/Kruskal-Wallis", "Multiple BoxPlots", "ATF3 Protein vs突变面板", "")
  add_row("scenario06", "cptac_correlation", "相关/关联场景6", "multiple continuous vs 1 categorical", "多组Wilcoxon/Kruskal-Wallis", "Multiple BoxPlots", "连续变量面板Protein vs KRAS Mutation", "")
  add_row("scenario07", "cptac_correlation", "相关/关联场景7", "categorical vs categorical", "Fisher或Chi-square+Odds Ratio", "BarPlot或Heatmap", "共突变/临床-突变关联", "")
  add_row("scenario08", "cptac_enrichment", "富集场景8", "1 categorical vs genome-wide", "limma DEA", "NetworkPlot", "KRAS Mutation驱动蛋白组扫描", "")
  add_row("scenario09", "cptac_enrichment", "富集场景9", "1 categorical vs pathways", "limma DEA -> fgsea", "GSEA DotPlot", "KRAS Mutation驱动通路富集", "")
  add_row("scenario10", "cptac_enrichment", "富集场景10", "multiple categorical vs genome-wide", "多变量limma DEA", "DotPlot Paired", "突变面板驱动蛋白组扫描", "")
  add_row("scenario11", "cptac_enrichment", "富集场景11", "multiple categorical vs pathways", "多变量DEA -> fgsea", "GSEA Matrix", "突变面板通路富集", "")
  add_row("scenario12", "cptac_enrichment", "富集场景12", "1 continuous vs genome-wide", "相关排序", "NetworkPlot", "ATF3 Protein相关蛋白组扫描", "")
  add_row("scenario13", "cptac_enrichment", "富集场景13", "1 continuous vs pathways", "相关排序 -> fgsea", "GSEA DotPlot", "ATF3 Protein相关通路富集", "")
  add_row("scenario14", "cptac_enrichment", "富集场景14", "multiple continuous vs genome-wide", "多连续变量相关排序", "DotPlot Paired", "连续变量面板蛋白组扫描", "")
  add_row("scenario15", "cptac_enrichment", "富集场景15", "multiple continuous vs pathways", "多连续变量相关排序 -> fgsea", "GSEA Matrix", "连续变量面板通路富集", "")
  add_row("scenario16", "cptac_survival", "生存场景16", "1 variable vs survival", "KM + Cox", "KM+Cox combined", "ATF3单变量OS/PFS", "")
  add_row("scenario17", "cptac_survival", "生存场景17", "multiple variables vs survival", "Cox回归", "Forest plot", "ATF3泛癌或多基因OS/PFS", "")

  catalog <- do.call(rbind, rows)
  rownames(catalog) <- NULL
  catalog
}

write_analysis_catalog <- function(selected_analyses) {
  catalog <- make_analysis_catalog(selected_analyses)
  write_table(catalog, TABLE_ROOT, "000_slcptac_analysis_catalog")

  exported_functions <- sort(getNamespaceExports("SLCPTAC"))
  cataloged_functions <- unique(catalog$Function)
  missing_from_catalog <- setdiff(exported_functions, cataloged_functions)
  missing_table <- data.frame(
    Function = missing_from_catalog,
    Status = "SLCPTAC导出但不是单独分析场景函数，或已由上层场景间接覆盖",
    stringsAsFactors = FALSE
  )
  write_table(missing_table, TABLE_ROOT, "000_slcptac_exported_functions_not_in_catalog")

  invisible(list(catalog = catalog, missing = missing_table))
}

resolve_requested_analyses <- function() {
  requested <- ANALYSES_TO_RUN
  if (length(requested) == 1L && requested == "literature_quick_scan") {
    return(LITERATURE_QUICK_SCAN_ANALYSES)
  }
  if (length(requested) == 1L && requested == "all") {
    return(ALL_WRAPPED_ANALYSES)
  }

  unknown <- setdiff(requested, ALL_WRAPPED_ANALYSES)
  if (length(unknown) > 0L) {
    stop("Unknown SLCPTAC analyses requested: ", paste(unknown, collapse = ", "))
  }

  requested
}


# 5. 任务构建 ------------------------------------------------------------------

make_task_requirement <- function(modal, cancers, role = "var1") {
  data.frame(
    Role = role,
    Modal = modal,
    Cancers = paste(unique(cancers), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

make_task <- function(
    analysis,
    scenario_id,
    function_name,
    args,
    target,
    context,
    title,
    requirements,
    width = 6,
    height = 5) {
  list(
    analysis = analysis,
    scenario_id = scenario_id,
    function_name = function_name,
    args = args,
    target = target,
    context = context,
    title = title,
    requirements = requirements,
    width = width,
    height = height
  )
}

make_default_task_title <- function(task) {
  if (!is.null(task$title) && nzchar(task$title)) {
    return(task$title)
  }
  paste(gsub("_", " ", task$analysis, fixed = TRUE), task$target, task$context, sep = " | ")
}

estimate_title_extra_height <- function(title) {
  if (is.null(title) || !nzchar(title) || !ADD_SCRIPT_TITLES) {
    return(0)
  }
  line_count <- max(1L, ceiling(nchar(title, type = "width") / 70))
  0.30 + 0.18 * (line_count - 1L)
}

adjust_task_plot_size <- function(task) {
  task$height <- task$height + estimate_title_extra_height(make_default_task_title(task))

  if (task$scenario_id %in% c(7, 17)) {
    task$width <- max(task$width, 7.5)
    task$height <- max(task$height, 6.2)
  }
  if (task$scenario_id %in% c(8, 12)) {
    task$width <- max(task$width, 8.5)
    task$height <- max(task$height, 5.8)
  }
  if (task$scenario_id %in% c(9, 11, 13, 15)) {
    task$width <- max(task$width, 8.0)
    task$height <- max(task$height, 6.2)
  }
  if (task$scenario_id == 16) {
    task$width <- max(task$width, 8.5)
    task$height <- max(task$height, 4.8)
  }

  task
}

assign_task_indices <- function(tasks) {
  if (length(tasks) == 0L) {
    return(tasks)
  }
  for (i in seq_along(tasks)) {
    tasks[[i]]$task_id <- i
    tasks[[i]] <- adjust_task_plot_size(tasks[[i]])
  }
  tasks
}

task_file_stem <- function(task, suffix = NULL) {
  stem_parts <- c(
    sprintf("%03d", task$task_id %||% 0L),
    safe_name(task$analysis),
    safe_name(task$target),
    safe_name(task$context)
  )
  if (!is.null(suffix) && nzchar(suffix)) {
    stem_parts <- c(stem_parts, safe_name(suffix))
  }
  paste(stem_parts, collapse = "_")
}

add_correlation_task <- function(tasks, analysis, scenario_id, var1, var1_modal, var1_cancers,
                                 var2, var2_modal, var2_cancers, target, context, title,
                                 width = 6, height = 5) {
  requirements <- rbind(
    make_task_requirement(var1_modal, var1_cancers, "var1"),
    make_task_requirement(var2_modal, var2_cancers, "var2")
  )

  task <- make_task(
    analysis = analysis,
    scenario_id = scenario_id,
    function_name = "cptac_correlation",
    args = list(
      var1 = var1,
      var1_modal = var1_modal,
      var1_cancers = var1_cancers,
      var2 = var2,
      var2_modal = var2_modal,
      var2_cancers = var2_cancers,
      method = CORRELATION_METHOD,
      use = "pairwise.complete.obs",
      p_adjust_method = P_ADJUST_METHOD,
      alpha = ALPHA
    ),
    target = target,
    context = context,
    title = title,
    requirements = requirements,
    width = width,
    height = height
  )
  c(tasks, list(task))
}

add_enrichment_task <- function(tasks, analysis, scenario_id, var1, var1_modal, var1_cancers,
                                analysis_type, genome_modal, enrich_database, target,
                                context, title, width = 8, height = 6) {
  requirements <- rbind(
    make_task_requirement(var1_modal, var1_cancers, "var1"),
    make_task_requirement(genome_modal, var1_cancers, "genome")
  )

  task <- make_task(
    analysis = analysis,
    scenario_id = scenario_id,
    function_name = "cptac_enrichment",
    args = list(
      var1 = var1,
      var1_modal = var1_modal,
      var1_cancers = var1_cancers,
      analysis_type = analysis_type,
      enrich_database = enrich_database,
      enrich_ont = ifelse(enrich_database == "GO", GO_ONTOLOGIES[1], "BP"),
      genome_modal = genome_modal,
      method = CORRELATION_METHOD,
      top_n = TOP_N,
      n_workers = 1L,
      kegg_category = KEGG_CATEGORY,
      msigdb_category = MSIGDB_CATEGORY,
      hgdisease_source = HGD_SOURCE,
      mesh_method = MESH_METHOD,
      mesh_category = MESH_CATEGORY,
      enrichrdb_library = ENRICHRDB_LIBRARY
    ),
    target = target,
    context = context,
    title = title,
    requirements = requirements,
    width = width,
    height = height
  )
  c(tasks, list(task))
}

add_survival_task <- function(tasks, analysis, scenario_id, var1, var1_modal, var1_cancers,
                              surv_type, target, context, title, width = 8.5, height = 5) {
  requirements <- rbind(
    make_task_requirement(var1_modal, var1_cancers, "var1"),
    make_task_requirement("Survival", var1_cancers, "survival")
  )

  task <- make_task(
    analysis = analysis,
    scenario_id = scenario_id,
    function_name = "cptac_survival",
    args = list(
      var1 = var1,
      var1_modal = var1_modal,
      var1_cancers = var1_cancers,
      surv_type = surv_type,
      cutoff_type = SURVIVAL_CUTOFF_TYPE,
      minprop = SURVIVAL_MINPROP,
      percent = SURVIVAL_PERCENT,
      show_cindex = TRUE
    ),
    target = target,
    context = context,
    title = title,
    requirements = requirements,
    width = width,
    height = height
  )
  c(tasks, list(task))
}

build_modal_pair_grid_tasks <- function(selected_analyses) {
  if (!RUN_MODAL_PAIR_GRID || !"modal_pair_grid" %in% selected_analyses) {
    return(list())
  }

  tasks <- list()
  continuous_modals <- c("RNAseq", "Protein", "logCNA", "Methylation")
  categorical_modals <- c("Mutation")

  for (gene in TARGET_GENES) {
    for (cancer in TARGET_CANCERS) {
      available_cont <- continuous_modals[
        vapply(continuous_modals, modal_available_for_cancers, logical(1), cancers = cancer)
      ]
      if (length(available_cont) >= 2L) {
        pairs <- utils::combn(available_cont, 2, simplify = FALSE)
        for (pair in pairs) {
          tasks <- add_correlation_task(
            tasks,
            analysis = "modal_pair_grid",
            scenario_id = 1L,
            var1 = gene,
            var1_modal = pair[1],
            var1_cancers = cancer,
            var2 = gene,
            var2_modal = pair[2],
            var2_cancers = cancer,
            target = gene,
            context = paste(cancer, pair[1], "vs", pair[2], sep = "_"),
            title = paste0(gene, " ", pair[1], " vs ", pair[2], " in CPTAC-", cancer),
            width = 5.4,
            height = 4.8
          )
        }
      }

      for (cat_modal in categorical_modals) {
        if (!modal_available_for_cancers(cat_modal, cancer)) {
          next
        }
        for (con_modal in available_cont) {
          tasks <- add_correlation_task(
            tasks,
            analysis = "modal_pair_grid",
            scenario_id = 4L,
            var1 = gene,
            var1_modal = cat_modal,
            var1_cancers = cancer,
            var2 = gene,
            var2_modal = con_modal,
            var2_cancers = cancer,
            target = gene,
            context = paste(cancer, cat_modal, "vs", con_modal, sep = "_"),
            title = paste0(gene, " ", cat_modal, " vs ", con_modal, " in CPTAC-", cancer),
            width = 5.2,
            height = 4.8
          )
        }
      }

      for (clinical_var in CLINICAL_VARIABLES) {
        for (con_modal in available_cont) {
          tasks <- add_correlation_task(
            tasks,
            analysis = "modal_pair_grid",
            scenario_id = 4L,
            var1 = clinical_var,
            var1_modal = "Clinical",
            var1_cancers = cancer,
            var2 = gene,
            var2_modal = con_modal,
            var2_cancers = cancer,
            target = gene,
            context = paste(cancer, clinical_var, "vs", con_modal, sep = "_"),
            title = paste0(clinical_var, " vs ", gene, " ", con_modal, " in CPTAC-", cancer),
            width = 5.2,
            height = 4.8
          )
        }
      }
    }
  }

  tasks
}

build_slcptac_tasks <- function(selected_analyses) {
  tasks <- list()
  main_gene <- TARGET_GENES[1]
  target_cancer <- TARGET_CANCERS[1]
  pan_common <- PAN_CANCERS
  pan_phospho <- available_cancers_for_modal("Phospho", PAN_CANCERS)
  pan_methylation <- available_cancers_for_modal("Methylation", PAN_CANCERS)

  tasks <- c(tasks, build_modal_pair_grid_tasks(selected_analyses))

  for (gene in TARGET_GENES) {
    for (cancer in TARGET_CANCERS) {
      if ("scenario01" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario01", 1L,
          gene, "RNAseq", cancer,
          gene, "Protein", cancer,
          gene, cancer,
          paste0(gene, " RNAseq vs Protein in CPTAC-", cancer),
          width = 5.4, height = 4.8
        )
        tasks <- add_correlation_task(
          tasks, "scenario01", 1L,
          gene, "RNAseq", pan_common,
          gene, "Protein", pan_common,
          gene, "pan_cancer_RNAseq_vs_Protein",
          paste0(gene, " RNAseq vs Protein across CPTAC pan-cancer"),
          width = 8, height = 4.5
        )
      }

      if ("scenario02" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario02", 2L,
          gene, "Protein", cancer,
          CONTINUOUS_PANEL_GENES, "RNAseq", cancer,
          gene, paste0(cancer, "_Protein_vs_panel_RNAseq"),
          paste0(gene, " Protein vs RNAseq panel in CPTAC-", cancer),
          width = 7, height = 4.8
        )
      }

      if ("scenario03" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario03", 3L,
          CONTINUOUS_PANEL_GENES, "RNAseq", cancer,
          CONTINUOUS_PANEL_GENES, "Protein", cancer,
          gene, paste0(cancer, "_panel_RNAseq_vs_Protein"),
          paste0("RNAseq vs Protein panel correlation in CPTAC-", cancer),
          width = 7, height = 5.6
        )
      }

      if ("scenario04" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario04", 4L,
          PRIMARY_MUTATION_GENE, "Mutation", cancer,
          gene, "Protein", cancer,
          gene, paste0(cancer, "_", PRIMARY_MUTATION_GENE, "_Mutation_vs_Protein"),
          paste0(PRIMARY_MUTATION_GENE, " Mutation vs ", gene, " Protein in CPTAC-", cancer),
          width = 5.2, height = 4.8
        )
        tasks <- add_correlation_task(
          tasks, "scenario04", 4L,
          PRIMARY_CLINICAL_VARIABLE, "Clinical", cancer,
          gene, "Protein", cancer,
          gene, paste0(cancer, "_", PRIMARY_CLINICAL_VARIABLE, "_vs_Protein"),
          paste0(PRIMARY_CLINICAL_VARIABLE, " vs ", gene, " Protein in CPTAC-", cancer),
          width = 5.2, height = 4.8
        )
      }

      if ("scenario05" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario05", 5L,
          gene, "Protein", cancer,
          MUTATION_PANEL_GENES, "Mutation", cancer,
          gene, paste0(cancer, "_Protein_vs_mutation_panel"),
          paste0(gene, " Protein vs mutation panel in CPTAC-", cancer),
          width = 9, height = 4.8
        )
      }

      if ("scenario06" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario06", 6L,
          CONTINUOUS_PANEL_GENES, "Protein", cancer,
          PRIMARY_MUTATION_GENE, "Mutation", cancer,
          gene, paste0(cancer, "_protein_panel_vs_", PRIMARY_MUTATION_GENE, "_Mutation"),
          paste0("Protein panel vs ", PRIMARY_MUTATION_GENE, " Mutation in CPTAC-", cancer),
          width = 9, height = 4.8
        )
      }

      if ("scenario07" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario07", 7L,
          unique(c(PRIMARY_MUTATION_GENE, SECONDARY_MUTATION_GENES)), "Mutation", cancer,
          setdiff(MUTATION_PANEL_GENES, PRIMARY_MUTATION_GENE), "Mutation", cancer,
          gene, paste0(cancer, "_comutation"),
          paste0("Co-mutation landscape in CPTAC-", cancer),
          width = 7.5, height = 5.8
        )
        tasks <- add_correlation_task(
          tasks, "scenario07", 7L,
          PRIMARY_CLINICAL_VARIABLE, "Clinical", cancer,
          PRIMARY_MUTATION_GENE, "Mutation", cancer,
          gene, paste0(cancer, "_clinical_vs_mutation"),
          paste0(PRIMARY_CLINICAL_VARIABLE, " vs ", PRIMARY_MUTATION_GENE, " Mutation in CPTAC-", cancer),
          width = 6.5, height = 5.2
        )
      }

      for (genome_modal in GENOME_SCAN_MODALS) {
        if ("scenario08" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario08", 8L,
            PRIMARY_MUTATION_GENE, "Mutation", cancer,
            "genome", genome_modal, "MsigDB",
            PRIMARY_MUTATION_GENE, paste0(cancer, "_Mutation_vs_", genome_modal, "_genome"),
            paste0(PRIMARY_MUTATION_GENE, " Mutation-driven ", genome_modal, " genome scan in CPTAC-", cancer),
            width = 8.5, height = 5.8
          )
        }
        if ("scenario10" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario10", 10L,
            MUTATION_PANEL_GENES, "Mutation", cancer,
            "genome", genome_modal, "MsigDB",
            paste_compact(MUTATION_PANEL_GENES), paste0(cancer, "_mutation_panel_vs_", genome_modal, "_genome"),
            paste0("Mutation panel-driven ", genome_modal, " genome scan in CPTAC-", cancer),
            width = 9, height = 6
          )
        }
        if ("scenario12" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario12", 12L,
            gene, "Protein", cancer,
            "genome", genome_modal, "MsigDB",
            gene, paste0(cancer, "_Protein_vs_", genome_modal, "_genome"),
            paste0(gene, " Protein-related ", genome_modal, " genome scan in CPTAC-", cancer),
            width = 8.5, height = 5.8
          )
        }
        if ("scenario14" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario14", 14L,
            CONTINUOUS_PANEL_GENES, "Protein", cancer,
            "genome", genome_modal, "MsigDB",
            paste_compact(CONTINUOUS_PANEL_GENES), paste0(cancer, "_protein_panel_vs_", genome_modal, "_genome"),
            paste0("Protein panel-related ", genome_modal, " genome scan in CPTAC-", cancer),
            width = 9, height = 6
          )
        }
      }

      for (database in ENRICH_DATABASES) {
        if ("scenario09" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario09", 9L,
            PRIMARY_MUTATION_GENE, "Mutation", cancer,
            "enrichment", "Protein", database,
            PRIMARY_MUTATION_GENE, paste0(cancer, "_Mutation_GSEA_", database),
            paste0(PRIMARY_MUTATION_GENE, " Mutation-driven ", database, " GSEA in CPTAC-", cancer),
            width = 8, height = 6.2
          )
        }
        if ("scenario11" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario11", 11L,
            MUTATION_PANEL_GENES, "Mutation", cancer,
            "enrichment", "Protein", database,
            paste_compact(MUTATION_PANEL_GENES), paste0(cancer, "_mutation_panel_GSEA_", database),
            paste0("Mutation panel-driven ", database, " GSEA in CPTAC-", cancer),
            width = 9, height = 6.5
          )
        }
        if ("scenario13" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario13", 13L,
            gene, "Protein", cancer,
            "enrichment", "Protein", database,
            gene, paste0(cancer, "_Protein_GSEA_", database),
            paste0(gene, " Protein-related ", database, " GSEA in CPTAC-", cancer),
            width = 8, height = 6.2
          )
        }
        if ("scenario15" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario15", 15L,
            CONTINUOUS_PANEL_GENES, "Protein", cancer,
            "enrichment", "Protein", database,
            paste_compact(CONTINUOUS_PANEL_GENES), paste0(cancer, "_protein_panel_GSEA_", database),
            paste0("Protein panel-related ", database, " GSEA in CPTAC-", cancer),
            width = 9, height = 6.5
          )
        }
      }

      for (surv_type in SURVIVAL_TYPES) {
        if ("scenario16" %in% selected_analyses) {
          tasks <- add_survival_task(
            tasks, "scenario16", 16L,
            gene, "RNAseq", cancer,
            surv_type, gene, paste0(cancer, "_RNAseq_", surv_type),
            paste0(gene, " RNAseq ", surv_type, " survival in CPTAC-", cancer),
            width = 8.5, height = 4.8
          )
          tasks <- add_survival_task(
            tasks, "scenario16", 16L,
            gene, "Protein", cancer,
            surv_type, gene, paste0(cancer, "_Protein_", surv_type),
            paste0(gene, " Protein ", surv_type, " survival in CPTAC-", cancer),
            width = 8.5, height = 4.8
          )
        }
        if ("scenario17" %in% selected_analyses) {
          tasks <- add_survival_task(
            tasks, "scenario17", 17L,
            CONTINUOUS_PANEL_GENES, "RNAseq", cancer,
            surv_type, paste_compact(CONTINUOUS_PANEL_GENES), paste0(cancer, "_RNAseq_panel_", surv_type),
            paste0("RNAseq panel ", surv_type, " forest plot in CPTAC-", cancer),
            width = 8, height = 6
          )
          tasks <- add_survival_task(
            tasks, "scenario17", 17L,
            gene, "RNAseq", pan_common,
            surv_type, gene, paste0("pan_cancer_RNAseq_", surv_type),
            paste0(gene, " RNAseq ", surv_type, " pan-cancer forest plot"),
            width = 8, height = 6
          )
        }
      }
    }
  }

  if (RUN_AUXILIARY_AVAILABLE_SCENARIOS && length(AUXILIARY_CANCERS) > 0L) {
    aux_phospho_cancers <- available_cancers_for_modal("Phospho", AUXILIARY_CANCERS)
    aux_methylation_cancers <- available_cancers_for_modal("Methylation", AUXILIARY_CANCERS)

    if (length(aux_phospho_cancers) > 0L) {
      cancer <- aux_phospho_cancers[1]
      if ("scenario02" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario02", 2L,
          AUXILIARY_PHOSPHO_GENE, "Protein", cancer,
          AUXILIARY_PHOSPHO_GENE, "Phospho", cancer,
          AUXILIARY_PHOSPHO_GENE, paste0(cancer, "_aux_Protein_vs_Phospho"),
          paste0(AUXILIARY_PHOSPHO_GENE, " Protein vs Phospho in CPTAC-", cancer),
          width = 7, height = 4.8
        )
      }
      if ("scenario03" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario03", 3L,
          AUXILIARY_PHOSPHO_GENE, "Phospho", cancer,
          AUXILIARY_PHOSPHO_GENE, "Phospho", cancer,
          AUXILIARY_PHOSPHO_GENE, paste0(cancer, "_aux_Phospho_vs_Phospho"),
          paste0(AUXILIARY_PHOSPHO_GENE, " phospho-site correlation in CPTAC-", cancer),
          width = 7, height = 5.6
        )
      }
    }

    if (length(aux_methylation_cancers) > 0L) {
      cancer <- aux_methylation_cancers[1]
      if ("scenario01" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario01", 1L,
          main_gene, "Methylation", cancer,
          main_gene, "RNAseq", cancer,
          main_gene, paste0(cancer, "_aux_Methylation_vs_RNAseq"),
          paste0(main_gene, " Methylation vs RNAseq in CPTAC-", cancer),
          width = 5.4, height = 4.8
        )
      }
    }
  }

  tasks
}

assign_runtime_workers_to_tasks <- function(tasks, nested_workers) {
  nested_workers <- max(as.integer(nested_workers), 1L)
  if (!is.na(ENRICHMENT_WORKERS)) {
    nested_workers <- max(as.integer(ENRICHMENT_WORKERS), 1L)
  }

  for (i in seq_along(tasks)) {
    if (identical(tasks[[i]]$function_name, "cptac_enrichment")) {
      tasks[[i]]$args$n_workers <- nested_workers
    }
  }
  tasks
}

tasks_to_data_frame <- function(tasks) {
  if (length(tasks) == 0L) {
    return(data.frame())
  }

  do.call(rbind, lapply(tasks, function(task) {
    requirements <- task$requirements
    data.frame(
      Task_ID = task$task_id %||% NA_integer_,
      Analysis = task$analysis,
      Scenario_ID = task$scenario_id,
      Function = task$function_name,
      Target = task$target,
      Context = task$context,
      Title = make_default_task_title(task),
      Arguments = paste(names(task$args), collapse = ";"),
      Requirement_Modals = paste(requirements$Modal, collapse = ";"),
      Requirement_Cancers = paste(requirements$Cancers, collapse = " | "),
      Width = task$width,
      Height = task$height,
      stringsAsFactors = FALSE
    )
  }))
}

validate_slcptac_tasks <- function(tasks) {
  validation <- if (length(tasks) == 0L) {
    data.frame()
  } else {
    do.call(rbind, lapply(seq_along(tasks), function(i) {
      task <- tasks[[i]]
      fn <- get(task$function_name, envir = asNamespace("SLCPTAC"), mode = "function")
      formal_names <- names(formals(fn))
      arg_names <- names(task$args)
      unknown_args <- setdiff(arg_names, formal_names)

      data.frame(
        Task_ID = i,
        Analysis = task$analysis,
        Function = task$function_name,
        Target = task$target,
        Context = task$context,
        Argument_Names = paste(arg_names, collapse = ";"),
        Unknown_Arguments = paste(unknown_args, collapse = ";"),
        Status = ifelse(length(unknown_args) == 0L, "ok", "invalid_args"),
        stringsAsFactors = FALSE
      )
    }))
  }

  write_table(validation, TABLE_ROOT, "000_slcptac_task_argument_validation")
  invalid <- validation[validation$Status != "ok", , drop = FALSE]
  if (nrow(invalid) > 0L) {
    stop(
      "SLCPTAC task argument validation failed. See: ",
      file.path(TABLE_ROOT, "000_slcptac_task_argument_validation.csv")
    )
  }
  invisible(validation)
}


# 6. 数据审计 ------------------------------------------------------------------

bulk_gene_file <- function(gene) {
  file.path(SLCPTAC_BULK_DATA_ROOT, "CPTAC_Omics_Split", paste0(gene, "_cptac.qs"))
}

bulk_genome_file <- function(modal) {
  file_mapping <- c(
    RNAseq = "LinkedOmicsKB_PanCancer_RNAseq_RSEM.qs",
    Protein = "LinkedOmicsKB_PanCancer_Protein_Quantification.qs",
    Phospho = "LinkedOmicsKB_PanCancer_Phospho_Quantification.qs",
    Methylation = "LinkedOmicsKB_PanCancer_Methylation.qs",
    logCNA = "LinkedOmicsKB_PanCancer_CNV_logCNA.qs",
    Mutation = "LinkedOmicsKB_PanCancer_Mutation_Binary.qs",
    Clinical = "LinkedOmicsKB_PanCancer_Clin.qs",
    Survival = "LinkedOmicsKB_PanCancer_Clin.qs"
  )
  file.path(SLCPTAC_BULK_DATA_ROOT, unname(file_mapping[[modal]]))
}

extract_task_genes <- function(task) {
  args <- task$args
  genes <- character(0)
  if (!is.null(args$var1) && !identical(args$var1_modal, "Clinical")) {
    genes <- c(genes, args$var1)
  }
  if (!is.null(args$var2) && !identical(args$var2_modal, "Clinical")) {
    genes <- c(genes, args$var2)
  }
  genes <- unique(trimws(as.character(genes)))
  genes[nzchar(genes)]
}

make_required_file_table <- function(tasks) {
  required_files <- list()
  add_file <- function(task_id, task_name, role, file) {
    required_files[[length(required_files) + 1L]] <<- data.frame(
      Task_ID = task_id,
      Analysis = task_name,
      Role = role,
      File = normalizePath(file, winslash = "/", mustWork = FALSE),
      Exists = file.exists(file),
      stringsAsFactors = FALSE
    )
  }

  add_file(NA_integer_, "global", "bulk_root", SLCPTAC_BULK_DATA_ROOT)
  add_file(NA_integer_, "global", "clinical_reference_gene", bulk_gene_file("TP53"))
  add_file(NA_integer_, "global", "clinical_bulk", bulk_genome_file("Clinical"))

  for (task in tasks) {
    for (gene in extract_task_genes(task)) {
      add_file(task$task_id, task$analysis, paste0("gene_", gene), bulk_gene_file(gene))
    }
    for (modal in unique(task$requirements$Modal)) {
      if (modal %in% c("Phospho", "Methylation", "RNAseq", "Protein", "logCNA", "Mutation", "Clinical", "Survival")) {
        add_file(task$task_id, task$analysis, paste0("bulk_", modal), bulk_genome_file(modal))
      }
    }
  }

  table <- do.call(rbind, required_files)
  table <- table[!duplicated(table[, c("Role", "File")]), , drop = FALSE]
  rownames(table) <- NULL
  table
}

run_data_summary <- function(tasks) {
  dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
  dir.create(DATA_ROOT, recursive = TRUE, showWarnings = FALSE)
  dir.create(SLCPTAC_REFERENCE_CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)

  package_info <- data.frame(
    Package = "SLCPTAC",
    Version = as.character(utils::packageVersion("SLCPTAC")),
    Library_Path = find.package("SLCPTAC"),
    Official_GitHub = "https://github.com/SolvingLab/SLCPTAC",
    SL_BULK_DATA = Sys.getenv("SL_BULK_DATA"),
    Data_Root = DATA_ROOT,
    Temporary_Root = TEMP_ROOT,
    Result_Root = RESULT_ROOT,
    stringsAsFactors = FALSE
  )
  write_table(package_info, TABLE_ROOT, "000_slcptac_package_info")

  availability <- make_modal_availability_table()
  write_table(availability, TABLE_ROOT, "000_slcptac_modal_availability")

  task_design <- tasks_to_data_frame(tasks)
  write_table(task_design, TABLE_ROOT, "000_slcptac_task_design")

  required_files <- make_required_file_table(tasks)
  write_table(required_files, TABLE_ROOT, "000_slcptac_required_files")

  config <- data.frame(
    Parameter = c(
      "TARGET_GENES", "TARGET_CANCERS", "PAN_CANCERS", "CONTINUOUS_PANEL_GENES",
      "MUTATION_PANEL_GENES", "CLINICAL_VARIABLES", "ENRICH_DATABASES",
      "GENOME_SCAN_MODALS", "SURVIVAL_TYPES", "RUN_MODAL_PAIR_GRID",
      "RUN_AUXILIARY_AVAILABLE_SCENARIOS"
    ),
    Value = c(
      paste(TARGET_GENES, collapse = ";"),
      paste(TARGET_CANCERS, collapse = ";"),
      paste(PAN_CANCERS, collapse = ";"),
      paste(CONTINUOUS_PANEL_GENES, collapse = ";"),
      paste(MUTATION_PANEL_GENES, collapse = ";"),
      paste(CLINICAL_VARIABLES, collapse = ";"),
      paste(ENRICH_DATABASES, collapse = ";"),
      paste(GENOME_SCAN_MODALS, collapse = ";"),
      paste(SURVIVAL_TYPES, collapse = ";"),
      as.character(RUN_MODAL_PAIR_GRID),
      as.character(RUN_AUXILIARY_AVAILABLE_SCENARIOS)
    ),
    stringsAsFactors = FALSE
  )
  write_table(config, TABLE_ROOT, "000_slcptac_run_configuration")

  missing_required <- required_files[!required_files$Exists, , drop = FALSE]
  write_table(missing_required, TABLE_ROOT, "000_slcptac_missing_required_files")

  if (STOP_IF_BULK_DATA_MISSING && nrow(missing_required) > 0L) {
    stop(
      "SLCPTAC bulk data files are missing. See: ",
      file.path(TABLE_ROOT, "000_slcptac_missing_required_files.csv")
    )
  }

  package_info
}


# 7. 绘图保存与文件复制 --------------------------------------------------------

get_pdf_page_count <- function(pdf_file) {
  if (!file.exists(pdf_file) || is.na(file.info(pdf_file)$size) || file.info(pdf_file)$size <= 0) {
    return(NA_integer_)
  }

  if (requireNamespace("pdftools", quietly = TRUE)) {
    page_count <- tryCatch(pdftools::pdf_info(pdf_file)$pages, error = function(error) NA_integer_)
    if (!is.na(page_count)) {
      return(as.integer(page_count))
    }
  }

  pdfinfo <- Sys.which("pdfinfo")
  if (!nzchar(pdfinfo)) {
    return(1L)
  }

  info <- tryCatch(
    suppressWarnings(system2(pdfinfo, args = pdf_file, stdout = TRUE, stderr = TRUE)),
    error = function(error) character(0)
  )
  page_line <- grep("^Pages:", info, value = TRUE)
  if (length(page_line) == 0L) {
    return(1L)
  }

  suppressWarnings(as.integer(trimws(sub("^Pages:", "", page_line[1]))))
}

get_png_nonwhite_fraction <- function(png_file, white_tolerance = 0.02) {
  if (!file.exists(png_file) || !requireNamespace("png", quietly = TRUE)) {
    return(NA_real_)
  }

  img <- tryCatch(png::readPNG(png_file), error = function(error) NULL)
  if (is.null(img)) {
    return(NA_real_)
  }

  dims <- dim(img)
  if (length(dims) == 2L) {
    return(mean(abs(img - 1) > white_tolerance))
  }

  channels <- min(3L, dims[3])
  nonwhite <- matrix(FALSE, nrow = dims[1], ncol = dims[2])
  for (channel in seq_len(channels)) {
    nonwhite <- nonwhite | abs(img[, , channel] - 1) > white_tolerance
  }
  if (dims[3] >= 4L) {
    nonwhite <- nonwhite & img[, , 4] > 0.01
  }
  mean(nonwhite)
}

is_png_blank <- function(png_file, nonwhite_fraction_cutoff = 0.001) {
  nonwhite_fraction <- get_png_nonwhite_fraction(png_file)
  is.na(nonwhite_fraction) || nonwhite_fraction < nonwhite_fraction_cutoff
}

render_one_pdf_page_png <- function(pdf_file, page, png_file, dpi = PNG_DPI) {
  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)
  unlink(png_file, force = TRUE)

  if (requireNamespace("pdftools", quietly = TRUE) && requireNamespace("png", quietly = TRUE)) {
    converted <- tryCatch(
      {
        bitmap <- pdftools::pdf_render_page(pdf_file, page = page, dpi = dpi)
        png::writePNG(bitmap, target = png_file)
        file.exists(png_file)
      },
      error = function(error) FALSE
    )
    if (converted) {
      return(png_file)
    }
  }

  pdftoppm <- Sys.which("pdftoppm")
  if (nzchar(pdftoppm)) {
    tmp_dir <- tempfile("slcptac_pdftoppm_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
    tmp_prefix <- file.path(tmp_dir, "page")
    status <- tryCatch(
      suppressWarnings(system2(
        pdftoppm,
        args = c("-png", "-r", as.character(dpi), "-f", as.character(page), "-l", as.character(page), pdf_file, tmp_prefix),
        stdout = TRUE,
        stderr = TRUE
      )),
      error = function(error) structure(character(0), status = 1L)
    )
    rendered <- sort(list.files(tmp_dir, pattern = "^page-.*[.]png$", full.names = TRUE))
    if (length(rendered) > 0L && is.null(attr(status, "status"))) {
      file.copy(rendered[1], png_file, overwrite = TRUE)
      return(png_file)
    }
  }

  sips <- Sys.which("sips")
  if (nzchar(sips) && identical(as.integer(page), 1L)) {
    status <- tryCatch(
      suppressWarnings(system2(
        sips,
        args = c("-s", "format", "png", pdf_file, "--out", png_file),
        stdout = TRUE,
        stderr = TRUE
      )),
      error = function(error) structure(character(0), status = 1L)
    )
    if (file.exists(png_file) && is.null(attr(status, "status"))) {
      return(png_file)
    }
  }

  character(0)
}

convert_pdf_to_png_outputs <- function(pdf_file, png_file, dpi = PNG_DPI) {
  page_count <- get_pdf_page_count(pdf_file)
  if (is.na(page_count) || page_count < 1L) {
    return(character(0))
  }

  output_stem <- tools::file_path_sans_ext(basename(png_file))
  if (page_count == 1L) {
    rendered <- render_one_pdf_page_png(pdf_file, page = 1L, png_file = png_file, dpi = dpi)
    if (length(rendered) == 0L || is_png_blank(png_file)) {
      return(character(0))
    }
    return(png_file)
  }

  output_dir <- file.path(dirname(png_file), output_stem)
  unlink(output_dir, recursive = TRUE, force = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  png_files <- character(0)
  for (page in seq_len(page_count)) {
    page_png <- file.path(output_dir, sprintf("%s_page_%03d.png", output_stem, page))
    rendered <- render_one_pdf_page_png(pdf_file, page = page, png_file = page_png, dpi = dpi)
    if (length(rendered) > 0L && !is_png_blank(page_png)) {
      png_files <- c(png_files, page_png)
    }
  }

  png_files
}

plot_has_title <- function(plot) {
  if (!inherits(plot, "ggplot")) {
    return(FALSE)
  }
  title <- plot$labels$title
  !is.null(title) && nzchar(as.character(title))
}

add_title_to_plot_if_needed <- function(plot, title) {
  if (!ADD_SCRIPT_TITLES || is.null(title) || !nzchar(title)) {
    return(plot)
  }

  if (inherits(plot, "patchwork") && requireNamespace("patchwork", quietly = TRUE)) {
    return(plot + patchwork::plot_annotation(title = title))
  }

  if (inherits(plot, "ggplot") && !plot_has_title(plot)) {
    return(plot + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }

  plot
}

draw_outer_title <- function(title) {
  if (!ADD_SCRIPT_TITLES || is.null(title) || !nzchar(title)) {
    return(invisible(FALSE))
  }
  grid::grid.text(
    label = title,
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.985, "npc"),
    gp = grid::gpar(fontface = "bold", fontsize = 13)
  )
  invisible(TRUE)
}

draw_in_title_viewport <- function(title, draw_fun) {
  if (!ADD_SCRIPT_TITLES || is.null(title) || !nzchar(title)) {
    draw_fun()
    return(invisible(NULL))
  }

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.0, "npc"),
    width = grid::unit(1, "npc"),
    height = grid::unit(0.92, "npc"),
    just = c("center", "bottom")
  ))
  on.exit({
    try(grid::popViewport(), silent = TRUE)
  }, add = TRUE)
  draw_fun()
  grid::popViewport()
  draw_outer_title(title)
  invisible(NULL)
}

render_slcptac_plot <- function(plot, title = NULL) {
  if (is.null(plot)) {
    draw_outer_title(title)
    return(invisible(NULL))
  }

  if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) {
    print(add_title_to_plot_if_needed(plot, title))
    return(invisible(NULL))
  }

  if (inherits(plot, "Heatmap") && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    draw_in_title_viewport(title, function() {
      ComplexHeatmap::draw(plot, newpage = FALSE)
    })
    return(invisible(NULL))
  }

  if (inherits(plot, "grob") || inherits(plot, "gTree") || inherits(plot, "gtable")) {
    draw_in_title_viewport(title, function() {
      grid::grid.draw(plot)
    })
    return(invisible(NULL))
  }

  if (is.list(plot) && "plot" %in% names(plot)) {
    render_slcptac_plot(plot$plot, title = title)
    return(invisible(NULL))
  }

  try(print(plot), silent = TRUE)
  draw_outer_title(title)
  invisible(NULL)
}

get_plot_size <- function(plot, task) {
  width <- suppressWarnings(as.numeric(attr(plot, "width")))
  height <- suppressWarnings(as.numeric(attr(plot, "height")))
  if (is.na(width) || width <= 0) {
    width <- task$width
  }
  if (is.na(height) || height <= 0) {
    height <- task$height
  }

  title_extra <- estimate_title_extra_height(make_default_task_title(task))
  list(
    width = max(width, task$width),
    height = max(height + title_extra, task$height)
  )
}

save_slcptac_pdf_png <- function(task, plot) {
  plot_size <- get_plot_size(plot, task)
  task_stem <- task_file_stem(task)
  pdf_file <- file.path(PLOT_PDF_DIR, paste0(task_stem, ".pdf"))
  png_file <- file.path(PLOT_PNG_DIR, paste0(task_stem, ".png"))

  dir.create(PLOT_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(PLOT_PNG_DIR, recursive = TRUE, showWarnings = FALSE)
  unlink(c(pdf_file, png_file), force = TRUE)
  unlink(file.path(PLOT_PNG_DIR, task_stem), recursive = TRUE, force = TRUE)

  device_opened <- FALSE
  Cairo::CairoPDF(
    file = pdf_file,
    width = plot_size$width,
    height = plot_size$height,
    bg = "white"
  )
  device_opened <- TRUE
  on.exit({
    if (device_opened && grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  }, add = TRUE)

  render_slcptac_plot(plot, title = make_default_task_title(task))
  invisible(grDevices::dev.off())
  device_opened <- FALSE

  if (!file.exists(pdf_file) || is.na(file.info(pdf_file)$size) || file.info(pdf_file)$size <= 0) {
    stop("PDF output is missing or empty: ", pdf_file)
  }

  png_files <- convert_pdf_to_png_outputs(pdf_file, png_file, dpi = PNG_DPI)
  list(pdf_file = pdf_file, png_file = png_files)
}

copy_generated_files <- function(task) {
  generated_dir <- file.path(getwd(), "slcptac_output")
  if (!dir.exists(generated_dir)) {
    return(character(0))
  }

  generated_files <- list.files(generated_dir, all.files = FALSE, recursive = FALSE, full.names = TRUE)
  generated_files <- generated_files[file.info(generated_files)$isdir == FALSE]
  if (length(generated_files) == 0L) {
    return(character(0))
  }

  copied_files <- character(0)
  for (source_file in generated_files) {
    extension <- tolower(tools::file_ext(source_file))
    source_stem <- tools::file_path_sans_ext(basename(source_file))
    output_stem <- task_file_stem(task, paste0("package_", source_stem))

    destination_dir <- if (extension %in% c("png", "jpg", "jpeg")) {
      PLOT_PNG_DIR
    } else if (extension == "pdf") {
      PLOT_PDF_DIR
    } else {
      TABLE_ROOT
    }
    dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
    destination_file <- file.path(destination_dir, paste0(output_stem, ".", extension))
    file.copy(source_file, destination_file, overwrite = TRUE)
    copied_files <- c(copied_files, destination_file)
  }

  copied_files
}

validate_plot_output_files <- function(files) {
  files <- unique(files[nzchar(files)])
  if (length(files) == 0L) {
    return(invisible(TRUE))
  }

  invalid <- character(0)
  for (file in files) {
    if (!file.exists(file) || is.na(file.info(file)$size) || file.info(file)$size <= 0) {
      invalid <- c(invalid, paste0(file, " [missing_or_empty]"))
      next
    }
    if (tolower(tools::file_ext(file)) == "png" && is_png_blank(file)) {
      invalid <- c(invalid, paste0(file, " [blank_png]"))
    }
  }

  if (length(invalid) > 0L) {
    stop("Invalid plot output detected: ", paste(invalid, collapse = "; "))
  }
  invisible(TRUE)
}

write_result_tables <- function(task, result) {
  output_files <- character(0)
  task_stem <- task_file_stem(task)

  if (is.list(result) && "stats" %in% names(result)) {
    stats_file <- file.path(TABLE_ROOT, paste0(task_stem, "_stats.csv"))
    write_csv_with_report_previews(normalize_result_table(result$stats), stats_file, na = "NA")
    output_files <- c(output_files, stats_file)
  }

  if (SAVE_RAW_DATA_TABLES && is.list(result) && "raw_data" %in% names(result)) {
    raw_data <- result$raw_data
    if (is.data.frame(raw_data) || is.matrix(raw_data)) {
      raw_file <- file.path(TABLE_ROOT, paste0(task_stem, "_raw_data.csv"))
      write_csv_with_report_previews(normalize_result_table(raw_data), raw_file, na = "NA")
      output_files <- c(output_files, raw_file)
    } else if (is.list(raw_data)) {
      summary_file <- file.path(TABLE_ROOT, paste0(task_stem, "_raw_data_list_summary.csv"))
      write_csv_with_report_previews(list_to_data_frame(raw_data), summary_file, na = "NA")
      output_files <- c(output_files, summary_file)

      for (i in seq_along(raw_data)) {
        item <- raw_data[[i]]
        if (is.data.frame(item) || is.matrix(item)) {
          item_name <- names(raw_data)[i] %||% paste0("item_", i)
          item_file <- file.path(
            TABLE_ROOT,
            paste0(task_stem, "_raw_data_", sprintf("%03d", i), "_", safe_name(item_name), ".csv")
          )
          write_csv_with_report_previews(normalize_result_table(item), item_file, na = "NA")
          output_files <- c(output_files, item_file)
        }
      }
    }
  }

  if (SAVE_RESULT_OBJECTS) {
    rds_file <- file.path(TABLE_ROOT, paste0(task_stem, "_result_object.rds"))
    saveRDS(result, rds_file)
    output_files <- c(output_files, rds_file)
  }

  output_files
}


# 8. 任务缓存与执行器 ----------------------------------------------------------

CACHEABLE_SLCPTAC_ANALYSES <- c(
  "scenario08", "scenario09", "scenario10", "scenario11",
  "scenario12", "scenario13", "scenario14", "scenario15"
)

split_manifest_files <- function(x) {
  x <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

manifest_files_exist <- function(manifest) {
  if (!is.data.frame(manifest) || nrow(manifest) == 0L) {
    return(FALSE)
  }
  files <- unique(c(
    split_manifest_files(manifest$Output_Files[1]),
    split_manifest_files(manifest$Generated_Files[1])
  ))
  length(files) > 0L && all(file.exists(files))
}

is_cacheable_slcptac_task <- function(task) {
  USE_SLCPTAC_TASK_CACHE && task$analysis %in% CACHEABLE_SLCPTAC_ANALYSES
}

get_slcptac_task_cache_file <- function(task) {
  file.path(
    SLCPTAC_TASK_CACHE_ROOT,
    safe_name(task$analysis),
    safe_name(task$target),
    paste0(safe_name(task$context), ".rds")
  )
}

make_slcptac_task_cache_metadata <- function(task) {
  make_resource_cache_metadata(
    resource_name = paste0("SLCPTAC::", task$function_name),
    species = "human",
    extra = list(
      package_version = as.character(utils::packageVersion("SLCPTAC")),
      output_layout = "flat_pdf_png_tables_v1",
      task_id = task$task_id %||% NA_integer_,
      analysis = task$analysis,
      scenario_id = task$scenario_id,
      target = task$target,
      context = task$context,
      width = task$width,
      height = task$height,
      args = task$args
    )
  )
}

read_slcptac_task_cache <- function(task) {
  if (!is_cacheable_slcptac_task(task)) {
    return(NULL)
  }

  cache_file <- get_slcptac_task_cache_file(task)
  metadata <- make_slcptac_task_cache_metadata(task)
  cached <- read_reference_cache(
    cache_file = cache_file,
    expected_metadata = metadata,
    max_age_days = SLCPTAC_TASK_CACHE_MAX_AGE_DAYS,
    use_cache = USE_SLCPTAC_TASK_CACHE
  )

  if (!cached$found || !manifest_files_exist(cached$result)) {
    return(NULL)
  }

  cached$result$Status <- "cached"
  cached$result$Cache_File <- cache_file
  cached$result
}

write_slcptac_task_cache <- function(task, manifest) {
  if (!is_cacheable_slcptac_task(task) || !identical(manifest$Status[1], "success")) {
    return(invisible(FALSE))
  }

  cache_file <- get_slcptac_task_cache_file(task)
  metadata <- make_slcptac_task_cache_metadata(task)
  write_reference_cache(cache_file = cache_file, metadata = metadata, result = manifest)
  invisible(TRUE)
}

get_slcptac_task_skip_reason <- function(task) {
  for (i in seq_len(nrow(task$requirements))) {
    requirement <- task$requirements[i, , drop = FALSE]
    cancers <- unlist(strsplit(requirement$Cancers, ";", fixed = TRUE), use.names = FALSE)
    cancers <- cancers[nzchar(cancers)]
    if (!modal_available_for_cancers(requirement$Modal, cancers)) {
      unavailable <- setdiff(cancers, MODAL_AVAILABILITY[[requirement$Modal]] %||% character(0))
      return(paste0(
        "Skipped because ", requirement$Modal, " data is unavailable for ",
        paste(unavailable, collapse = ", "),
        ". This follows SLCPTAC official data coverage."
      ))
    }
  }

  if (!dir.exists(SLCPTAC_BULK_DATA_ROOT)) {
    return(paste0(
      "Skipped because SL_BULK_DATA does not exist: ",
      normalizePath(SLCPTAC_BULK_DATA_ROOT, winslash = "/", mustWork = FALSE)
    ))
  }

  ""
}

with_task_workspace <- function(task, expr) {
  task_temp_dir <- file.path(
    TEMP_ROOT,
    "tasks",
    safe_name(task$analysis),
    safe_name(task$target),
    safe_name(task$context)
  )
  dir.create(task_temp_dir, recursive = TRUE, showWarnings = FALSE)
  if (CLEAN_TASK_OUTPUT_DIR) {
    unlink(list.files(task_temp_dir, all.files = FALSE, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)
  }

  old_wd <- getwd()
  setwd(task_temp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  force(expr)
}

make_task_result <- function(
    task,
    status,
    output_files = character(0),
    generated_files = character(0),
    warnings = character(0),
    messages = character(0),
    error = "",
    start_time = NULL,
    end_time = Sys.time(),
    cache_file = "") {
  runtime_seconds <- if (is.null(start_time)) {
    NA_real_
  } else {
    as.numeric(difftime(end_time, start_time, units = "secs"))
  }

  data.frame(
    Analysis = task$analysis,
    Scenario_ID = task$scenario_id,
    Function = task$function_name,
    Target = task$target,
    Context = task$context,
    Status = status,
    Start_Time = format_task_timestamp(start_time),
    End_Time = format_task_timestamp(end_time),
    Runtime_Seconds = runtime_seconds,
    Runtime = ifelse(is.na(runtime_seconds), "", format_runtime_seconds(runtime_seconds)),
    Output_Files = paste(output_files, collapse = ";"),
    Generated_Files = paste(generated_files, collapse = ";"),
    Cache_File = cache_file,
    Warning = paste(warnings, collapse = " | "),
    Message = paste(messages, collapse = " | "),
    Error = error,
    stringsAsFactors = FALSE
  )
}

run_one_slcptac_task <- function(task) {
  task_start_time <- Sys.time()
  cache_file <- if (is_cacheable_slcptac_task(task)) {
    get_slcptac_task_cache_file(task)
  } else {
    ""
  }

  skip_reason <- get_slcptac_task_skip_reason(task)
  if (nzchar(skip_reason)) {
    return(make_task_result(
      task = task,
      status = "skipped",
      messages = skip_reason,
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  cached_manifest <- read_slcptac_task_cache(task)
  if (!is.null(cached_manifest)) {
    return(make_task_result(
      task = task,
      status = "cached",
      output_files = split_manifest_files(cached_manifest$Output_Files[1]),
      generated_files = split_manifest_files(cached_manifest$Generated_Files[1]),
      messages = "Reused SLCPTAC task cache because output files already exist",
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  captured <- with_task_workspace(task, capture_task({
    fn <- get(task$function_name, envir = asNamespace("SLCPTAC"), mode = "function")
    result <- do.call(fn, task$args)

    generated_files <- copy_generated_files(task)
    table_files <- write_result_tables(task, result)

    plot_files <- character(0)
    if (is.list(result) && "plot" %in% names(result)) {
      saved_plot <- save_slcptac_pdf_png(task, result$plot)
      plot_files <- unname(unlist(saved_plot))
    }

    list(
      output_files = c(plot_files, table_files),
      generated_files = generated_files
    )
  }))

  if (inherits(captured$result, "error")) {
    return(make_task_result(
      task = task,
      status = "failed",
      warnings = captured$warnings,
      messages = captured$messages,
      error = conditionMessage(captured$result),
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  validation_error <- tryCatch(
    {
      validate_plot_output_files(c(captured$result$output_files, captured$result$generated_files))
      NULL
    },
    error = function(error) error
  )
  if (inherits(validation_error, "error")) {
    return(make_task_result(
      task = task,
      status = "failed",
      warnings = captured$warnings,
      messages = captured$messages,
      error = conditionMessage(validation_error),
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  manifest <- make_task_result(
    task = task,
    status = "success",
    output_files = captured$result$output_files,
    generated_files = captured$result$generated_files,
    warnings = captured$warnings,
    messages = captured$messages,
    start_time = task_start_time,
    cache_file = cache_file
  )
  write_slcptac_task_cache(task, manifest)
  manifest
}

normalize_parallel_slcptac_results <- function(results, tasks) {
  normalized <- vector("list", length(results))
  for (i in seq_along(results)) {
    result <- results[[i]]
    if (inherits(result, "try-error") || is.null(result)) {
      normalized[[i]] <- make_task_result(
        task = tasks[[i]],
        status = "failed",
        error = paste(as.character(result), collapse = "\n")
      )
    } else {
      normalized[[i]] <- result
    }
  }
  normalized
}

is_success_status <- function(status) {
  status %in% c("success", "cached", "skipped")
}

make_runtime_summary <- function(task_summary, start_time) {
  total_runtime_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  task_status <- if (nrow(task_summary) > 0L) task_summary$Status else character(0)

  data.frame(
    Section = c("slcptac_tasks", "total"),
    Tasks = c(nrow(task_summary), nrow(task_summary)),
    Completed_Without_Error = c(sum(is_success_status(task_status)), sum(is_success_status(task_status))),
    Failed = c(sum(!is_success_status(task_status)), sum(!is_success_status(task_status))),
    Runtime_Seconds = c(
      sum(suppressWarnings(as.numeric(task_summary$Runtime_Seconds)), na.rm = TRUE),
      total_runtime_seconds
    ),
    Runtime = c(
      format_runtime_seconds(sum(suppressWarnings(as.numeric(task_summary$Runtime_Seconds)), na.rm = TRUE)),
      format_runtime_seconds(total_runtime_seconds)
    ),
    stringsAsFactors = FALSE
  )
}

clear_previous_slcptac_run_outputs <- function() {
  if (!CLEAR_PREVIOUS_RUN_OUTPUTS) {
    return(invisible(FALSE))
  }

  unlink(PLOT_ROOT, recursive = TRUE, force = TRUE)
  unlink(TABLE_ROOT, recursive = TRUE, force = TRUE)
  unlink(TEMP_ROOT, recursive = TRUE, force = TRUE)
  unlink(SLCPTAC_TASK_CACHE_ROOT, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}


# 9. 主运行入口 ----------------------------------------------------------------

clear_previous_slcptac_run_outputs()

dir.create(RESULT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_PNG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(SLCPTAC_REFERENCE_CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(SLCPTAC_TASK_CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)

selected_analyses <- resolve_requested_analyses()
slcptac_tasks <- build_slcptac_tasks(selected_analyses)
slcptac_tasks <- assign_task_indices(slcptac_tasks)

parallel_strategy <- setup_parallel_strategy(
  total_tasks = max(length(slcptac_tasks), 1L),
  max_workers = PARALLEL_WORKERS,
  inner_label = "SLCPTAC inner workers",
  nested_label = "SLCPTAC enrichment workers"
)
slcptac_tasks <- assign_runtime_workers_to_tasks(
  tasks = slcptac_tasks,
  nested_workers = parallel_strategy$nested_workers
)
validate_slcptac_tasks(slcptac_tasks)

cat("\nSLCPTAC quick analysis configuration:\n")
cat("SLCPTAC version: ", as.character(utils::packageVersion("SLCPTAC")), "\n", sep = "")
cat("Target genes: ", paste(TARGET_GENES, collapse = ", "), "\n", sep = "")
cat("Target cancers: ", paste(TARGET_CANCERS, collapse = ", "), "\n", sep = "")
cat("Pan-cancer cancers: ", paste(PAN_CANCERS, collapse = ", "), "\n", sep = "")
cat("Selected analyses: ", paste(selected_analyses, collapse = ", "), "\n", sep = "")
cat("Prepared tasks: ", length(slcptac_tasks), "\n", sep = "")
cat("Result root: ", RESULT_ROOT, "\n", sep = "")
cat("Temporary root: ", TEMP_ROOT, "\n", sep = "")
cat("SLCPTAC bulk data root: ", Sys.getenv("SL_BULK_DATA"), "\n", sep = "")
cat("SLCPTAC data/cache root: ", DATA_ROOT, "\n", sep = "")

cat("\nWriting SLCPTAC analysis catalog and data audit...\n")
write_analysis_catalog(selected_analyses)
run_data_summary(slcptac_tasks)

cat("\nSLCPTAC plotting/statistical tasks: ", length(slcptac_tasks), "\n", sep = "")

task_summary <- if (length(slcptac_tasks) > 0L) {
  raw_task_results <- run_indexed_tasks_with_progress(
    total_tasks = length(slcptac_tasks),
    task_function = function(task_id) {
      run_one_slcptac_task(slcptac_tasks[[task_id]])
    },
    workers = min(parallel_strategy$task_workers, length(slcptac_tasks)),
    progress_label = "SLCPTAC tasks"
  )
  normalized_task_results <- normalize_parallel_slcptac_results(
    results = raw_task_results,
    tasks = slcptac_tasks
  )
  do.call(rbind, normalized_task_results)
} else {
  data.frame()
}

if (nrow(task_summary) > 0L) {
  write_table(task_summary, TABLE_ROOT, "000_slcptac_task_summary")
  failed_tasks <- task_summary[!is_success_status(task_summary$Status), , drop = FALSE]
  write_table(failed_tasks, TABLE_ROOT, "000_slcptac_failed_tasks")
}

runtime_summary <- make_runtime_summary(task_summary = task_summary, start_time = SCRIPT_START_TIME)
write_table(runtime_summary, TABLE_ROOT, "000_slcptac_runtime_summary")

cat("\nSLCPTAC quick analysis finished.\n")
if (nrow(task_summary) > 0L) {
  cat("Completed without error: ", sum(is_success_status(task_summary$Status)), "\n", sep = "")
  cat("Failed tasks: ", sum(!is_success_status(task_summary$Status)), "\n", sep = "")
  cat("Task summary: ", file.path(TABLE_ROOT, "000_slcptac_task_summary.csv"), "\n", sep = "")
  cat("Failed task summary: ", file.path(TABLE_ROOT, "000_slcptac_failed_tasks.csv"), "\n", sep = "")
}
cat("Analysis catalog: ", file.path(TABLE_ROOT, "000_slcptac_analysis_catalog.csv"), "\n", sep = "")
cat("Task design: ", file.path(TABLE_ROOT, "000_slcptac_task_design.csv"), "\n", sep = "")
cat("Required file audit: ", file.path(TABLE_ROOT, "000_slcptac_required_files.csv"), "\n", sep = "")
cat("Runtime summary: ", file.path(TABLE_ROOT, "000_slcptac_runtime_summary.csv"), "\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
