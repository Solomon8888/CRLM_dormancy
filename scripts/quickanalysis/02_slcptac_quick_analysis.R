# SLCPTAC快速CPTAC多组学全套分析脚本
#
# 设计目的：
# 1. 使用SLCPTAC包快速完成ATF3在CPTAC-COAD及泛癌中的多组学探索；
# 2. 覆盖SLCPTAC官方说明中的17个分析场景，并把每个场景登记到结果目录中的catalog表；
# 3. 正式结果统一保存到results/quickanalysis/cptac，不让SLCPTAC自动生成的slcptac_output散落在项目根目录；
# 4. CPTAC bulk数据默认放在data/cptac/bulk_data，并通过SL_BULK_DATA传给SLCPTAC；
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
# SLCPTAC_TARGET_GENES=ATF3,MYC Rscript scripts/quickanalysis/02_slcptac_quick_analysis.R
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
SCRIPT_FILE <- file.path(PROJECT_ROOT, "scripts", "quickanalysis", "02_slcptac_quick_analysis.R")
DATASET_ID <- "cptac"
RESULT_ROOT <- file.path(PROJECT_ROOT, "results", "quickanalysis", DATASET_ID)
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "SLCPTAC")
TABLE_ROOT <- file.path(RESULT_ROOT, "tables", "SLCPTAC")
PLOT_PDF_DIR <- file.path(PLOT_ROOT, "pdf")
PLOT_PNG_DIR <- file.path(PLOT_ROOT, "png")
DATA_ROOT <- file.path(PROJECT_ROOT, "data", "cptac")
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
# 默认在运行前清空本脚本前次运行产生的全部结果与中间文件。
# 这会删除results/quickanalysis/cptac、temporary/slcptac，
# 以及本脚本的任务manifest缓存；不会删除data/cptac中的SLCPTAC输入数据。
CLEAR_PREVIOUS_RUN_OUTPUTS <- parse_env_logical("SLCPTAC_CLEAR_PREVIOUS_OUTPUTS", TRUE)
CLEAN_TASK_OUTPUT_DIR <- parse_env_logical("SLCPTAC_CLEAN_OUTPUT", TRUE)
SAVE_RAW_DATA_TABLES <- parse_env_logical("SLCPTAC_SAVE_RAW_DATA_TABLES", TRUE)
SAVE_RESULT_OBJECTS <- parse_env_logical("SLCPTAC_SAVE_RESULT_OBJECTS", FALSE)
ADD_SCRIPT_TITLES <- parse_env_logical("SLCPTAC_ADD_SCRIPT_TITLES", TRUE)
SMOKE_TEST_MODE <- parse_env_logical("SLCPTAC_SMOKE_TEST", FALSE)
SMOKE_TEST_SAMPLES_PER_CANCER <- parse_env_integer("SLCPTAC_SMOKE_SAMPLES_PER_CANCER", 24L)
if (SMOKE_TEST_MODE &&
    !nzchar(Sys.getenv("SL_BULK_DATA", unset = "")) &&
    !nzchar(Sys.getenv("SLCPTAC_BULK_DATA", unset = ""))) {
  SLCPTAC_BULK_DATA_ROOT <- file.path(TEMP_ROOT, "smoke_bulk_data")
}
STOP_IF_BULK_DATA_MISSING <- parse_env_logical("SLCPTAC_STOP_IF_DATA_MISSING", TRUE)
AUTO_PREPARE_CPTAC_BULK <- parse_env_logical("SLCPTAC_AUTO_PREPARE_BULK", TRUE)
LINKEDOMICS_RAW_ROOT <- file.path(DATA_ROOT, "raw_linkedomics")
LINKEDOMICS_BULK_CONVERTER_VERSION <- "2026-06-08_sample_named_unlist_v2"
LINKEDOMICS_BULK_MANIFEST_FILE <- file.path(DATA_ROOT, "bulk_data_manifest.csv")

# 包安装控制。SLCPTAC 1.2.0 官方DESCRIPTION仍写着依赖qs，但本项目统一使用qs2；
# 因此缺少SLCPTAC、或检测到已安装的SLCPTAC仍不是qs2版时，默认自动安装一个
# 仅替换读取依赖的patched版。
# 如需禁用自动安装：SLCPTAC_AUTO_INSTALL=0 Rscript scripts/quickanalysis/02_slcptac_quick_analysis.R
slcptac_is_qs2_patched <- function() {
  package_dir <- tryCatch(
    find.package("SLCPTAC", quiet = TRUE),
    error = function(error) character(0)
  )
  if (length(package_dir) == 0L || !nzchar(package_dir[1])) {
    return(FALSE)
  }

  description_file <- file.path(package_dir[1], "DESCRIPTION")
  if (!file.exists(description_file)) {
    return(FALSE)
  }

  description_text <- paste(readLines(description_file, warn = FALSE), collapse = "\n")
  grepl("Imports:.*qs2|\\n[[:space:]]+.*qs2", description_text)
}

AUTO_INSTALL_SLCPTAC <- parse_env_logical(
  "SLCPTAC_AUTO_INSTALL",
  !requireNamespace("SLCPTAC", quietly = TRUE) || !slcptac_is_qs2_patched()
)

# 并行配置。外层任务并行时，SLCPTAC内部GSEA worker会自动限制，避免过度抢核。
MAX_PARALLEL_WORKERS <- parse_env_integer(
  "SLCPTAC_PARALLEL_WORKERS",
  parse_env_integer(
    "QUICKANALYSIS_PARALLEL_WORKERS",
    parse_env_integer("PARALLEL_RUNTIME_WORKERS", NA_integer_)
  )
)
ENRICHMENT_WORKERS <- parse_env_integer("SLCPTAC_ENRICHMENT_WORKERS", NA_integer_)
ALLOW_FORK_PARALLEL <- parse_env_logical(
  "SLCPTAC_ALLOW_FORK_PARALLEL",
  !identical(Sys.info()[["sysname"]], "Darwin")
)
QUICKANALYSIS_VERBOSE <- parse_env_logical("SLCPTAC_VERBOSE", FALSE)
DISABLE_FORK_PARALLEL <- parse_env_logical("SLCPTAC_DISABLE_FORK", interactive() || !ALLOW_FORK_PARALLEL)
PARALLEL_BACKEND <- Sys.getenv(
  "SLCPTAC_PARALLEL_BACKEND",
  unset = Sys.getenv(
    "QUICKANALYSIS_PARALLEL_BACKEND",
    unset = Sys.getenv("PARALLEL_RUNTIME_BACKEND", unset = "auto")
  )
)
USE_SLCPTAC_TASK_CACHE <- parse_env_logical("SLCPTAC_USE_TASK_CACHE", TRUE)
SLCPTAC_TASK_CACHE_MAX_AGE_DAYS <- parse_env_integer("SLCPTAC_TASK_CACHE_MAX_AGE_DAYS", 30L)

options(width = 200)
options(bitmapType = "cairo")
options(quickanalysis_verbose = QUICKANALYSIS_VERBOSE)
options(parallel_runtime_force_single_line_progress = TRUE)
options(parallel_runtime_quiet_strategy = !QUICKANALYSIS_VERBOSE)
options(parallel_runtime_disable_fork = DISABLE_FORK_PARALLEL)
options(parallel_runtime_backend = PARALLEL_BACKEND)


# 2. 加载R包和项目共用函数 ----------------------------------------------------

install_slcptac_if_requested <- function() {
  if (!AUTO_INSTALL_SLCPTAC) {
    return(invisible(FALSE))
  }

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("SLCPTAC is not installed, and package 'remotes' could not be installed.")
  }

  install_cran_if_missing <- function(packages, repos = getOption("repos")) {
    packages <- setdiff(packages, rownames(installed.packages()))
    if (length(packages) == 0) {
      return(invisible(TRUE))
    }
    install.packages(packages, repos = repos)
    invisible(TRUE)
  }

  ensure_qs2_dependency <- function() {
    if (requireNamespace("qs2", quietly = TRUE)) {
      return(invisible(TRUE))
    }

    # 本项目缓存和大对象读写统一使用qs2。先走CRAN，再走qsbase r-universe。
    install_attempts <- list(
      function() install.packages(
        "qs2",
        repos = c(CRAN = "https://cloud.r-project.org")
      ),
      function() install.packages(
        "qs2",
        repos = c(qsbase = "https://qsbase.r-universe.dev", CRAN = "https://cloud.r-project.org")
      )
    )

    for (attempt in install_attempts) {
      try(suppressWarnings(attempt()), silent = TRUE)
      if (requireNamespace("qs2", quietly = TRUE)) {
        return(invisible(TRUE))
      }
    }

    stop("Dependency package 'qs2' could not be installed.")
  }

  prepare_slcptac_qs2_source <- function() {
    source_candidates <- c(
      file.path(PROJECT_ROOT, "temporary", "slcptac_source"),
      file.path(TEMP_ROOT, "slcptac_source")
    )
    existing_source <- source_candidates[
      file.exists(file.path(source_candidates, "DESCRIPTION"))
    ][1]

    package_build_root <- file.path(TEMP_ROOT, "package_build")
    patched_source <- file.path(package_build_root, "SLCPTAC_qs2")
    unlink(package_build_root, recursive = TRUE, force = TRUE)
    dir.create(package_build_root, recursive = TRUE, showWarnings = FALSE)

    if (!is.na(existing_source) && nzchar(existing_source)) {
      file.copy(existing_source, dirname(patched_source), recursive = TRUE)
      copied_source <- file.path(dirname(patched_source), basename(existing_source))
      if (!identical(normalizePath(copied_source, winslash = "/", mustWork = FALSE),
                     normalizePath(patched_source, winslash = "/", mustWork = FALSE))) {
        if (dir.exists(patched_source)) {
          unlink(patched_source, recursive = TRUE, force = TRUE)
        }
        file.rename(copied_source, patched_source)
      }
    } else {
      git_bin <- Sys.which("git")
      if (!nzchar(git_bin)) {
        stop(
          "SLCPTAC source is not available locally and git is unavailable. ",
          "Please clone https://github.com/SolvingLab/SLCPTAC to temporary/slcptac_source."
        )
      }
      status <- system2(
        git_bin,
        c("clone", "--depth", "1", "https://github.com/SolvingLab/SLCPTAC.git", patched_source),
        stdout = TRUE,
        stderr = TRUE
      )
      if (!file.exists(file.path(patched_source, "DESCRIPTION"))) {
        stop("Failed to clone SLCPTAC source: ", paste(status, collapse = "\n"))
      }
    }

    description_file <- file.path(patched_source, "DESCRIPTION")
    description_lines <- readLines(description_file, warn = FALSE)
    description_lines <- gsub("qs \\(>= 0\\.25\\.0\\),\\s*", "qs2 (>= 0.2.0), ", description_lines)
    description_lines <- gsub("Imports:\\s*qs,", "Imports: qs2,", description_lines)
    writeLines(description_lines, description_file, useBytes = TRUE)

    r_files <- list.files(file.path(patched_source, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
    for (r_file in r_files) {
      code <- readLines(r_file, warn = FALSE)
      patched_code <- gsub("qs::qread\\(", "qs2::qs_read(", code, fixed = FALSE)
      if (!identical(code, patched_code)) {
        writeLines(patched_code, r_file, useBytes = TRUE)
      }
    }

    patched_source
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(
    c("fgsea", "ComplexHeatmap", "limma", "AnnotationDbi", "org.Hs.eg.db"),
    ask = FALSE,
    update = FALSE
  )

  install_cran_if_missing(
    c("dplyr", "ggplot2", "survival", "patchwork", "jsonlite",
      "Hmisc", "tidyr", "ggnewscale", "ggrepel", "scales", "circlize",
      "igraph", "ggraph", "Cairo", "png", "pdftools", "qs2")
  )
  ensure_qs2_dependency()

  remotes::install_github("GangLiLab/geneset", upgrade = "never")
  if (!requireNamespace("genekitr", quietly = TRUE)) {
    remotes::install_github("GangLiLab/genekitr", upgrade = "never")
  }
  remotes::install_github("Zaoqu-Liu/genekitr2", upgrade = "never")
  patched_source <- prepare_slcptac_qs2_source()
  unlink(file.path(.libPaths()[1], "00LOCK-SLCPTAC"), recursive = TRUE, force = TRUE)
  remotes::install_local(patched_source, upgrade = "never", dependencies = FALSE)
  invisible(TRUE)
}

if (AUTO_INSTALL_SLCPTAC) {
  install_slcptac_if_requested()
}

required_packages <- c(
  "SLCPTAC", "ggplot2", "grid", "Cairo", "png", "qs2", "genekitr",
  "AnnotationDbi", "org.Hs.eg.db"
)
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
    "\nThis project uses a qs2-patched SLCPTAC install path because qs is not compatible with this R 4.6 build."
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

create_slcptac_smoke_bulk_data <- function(output_root) {
  if (!SMOKE_TEST_MODE) {
    return(invisible(FALSE))
  }

  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop("SLCPTAC smoke test requires package 'qs2'.")
  }

  set.seed(20260608)
  output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  split_root <- file.path(output_root, "CPTAC_Omics_Split")

  unlink(output_root, recursive = TRUE, force = TRUE)
  dir.create(split_root, recursive = TRUE, showWarnings = FALSE)

  cancers <- unique(PAN_CANCERS)
  n_per_cancer <- max(SMOKE_TEST_SAMPLES_PER_CANCER, 16L)
  sample_table <- do.call(rbind, lapply(cancers, function(cancer) {
    data.frame(
      sampleid = sprintf("%s_SMK%03d", cancer, seq_len(n_per_cancer)),
      cancer = cancer,
      sample_index = seq_len(n_per_cancer),
      stringsAsFactors = FALSE
    )
  }))
  rownames(sample_table) <- sample_table$sampleid

  sample_table$type <- "Tumor"
  sample_table$age <- 45 + (sample_table$sample_index %% 28) + match(sample_table$cancer, cancers)
  sample_table$bmi <- 20 + (sample_table$sample_index %% 11) / 1.5
  sample_table$gender <- ifelse(sample_table$sample_index %% 2 == 0, "Female", "Male")
  sample_table$tumor_stage <- paste0("Stage ", c("I", "II", "III", "IV")[(sample_table$sample_index %% 4) + 1L])
  sample_table$os_time <- 240 + sample_table$sample_index * 18 + match(sample_table$cancer, cancers) * 7
  sample_table$os_event <- as.integer(sample_table$sample_index %% 3 == 0)
  sample_table$pfs_time <- 160 + sample_table$sample_index * 13 + match(sample_table$cancer, cancers) * 5
  sample_table$pfs_event <- as.integer(sample_table$sample_index %% 4 == 0)

  required_gene_symbols <- unique(c(
    TARGET_GENES,
    CONTINUOUS_PANEL_GENES,
    MUTATION_PANEL_GENES,
    PRIMARY_MUTATION_GENE,
    SECONDARY_MUTATION_GENES,
    AUXILIARY_PHOSPHO_GENE,
    "TP53"
  ))
  required_gene_symbols <- required_gene_symbols[nzchar(required_gene_symbols)]

  geneset_symbols <- tryCatch({
    gs <- SLCPTAC:::.load_geneset_df(enrich_type = "MsigDB", msigdb_category = MSIGDB_CATEGORY)
    unique(as.character(gs$gene))
  }, error = function(error) {
    c(
      "ABCA1", "ACADM", "ACLY", "AKT1", "ALDOA", "APOE", "ARAF", "BCL2",
      "BCL2L13", "BRAF", "CASP3", "CCND1", "CDK1", "CDK2", "CDKN1A",
      "CEBPB", "DDIT3", "EGFR", "EGR1", "EIF4E", "FOS", "GAPDH",
      "HIF1A", "JUN", "KRAS", "MAPK1", "MDM2", "MTOR", "MYC", "NFKB1",
      "PIK3CA", "PTEN", "RPS6", "STAT3", "TGFB1", "TP53", "VEGFA"
    )
  })
  genome_gene_symbols <- unique(c(required_gene_symbols, geneset_symbols))
  genome_gene_symbols <- genome_gene_symbols[nzchar(genome_gene_symbols)]

  make_numeric_signal <- function(gene, modal_shift = 0) {
    gene_id <- match(gene, genome_gene_symbols)
    cancer_id <- match(sample_table$cancer, cancers)
    index <- sample_table$sample_index
    signal <- sin(index / 3 + gene_id / 7) +
      cos(cancer_id / 2 + gene_id / 11) +
      modal_shift +
      rnorm(nrow(sample_table), sd = 0.35)
    as.numeric(scale(signal))
  }

  make_mutation_status <- function(gene) {
    gene_id <- match(gene, genome_gene_symbols)
    is_mut <- (sample_table$sample_index + gene_id + match(sample_table$cancer, cancers)) %% 4 == 0
    ifelse(is_mut, "Mutation", "WildType")
  }

  clinical_cols <- c(
    "sampleid", "type", "age", "bmi", "gender", "tumor_stage",
    "os_time", "os_event", "pfs_time", "pfs_event"
  )

  for (gene in required_gene_symbols) {
    gene_data <- sample_table[, clinical_cols, drop = FALSE]
    gene_data[[paste0(gene, "_RNAseq")]] <- make_numeric_signal(gene, modal_shift = 0)
    gene_data[[paste0(gene, "_Protein")]] <- 0.55 * gene_data[[paste0(gene, "_RNAseq")]] +
      make_numeric_signal(gene, modal_shift = 0.35) * 0.45
    gene_data[[paste0(gene, "_logCNA")]] <- make_numeric_signal(gene, modal_shift = -0.15)
    gene_data[[paste0(gene, "_Methylation")]] <- -0.4 * gene_data[[paste0(gene, "_RNAseq")]] +
      make_numeric_signal(gene, modal_shift = 0.2) * 0.6
    gene_data[[paste0(gene, "_Mutation")]] <- make_mutation_status(gene)

    qs2::qs_save(
      gene_data,
      file.path(split_root, paste0(gene, "_cptac.qs"))
    )
  }

  make_genome_matrix <- function(modal, modal_shift = 0) {
    mat <- vapply(genome_gene_symbols, function(gene) {
      make_numeric_signal(gene, modal_shift = modal_shift)
    }, numeric(nrow(sample_table)))
    mat <- as.data.frame(mat, check.names = FALSE)
    rownames(mat) <- sample_table$sampleid
    mat
  }

  protein_matrix <- make_genome_matrix("Protein", modal_shift = 0.3)
  rna_matrix <- make_genome_matrix("RNAseq", modal_shift = 0)
  logcna_matrix <- make_genome_matrix("logCNA", modal_shift = -0.2)
  methylation_matrix <- make_genome_matrix("Methylation", modal_shift = 0.1)

  kras_mut <- make_mutation_status(PRIMARY_MUTATION_GENE)
  affected_genes <- intersect(colnames(protein_matrix), head(geneset_symbols, 80))
  protein_matrix[kras_mut == "Mutation", affected_genes] <-
    protein_matrix[kras_mut == "Mutation", affected_genes, drop = FALSE] + 1.2
  rna_matrix[kras_mut == "Mutation", affected_genes] <-
    rna_matrix[kras_mut == "Mutation", affected_genes, drop = FALSE] + 0.8

  mutation_matrix <- as.data.frame(vapply(genome_gene_symbols, function(gene) {
    make_mutation_status(gene)
  }, character(nrow(sample_table))), check.names = FALSE)
  rownames(mutation_matrix) <- sample_table$sampleid

  phospho_matrix <- data.frame(.row_id = seq_len(nrow(sample_table)))
  rownames(phospho_matrix) <- sample_table$sampleid
  phospho_matrix$.row_id <- NULL
  for (site in c("S124", "S126", "S129", "T308", "S473")) {
    phospho_matrix[[paste0(site, "_", AUXILIARY_PHOSPHO_GENE)]] <-
      make_numeric_signal(AUXILIARY_PHOSPHO_GENE, modal_shift = runif(1, -0.2, 0.2))
  }

  clinical_bulk <- sample_table[, clinical_cols, drop = FALSE]
  qs2::qs_save(rna_matrix, file.path(output_root, "LinkedOmicsKB_PanCancer_RNAseq_RSEM.qs"))
  qs2::qs_save(protein_matrix, file.path(output_root, "LinkedOmicsKB_PanCancer_Protein_Quantification.qs"))
  qs2::qs_save(phospho_matrix, file.path(output_root, "LinkedOmicsKB_PanCancer_Phospho_Quantification.qs"))
  qs2::qs_save(methylation_matrix, file.path(output_root, "LinkedOmicsKB_PanCancer_Methylation.qs"))
  qs2::qs_save(logcna_matrix, file.path(output_root, "LinkedOmicsKB_PanCancer_CNV_logCNA.qs"))
  qs2::qs_save(mutation_matrix, file.path(output_root, "LinkedOmicsKB_PanCancer_Mutation_Binary.qs"))
  qs2::qs_save(clinical_bulk, file.path(output_root, "LinkedOmicsKB_PanCancer_Clin.qs"))

  manifest <- data.frame(
    Smoke_Test = TRUE,
    Bulk_Data_Root = output_root,
    Cancers = paste(cancers, collapse = ";"),
    Samples = nrow(sample_table),
    Gene_Files = length(required_gene_symbols),
    Genome_Genes = length(genome_gene_symbols),
    Generated_At = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    manifest,
    file = file.path(TABLE_ROOT, "000_slcptac_smoke_test_bulk_manifest.csv"),
    row.names = FALSE,
    na = "NA"
  )

  invisible(TRUE)
}

LINKEDOMICS_BASE_URL <- "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized"
LINKEDOMICS_CANCER_CODES <- c(
  BRCA = "BRCA", CCRCC = "CCRCC", COAD = "COAD", GBM = "GBM", HNSCC = "HNSCC",
  LUAD = "LUAD", LUSC = "LSCC", OV = "OV", PDAC = "PDAC", UCEC = "UCEC"
)
LINKEDOMICS_RESOURCE_SUFFIX <- c(
  RNAseq = "RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
  Protein = "proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
  logCNA = "WES_CNV_gene_ratio_log2.txt",
  Mutation = "somatic_mutation_gene_level_binary.txt",
  Meta = "meta.txt",
  Survival = "survival.txt"
)

get_required_gene_symbols_for_linkedomics <- function(tasks) {
  genes <- unique(c(
    TARGET_GENES,
    CONTINUOUS_PANEL_GENES,
    MUTATION_PANEL_GENES,
    PRIMARY_MUTATION_GENE,
    SECONDARY_MUTATION_GENES,
    "TP53",
    unlist(lapply(tasks, extract_task_genes), use.names = FALSE)
  ))
  genes <- unique(trimws(as.character(genes)))
  genes[nzchar(genes)]
}

linkedomics_code_for_cancer <- function(cancer) {
  code <- unname(LINKEDOMICS_CANCER_CODES[[cancer]])
  if (is.null(code) || is.na(code)) {
    stop("No LinkedOmicsKB download code configured for cancer: ", cancer)
  }
  code
}

linkedomics_url <- function(cancer, resource) {
  code <- linkedomics_code_for_cancer(cancer)
  suffix <- unname(LINKEDOMICS_RESOURCE_SUFFIX[[resource]])
  if (is.null(suffix) || is.na(suffix)) {
    stop("No LinkedOmicsKB resource suffix configured for: ", resource)
  }
  sprintf("%s/%s/%s_%s", LINKEDOMICS_BASE_URL, code, code, suffix)
}

linkedomics_raw_file <- function(cancer, resource) {
  code <- linkedomics_code_for_cancer(cancer)
  file.path(LINKEDOMICS_RAW_ROOT, code, basename(linkedomics_url(cancer, resource)))
}

download_linkedomics_file <- function(cancer, resource, overwrite = FALSE) {
  destination <- linkedomics_raw_file(cancer, resource)
  if (file.exists(destination) && file.info(destination)$size > 0 && !overwrite) {
    return(destination)
  }

  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  url <- linkedomics_url(cancer, resource)
  if (QUICKANALYSIS_VERBOSE) {
    cat("Downloading LinkedOmicsKB ", cancer, " ", resource, "...\n", sep = "")
  }

  ok <- FALSE
  last_error <- NULL
  for (attempt in seq_len(4L)) {
    last_error <- tryCatch({
      utils::download.file(
        url = url,
        destfile = destination,
        mode = "wb",
        method = "libcurl",
        quiet = TRUE,
        headers = c("User-Agent" = "Mozilla/5.0")
      )
      NULL
    }, error = function(error) error)
    ok <- is.null(last_error) && file.exists(destination) && file.info(destination)$size > 0
    if (ok) {
      break
    }
    Sys.sleep(1 + attempt)
  }

  if (!ok) {
    unlink(destination, force = TRUE)
    stop(
      "Failed to download LinkedOmicsKB file: ", url,
      if (!is.null(last_error)) paste0("\n", conditionMessage(last_error)) else ""
    )
  }

  destination
}

.ensembl_symbol_cache <- new.env(parent = emptyenv())

map_ensembl_to_symbol <- function(ensembl_ids) {
  ensembl_base <- sub("\\..*$", "", as.character(ensembl_ids))
  unique_ids <- unique(ensembl_base[nzchar(ensembl_base)])
  missing_ids <- setdiff(unique_ids, ls(.ensembl_symbol_cache, all.names = TRUE))

  if (length(missing_ids) > 0L) {
    mapped <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = missing_ids,
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    ))
    mapped <- mapped[!is.na(mapped$SYMBOL) & nzchar(mapped$SYMBOL), , drop = FALSE]
    mapped <- mapped[!duplicated(mapped$ENSEMBL), , drop = FALSE]
    symbol_map <- setNames(mapped$SYMBOL, mapped$ENSEMBL)

    for (ensembl_id in missing_ids) {
      symbol <- unname(symbol_map[ensembl_id])
      if (length(symbol) == 0L || is.na(symbol)) {
        symbol <- NA_character_
      }
      assign(
        ensembl_id,
        symbol,
        envir = .ensembl_symbol_cache
      )
    }
  }

  unname(vapply(
    ensembl_base,
    function(ensembl_id) {
      if (!nzchar(ensembl_id) || !exists(ensembl_id, envir = .ensembl_symbol_cache, inherits = FALSE)) {
        return(NA_character_)
      }
      get(ensembl_id, envir = .ensembl_symbol_cache, inherits = FALSE)
    },
    character(1)
  ))
}

read_linkedomics_matrix <- function(file, cancer, wanted_symbols = NULL, binary = FALSE) {
  dat <- utils::read.delim(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (ncol(dat) < 2L) {
    stop("Downloaded matrix has too few columns: ", file)
  }

  row_ids <- dat[[1]]
  symbols <- map_ensembl_to_symbol(row_ids)
  keep <- !is.na(symbols) & nzchar(symbols)
  if (!is.null(wanted_symbols)) {
    keep <- keep & symbols %in% wanted_symbols
  }
  dat <- dat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (nrow(dat) == 0L) {
    return(data.frame(row.names = character(0)))
  }

  value_df <- dat[, -1, drop = FALSE]
  value_df[] <- lapply(value_df, function(x) suppressWarnings(as.numeric(x)))
  value_df$SYMBOL <- symbols
  value_df <- stats::aggregate(. ~ SYMBOL, data = value_df, FUN = function(x) {
    if (binary) {
      as.numeric(any(x > 0, na.rm = TRUE))
    } else {
      mean(x, na.rm = TRUE)
    }
  })
  rownames(value_df) <- value_df$SYMBOL
  value_df$SYMBOL <- NULL

  matrix_df <- as.data.frame(t(as.matrix(value_df)), check.names = FALSE)
  rownames(matrix_df) <- paste0(cancer, "_", rownames(matrix_df))
  matrix_df
}

read_linkedomics_meta <- function(file, cancer) {
  meta <- utils::read.delim(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if ("case_id" %in% names(meta)) {
    meta <- meta[meta$case_id != "data_type", , drop = FALSE]
  }
  sample_raw <- meta$case_id
  out <- data.frame(
    sampleid = paste0(cancer, "_", sample_raw),
    type = "Tumor",
    age = suppressWarnings(as.numeric(meta[["Age"]] %||% NA)),
    bmi = suppressWarnings(as.numeric(meta[["BMI"]] %||% NA)),
    gender = as.character(meta[["Sex"]] %||% NA),
    tumor_stage = as.character((meta[["Stage"]] %||% meta[["tumor_stage"]] %||% NA)),
    stringsAsFactors = FALSE
  )
  out
}

read_linkedomics_survival <- function(file, cancer) {
  survival <- utils::read.delim(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  data.frame(
    sampleid = paste0(cancer, "_", survival$case_id),
    os_time = suppressWarnings(as.numeric(survival$OS_days)),
    os_event = suppressWarnings(as.integer(survival$OS_event)),
    pfs_time = suppressWarnings(as.numeric(survival$PFS_days)),
    pfs_event = suppressWarnings(as.integer(survival$PFS_event)),
    stringsAsFactors = FALSE
  )
}

merge_named_vector <- function(target, values, column) {
  if (length(values) == 0L || is.null(names(values))) {
    target[[column]] <- NA_real_
    return(target)
  }
  target[[column]] <- unname(values[target$sampleid])
  target
}

flatten_sample_named_values <- function(modal_value_list) {
  # unlist() 会把外层癌种名拼进 names，导致 COAD_01CO001 变成 COAD.COAD_01CO001。
  # SLCPTAC 原包按 sampleid 精确匹配，所以这里必须只保留内层样本名。
  modal_value_list <- modal_value_list[lengths(modal_value_list) > 0L]
  if (length(modal_value_list) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }

  values <- unlist(unname(modal_value_list), use.names = TRUE)
  sample_names <- names(values)
  if (is.null(sample_names)) {
    return(values)
  }

  duplicated_samples <- duplicated(sample_names)
  if (any(duplicated_samples)) {
    values <- values[!duplicated_samples]
  }
  values
}

empty_omics_frame <- function() {
  data.frame(check.names = FALSE)
}

prepare_linkedomics_bulk_data <- function(tasks) {
  if (SMOKE_TEST_MODE || !AUTO_PREPARE_CPTAC_BULK) {
    return(invisible(FALSE))
  }

  required_files <- make_required_file_table(tasks)
  manifest_ok <- FALSE
  if (file.exists(LINKEDOMICS_BULK_MANIFEST_FILE)) {
    old_manifest <- tryCatch(
      utils::read.csv(LINKEDOMICS_BULK_MANIFEST_FILE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(error) data.frame()
    )
    manifest_ok <- nrow(old_manifest) > 0L &&
      "Converter_Version" %in% colnames(old_manifest) &&
      identical(old_manifest$Converter_Version[1], LINKEDOMICS_BULK_CONVERTER_VERSION)
  }
  if (nrow(required_files) > 0L && all(required_files$Exists) && manifest_ok) {
    return(invisible(FALSE))
  }

  if (QUICKANALYSIS_VERBOSE) {
    cat("\nPreparing SLCPTAC bulk data from LinkedOmicsKB public downloads...\n")
  }
  dir.create(SLCPTAC_BULK_DATA_ROOT, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(SLCPTAC_BULK_DATA_ROOT, "CPTAC_Omics_Split"), recursive = TRUE, showWarnings = FALSE)
  dir.create(LINKEDOMICS_RAW_ROOT, recursive = TRUE, showWarnings = FALSE)

  cancers_for_gene_files <- unique(c(PAN_CANCERS, TARGET_CANCERS))
  cancers_for_gene_files <- intersect(cancers_for_gene_files, names(LINKEDOMICS_CANCER_CODES))
  target_cancers_for_genome <- unique(TARGET_CANCERS)
  required_genes <- get_required_gene_symbols_for_linkedomics(tasks)

  clinical_list <- list()
  survival_list <- list()
  gene_modal_values <- setNames(vector("list", length(required_genes)), required_genes)
  for (gene in required_genes) {
    gene_modal_values[[gene]] <- list(RNAseq = list(), Protein = list(), logCNA = list(), Mutation = list())
  }

  coad_genome_matrices <- list()

  for (cancer in cancers_for_gene_files) {
    meta_file <- download_linkedomics_file(cancer, "Meta")
    survival_file <- download_linkedomics_file(cancer, "Survival")
    clinical_list[[cancer]] <- read_linkedomics_meta(meta_file, cancer)
    survival_list[[cancer]] <- read_linkedomics_survival(survival_file, cancer)

    for (modal in c("RNAseq", "Protein")) {
      raw_file <- download_linkedomics_file(cancer, modal)
      subset_matrix <- read_linkedomics_matrix(raw_file, cancer, wanted_symbols = required_genes)
      for (gene in intersect(required_genes, colnames(subset_matrix))) {
        gene_modal_values[[gene]][[modal]][[cancer]] <- setNames(subset_matrix[[gene]], rownames(subset_matrix))
      }
      if (cancer %in% target_cancers_for_genome) {
        coad_genome_matrices[[modal]] <- read_linkedomics_matrix(raw_file, cancer)
      }
    }
  }

  for (cancer in target_cancers_for_genome) {
    for (modal in c("logCNA", "Mutation")) {
      raw_file <- download_linkedomics_file(cancer, modal)
      matrix <- read_linkedomics_matrix(raw_file, cancer, wanted_symbols = if (modal == "Mutation") required_genes else NULL, binary = modal == "Mutation")
      for (gene in intersect(required_genes, colnames(matrix))) {
        gene_modal_values[[gene]][[modal]][[cancer]] <- setNames(matrix[[gene]], rownames(matrix))
      }
      coad_genome_matrices[[modal]] <- matrix
    }
  }

  clinical <- do.call(rbind, unname(clinical_list))
  survival <- do.call(rbind, unname(survival_list))
  clinical <- merge(clinical, survival, by = "sampleid", all = TRUE)
  clinical$type <- ifelse(is.na(clinical$type), "Tumor", clinical$type)

  for (gene in required_genes) {
    gene_data <- clinical
    rna_values <- flatten_sample_named_values(gene_modal_values[[gene]]$RNAseq)
    protein_values <- flatten_sample_named_values(gene_modal_values[[gene]]$Protein)
    logcna_values <- flatten_sample_named_values(gene_modal_values[[gene]]$logCNA)
    mutation_values <- flatten_sample_named_values(gene_modal_values[[gene]]$Mutation)

    gene_data <- merge_named_vector(gene_data, rna_values, paste0(gene, "_RNAseq"))
    gene_data <- merge_named_vector(gene_data, protein_values, paste0(gene, "_Protein"))
    gene_data <- merge_named_vector(gene_data, logcna_values, paste0(gene, "_logCNA"))
    gene_data[[paste0(gene, "_Methylation")]] <- NA_real_
    mutation_status <- rep(NA_real_, nrow(gene_data))
    if (length(mutation_values) > 0L && !is.null(names(mutation_values))) {
      mutation_status <- unname(mutation_values[gene_data$sampleid])
    }
    mutation_label <- ifelse(
      is.na(mutation_status),
      NA_character_,
      ifelse(as.numeric(mutation_status) > 0, "Mutation", "WildType")
    )
    target_sample <- grepl(
      paste0("^(", paste(TARGET_CANCERS, collapse = "|"), ")_"),
      gene_data$sampleid
    )
    mutation_label[target_sample & is.na(mutation_label)] <- "WildType"
    gene_data[[paste0(gene, "_Mutation")]] <- mutation_label

    qs2::qs_save(
      gene_data,
      file.path(SLCPTAC_BULK_DATA_ROOT, "CPTAC_Omics_Split", paste0(gene, "_cptac.qs"))
    )
  }

  qs2::qs_save(coad_genome_matrices$RNAseq %||% empty_omics_frame(), bulk_genome_file("RNAseq"))
  qs2::qs_save(coad_genome_matrices$Protein %||% empty_omics_frame(), bulk_genome_file("Protein"))
  qs2::qs_save(empty_omics_frame(), bulk_genome_file("Phospho"))
  qs2::qs_save(empty_omics_frame(), bulk_genome_file("Methylation"))
  qs2::qs_save(coad_genome_matrices$logCNA %||% empty_omics_frame(), bulk_genome_file("logCNA"))
  qs2::qs_save(coad_genome_matrices$Mutation %||% empty_omics_frame(), bulk_genome_file("Mutation"))
  qs2::qs_save(clinical, bulk_genome_file("Clinical"))

  manifest <- data.frame(
    Source = "LinkedOmicsKB CPTAC pan-cancer public S3 downloads",
    Source_Page = "https://kb.linkedomics.org/download",
    Converter_Version = LINKEDOMICS_BULK_CONVERTER_VERSION,
    Raw_Root = LINKEDOMICS_RAW_ROOT,
    Bulk_Data_Root = SLCPTAC_BULK_DATA_ROOT,
    Cancers = paste(cancers_for_gene_files, collapse = ";"),
    Target_Cancers_For_Genome = paste(target_cancers_for_genome, collapse = ";"),
    Gene_Files = length(required_genes),
    Generated_At = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    manifest,
    file = LINKEDOMICS_BULK_MANIFEST_FILE,
    row.names = FALSE,
    na = "NA"
  )
  write_table(manifest, TABLE_ROOT, "000_slcptac_linkedomics_bulk_manifest")
  invisible(TRUE)
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
  add_row("scenario01", "cptac_correlation", "相关/关联场景1", "1 continuous vs 1 continuous", "Pearson/Spearman/Kendall", "CorPlot或LollipopPlot", "ATF3 RNAseq vs logCNA", "当前LinkedOmicsKB COAD蛋白矩阵缺少ATF3蛋白，默认用ATF3 RNAseq/logCNA完成主线")
  add_row("scenario02", "cptac_correlation", "相关/关联场景2", "1 continuous vs multiple continuous", "相关分析+FDR", "LollipopPlot或DotPlot", "ATF3 RNAseq vs连续变量logCNA面板", "")
  add_row("scenario03", "cptac_correlation", "相关/关联场景3", "multiple continuous vs multiple continuous", "相关矩阵+FDR", "DotPlot", "连续变量面板RNAseq vs logCNA", "")
  add_row("scenario04", "cptac_correlation", "相关/关联场景4", "1 categorical vs 1 continuous", "Wilcoxon或Kruskal-Wallis", "BoxPlot", "KRAS Mutation/Stage vs ATF3 RNAseq", "")
  add_row("scenario05", "cptac_correlation", "相关/关联场景5", "1 continuous vs multiple categorical", "多组Wilcoxon/Kruskal-Wallis", "Multiple BoxPlots", "ATF3 RNAseq vs突变面板", "")
  add_row("scenario06", "cptac_correlation", "相关/关联场景6", "multiple continuous vs 1 categorical", "多组Wilcoxon/Kruskal-Wallis", "Multiple BoxPlots", "连续变量RNAseq面板 vs KRAS Mutation", "")
  add_row("scenario07", "cptac_correlation", "相关/关联场景7", "categorical vs categorical", "Fisher或Chi-square+Odds Ratio", "BarPlot或Heatmap", "共突变/临床-突变关联", "")
  add_row("scenario08", "cptac_enrichment", "富集场景8", "1 categorical vs genome-wide", "limma DEA", "NetworkPlot", "KRAS Mutation驱动蛋白组扫描", "")
  add_row("scenario09", "cptac_enrichment", "富集场景9", "1 categorical vs pathways", "limma DEA -> fgsea", "GSEA DotPlot", "KRAS Mutation驱动通路富集", "")
  add_row("scenario10", "cptac_enrichment", "富集场景10", "multiple categorical vs genome-wide", "多变量limma DEA", "DotPlot Paired", "突变面板驱动蛋白组扫描", "")
  add_row("scenario11", "cptac_enrichment", "富集场景11", "multiple categorical vs pathways", "多变量DEA -> fgsea", "GSEA Matrix", "突变面板通路富集", "")
  add_row("scenario12", "cptac_enrichment", "富集场景12", "1 continuous vs genome-wide", "相关排序", "NetworkPlot", "ATF3 RNAseq相关蛋白组扫描", "var1使用ATF3 RNAseq，genome_modal仍使用CPTAC COAD Protein矩阵")
  add_row("scenario13", "cptac_enrichment", "富集场景13", "1 continuous vs pathways", "相关排序 -> fgsea", "GSEA DotPlot", "ATF3 RNAseq相关通路富集", "var1使用ATF3 RNAseq，GSEA排名来自Protein genome")
  add_row("scenario14", "cptac_enrichment", "富集场景14", "multiple continuous vs genome-wide", "多连续变量相关排序", "DotPlot Paired", "连续变量RNAseq面板相关蛋白组扫描", "")
  add_row("scenario15", "cptac_enrichment", "富集场景15", "multiple continuous vs pathways", "多连续变量相关排序 -> fgsea", "GSEA Matrix", "连续变量RNAseq面板通路富集", "")
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
  missing_table <- if (length(missing_from_catalog) > 0L) {
    data.frame(
      Function = missing_from_catalog,
      Status = rep(
        "SLCPTAC导出但不是单独分析场景函数，或已由上层场景间接覆盖",
        length(missing_from_catalog)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Function = character(0),
      Status = character(0),
      stringsAsFactors = FALSE
    )
  }
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
  line_count <- max(1L, ceiling(nchar(title, type = "width") / 68))
  0.55 + 0.24 * (line_count - 1L)
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
          gene, "logCNA", cancer,
          gene, cancer,
          paste0(gene, " RNAseq vs logCNA in CPTAC-", cancer),
          width = 5.4, height = 4.8
        )
        tasks <- add_correlation_task(
          tasks, "scenario01", 1L,
          gene, "RNAseq", pan_common,
          PRIMARY_MUTATION_GENE, "RNAseq", pan_common,
          gene, "pan_cancer_RNAseq_vs_KRAS_RNAseq",
          paste0(gene, " RNAseq vs ", PRIMARY_MUTATION_GENE, " RNAseq across CPTAC pan-cancer"),
          width = 8, height = 4.5
        )
      }

      if ("scenario02" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario02", 2L,
          gene, "RNAseq", cancer,
          CONTINUOUS_PANEL_GENES, "logCNA", cancer,
          gene, paste0(cancer, "_RNAseq_vs_panel_logCNA"),
          paste0(gene, " RNAseq vs logCNA panel in CPTAC-", cancer),
          width = 7, height = 4.8
        )
      }

      if ("scenario03" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario03", 3L,
          CONTINUOUS_PANEL_GENES, "RNAseq", cancer,
          CONTINUOUS_PANEL_GENES, "logCNA", cancer,
          gene, paste0(cancer, "_panel_RNAseq_vs_logCNA"),
          paste0("RNAseq vs logCNA panel correlation in CPTAC-", cancer),
          width = 7, height = 5.6
        )
      }

      if ("scenario04" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario04", 4L,
          PRIMARY_MUTATION_GENE, "Mutation", cancer,
          gene, "RNAseq", cancer,
          gene, paste0(cancer, "_", PRIMARY_MUTATION_GENE, "_Mutation_vs_RNAseq"),
          paste0(PRIMARY_MUTATION_GENE, " Mutation vs ", gene, " RNAseq in CPTAC-", cancer),
          width = 5.2, height = 4.8
        )
        tasks <- add_correlation_task(
          tasks, "scenario04", 4L,
          PRIMARY_CLINICAL_VARIABLE, "Clinical", cancer,
          gene, "RNAseq", cancer,
          gene, paste0(cancer, "_", PRIMARY_CLINICAL_VARIABLE, "_vs_RNAseq"),
          paste0(PRIMARY_CLINICAL_VARIABLE, " vs ", gene, " RNAseq in CPTAC-", cancer),
          width = 5.2, height = 4.8
        )
      }

      if ("scenario05" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario05", 5L,
          gene, "RNAseq", cancer,
          MUTATION_PANEL_GENES, "Mutation", cancer,
          gene, paste0(cancer, "_RNAseq_vs_mutation_panel"),
          paste0(gene, " RNAseq vs mutation panel in CPTAC-", cancer),
          width = 9, height = 4.8
        )
      }

      if ("scenario06" %in% selected_analyses) {
        tasks <- add_correlation_task(
          tasks, "scenario06", 6L,
          CONTINUOUS_PANEL_GENES, "RNAseq", cancer,
          PRIMARY_MUTATION_GENE, "Mutation", cancer,
          gene, paste0(cancer, "_RNAseq_panel_vs_", PRIMARY_MUTATION_GENE, "_Mutation"),
          paste0("RNAseq panel vs ", PRIMARY_MUTATION_GENE, " Mutation in CPTAC-", cancer),
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
            gene, "RNAseq", cancer,
            "genome", genome_modal, "MsigDB",
            gene, paste0(cancer, "_RNAseq_vs_", genome_modal, "_genome"),
            paste0(gene, " RNAseq-related ", genome_modal, " genome scan in CPTAC-", cancer),
            width = 8.5, height = 5.8
          )
        }
        if ("scenario14" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario14", 14L,
            CONTINUOUS_PANEL_GENES, "RNAseq", cancer,
            "genome", genome_modal, "MsigDB",
            paste_compact(CONTINUOUS_PANEL_GENES), paste0(cancer, "_RNAseq_panel_vs_", genome_modal, "_genome"),
            paste0("RNAseq panel-related ", genome_modal, " genome scan in CPTAC-", cancer),
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
            gene, "RNAseq", cancer,
            "enrichment", "Protein", database,
            gene, paste0(cancer, "_RNAseq_GSEA_", database),
            paste0(gene, " RNAseq-related ", database, " GSEA in CPTAC-", cancer),
            width = 8, height = 6.2
          )
        }
        if ("scenario15" %in% selected_analyses) {
          tasks <- add_enrichment_task(
            tasks, "scenario15", 15L,
            CONTINUOUS_PANEL_GENES, "RNAseq", cancer,
            "enrichment", "Protein", database,
            paste_compact(CONTINUOUS_PANEL_GENES), paste0(cancer, "_RNAseq_panel_GSEA_", database),
            paste0("RNAseq panel-related ", database, " GSEA in CPTAC-", cancer),
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
            PRIMARY_MUTATION_GENE, "Protein", cancer,
            surv_type, PRIMARY_MUTATION_GENE, paste0(cancer, "_", PRIMARY_MUTATION_GENE, "_Protein_", surv_type),
            paste0(PRIMARY_MUTATION_GENE, " Protein ", surv_type, " survival in CPTAC-", cancer),
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

.bulk_data_cache <- new.env(parent = emptyenv())

read_bulk_gene_data <- function(gene) {
  cache_key <- paste0("gene::", gene)
  if (exists(cache_key, envir = .bulk_data_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .bulk_data_cache, inherits = FALSE))
  }

  file <- bulk_gene_file(gene)
  dat <- if (file.exists(file)) {
    tryCatch(qs2::qs_read(file), error = function(error) NULL)
  } else {
    NULL
  }
  assign(cache_key, dat, envir = .bulk_data_cache)
  dat
}

read_bulk_genome_data <- function(modal) {
  cache_key <- paste0("genome::", modal)
  if (exists(cache_key, envir = .bulk_data_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .bulk_data_cache, inherits = FALSE))
  }

  file <- bulk_genome_file(modal)
  dat <- if (file.exists(file)) {
    tryCatch(qs2::qs_read(file), error = function(error) NULL)
  } else {
    NULL
  }
  assign(cache_key, dat, envir = .bulk_data_cache)
  dat
}

get_cancer_prefix <- function(sample_ids) {
  sapply(strsplit(as.character(sample_ids), "_", fixed = TRUE), `[`, 1)
}

subset_gene_data_by_cancers <- function(dat, cancers) {
  if (is.null(dat) || !is.data.frame(dat) || !"sampleid" %in% colnames(dat)) {
    return(data.frame())
  }
  cancer_type <- get_cancer_prefix(dat$sampleid)
  keep <- cancer_type %in% cancers
  if ("type" %in% colnames(dat)) {
    keep <- keep & (is.na(dat$type) | dat$type == "Tumor")
  }
  dat[keep, , drop = FALSE]
}

has_enough_numeric_values <- function(values, min_non_na = 3L) {
  values <- suppressWarnings(as.numeric(as.character(values)))
  values <- values[!is.na(values)]
  length(values) >= min_non_na && length(unique(values)) >= 2L
}

has_enough_categorical_values <- function(values, min_non_na = 3L) {
  values <- values[!is.na(values) & nzchar(as.character(values))]
  length(values) >= min_non_na && length(unique(values)) >= 2L
}

match_clinical_column_for_script <- function(search, available) {
  search_lower <- tolower(search)
  available_lower <- tolower(available)

  direct <- which(available_lower == search_lower)
  if (length(direct) > 0L) {
    return(available[direct[1]])
  }

  aliases <- list(
    age = c("age"),
    bmi = c("bmi", "body_mass_index"),
    gender = c("gender", "sex"),
    sex = c("sex", "gender"),
    tumor_stage = c("tumor_stage", "stage", "ajcc_pathologic_stage"),
    stage = c("stage", "tumor_stage", "ajcc_pathologic_stage")
  )
  for (standard in names(aliases)) {
    if (search_lower %in% aliases[[standard]]) {
      hit <- which(available_lower %in% aliases[[standard]])
      if (length(hit) > 0L) {
        return(available[hit[1]])
      }
    }
  }

  NULL
}

has_usable_gene_modal <- function(gene, modal, cancers) {
  dat <- subset_gene_data_by_cancers(read_bulk_gene_data(gene), cancers)
  if (nrow(dat) == 0L) {
    return(FALSE)
  }

  column <- paste0(gene, "_", modal)
  if (!column %in% colnames(dat)) {
    return(FALSE)
  }

  if (identical(modal, "Mutation")) {
    return(has_enough_categorical_values(dat[[column]]))
  }

  has_enough_numeric_values(dat[[column]])
}

has_usable_clinical_variable <- function(variable, cancers) {
  dat <- subset_gene_data_by_cancers(read_bulk_gene_data("TP53"), cancers)
  if (nrow(dat) == 0L) {
    return(FALSE)
  }
  clinical_columns <- setdiff(
    colnames(dat),
    grep(
      "_RNAseq$|_Protein$|_Methylation$|_logCNA$|_Mutation$|^sampleid$|^type$|^os_|^pfs_",
      colnames(dat),
      value = TRUE
    )
  )
  matched <- match_clinical_column_for_script(variable, clinical_columns)
  if (is.null(matched)) {
    return(FALSE)
  }
  has_enough_categorical_values(dat[[matched]])
}

has_usable_survival_data <- function(cancers, surv_type = "OS") {
  dat <- subset_gene_data_by_cancers(read_bulk_gene_data("TP53"), cancers)
  if (nrow(dat) == 0L) {
    return(FALSE)
  }
  time_col <- if (identical(surv_type, "PFS")) "pfs_time" else "os_time"
  event_col <- if (identical(surv_type, "PFS")) "pfs_event" else "os_event"
  if (!all(c(time_col, event_col) %in% colnames(dat))) {
    return(FALSE)
  }
  valid <- !is.na(dat[[time_col]]) & !is.na(dat[[event_col]])
  sum(valid) >= 10L && length(unique(dat[[event_col]][valid])) >= 2L
}

has_usable_genome_modal <- function(modal, cancers) {
  dat <- read_bulk_genome_data(modal)
  if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0L || ncol(dat) == 0L) {
    return(FALSE)
  }
  cancer_type <- get_cancer_prefix(rownames(dat))
  dat <- dat[cancer_type %in% cancers, , drop = FALSE]
  if (nrow(dat) == 0L) {
    return(FALSE)
  }
  usable_columns <- vapply(dat, has_enough_numeric_values, logical(1), min_non_na = 3L)
  sum(usable_columns) >= 10L
}

has_usable_modal_set <- function(vars, modal, cancers, surv_type = "OS") {
  if (identical(modal, "Clinical")) {
    return(any(vapply(vars, has_usable_clinical_variable, logical(1), cancers = cancers)))
  }
  if (identical(modal, "Survival")) {
    return(has_usable_survival_data(cancers, surv_type = surv_type))
  }
  if (modal %in% c("RNAseq", "Protein", "logCNA", "Methylation", "Mutation")) {
    return(any(vapply(vars, has_usable_gene_modal, logical(1), modal = modal, cancers = cancers)))
  }
  if (identical(modal, "Phospho")) {
    return(FALSE)
  }
  FALSE
}

get_actual_data_skip_reason <- function(task) {
  args <- task$args
  surv_type <- args$surv_type %||% "OS"

  if (!has_usable_modal_set(args$var1, args$var1_modal, args$var1_cancers, surv_type = surv_type)) {
    return(paste0(
      "Skipped because local SLCPTAC bulk data has no usable ",
      args$var1_modal,
      " values for ",
      paste(args$var1, collapse = ", "),
      " in ",
      paste(args$var1_cancers, collapse = ", "),
      "."
    ))
  }

  if (identical(task$function_name, "cptac_correlation")) {
    if (!has_usable_modal_set(args$var2, args$var2_modal, args$var2_cancers, surv_type = surv_type)) {
      return(paste0(
        "Skipped because local SLCPTAC bulk data has no usable ",
        args$var2_modal,
        " values for ",
        paste(args$var2, collapse = ", "),
        " in ",
        paste(args$var2_cancers, collapse = ", "),
        "."
      ))
    }
  }

  if (identical(task$function_name, "cptac_enrichment")) {
    if (!has_usable_genome_modal(args$genome_modal, args$var1_cancers)) {
      return(paste0(
        "Skipped because local SLCPTAC bulk data has no usable genome-wide ",
        args$genome_modal,
        " matrix for ",
        paste(args$var1_cancers, collapse = ", "),
        "."
      ))
    }
  }

  if (identical(task$function_name, "cptac_survival")) {
    if (!has_usable_survival_data(args$var1_cancers, surv_type = surv_type)) {
      return(paste0(
        "Skipped because local SLCPTAC bulk data has no usable ",
        surv_type,
        " survival data for ",
        paste(args$var1_cancers, collapse = ", "),
        "."
      ))
    }
  }

  ""
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
    Smoke_Test_Mode = SMOKE_TEST_MODE,
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
	      "RUN_AUXILIARY_AVAILABLE_SCENARIOS", "SMOKE_TEST_MODE",
	      "SMOKE_TEST_SAMPLES_PER_CANCER", "STOP_IF_BULK_DATA_MISSING",
	      "ALLOW_FORK_PARALLEL", "DISABLE_FORK_PARALLEL",
	      "PARALLEL_BACKEND", "PARALLEL_WORKERS"
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
      as.character(RUN_AUXILIARY_AVAILABLE_SCENARIOS),
      as.character(SMOKE_TEST_MODE),
	      as.character(SMOKE_TEST_SAMPLES_PER_CANCER),
	      as.character(STOP_IF_BULK_DATA_MISSING),
	      as.character(ALLOW_FORK_PARALLEL),
	      as.character(DISABLE_FORK_PARALLEL),
	      as.character(PARALLEL_BACKEND),
	      as.character(PARALLEL_WORKERS)
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

format_outer_title <- function(title, width = 76L) {
  if (is.null(title) || !nzchar(title)) {
    return("")
  }
  paste(strwrap(title, width = width), collapse = "\n")
}

outer_title_line_count <- function(title) {
  formatted <- format_outer_title(title)
  if (!nzchar(formatted)) {
    return(0L)
  }
  length(strsplit(formatted, "\n", fixed = TRUE)[[1]])
}

outer_title_top_fraction <- function(title) {
  if (!ADD_SCRIPT_TITLES || is.null(title) || !nzchar(title)) {
    return(0)
  }
  line_count <- outer_title_line_count(title)
  min(0.25, 0.115 + 0.045 * max(0L, line_count - 1L))
}

draw_outer_title <- function(title) {
  if (!ADD_SCRIPT_TITLES || is.null(title) || !nzchar(title)) {
    return(invisible(FALSE))
  }
  title_text <- format_outer_title(title)
  top_fraction <- outer_title_top_fraction(title)
  grid::grid.text(
    label = title_text,
    x = grid::unit(0.5, "npc"),
    y = grid::unit(1 - top_fraction / 2, "npc"),
    just = "center",
    gp = grid::gpar(fontface = "bold", fontsize = 12, lineheight = 0.95)
  )
  invisible(TRUE)
}

draw_in_title_viewport <- function(title, draw_fun) {
  if (!ADD_SCRIPT_TITLES || is.null(title) || !nzchar(title)) {
    draw_fun()
    return(invisible(NULL))
  }

  grid::grid.newpage()
  top_fraction <- outer_title_top_fraction(title)
  body_height <- max(0.70, 1 - top_fraction - 0.025)
  grid::pushViewport(grid::viewport(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.0, "npc"),
    width = grid::unit(1, "npc"),
    height = grid::unit(body_height, "npc"),
    just = c("center", "bottom")
  ))
  draw_fun()
  try(grid::popViewport(), silent = TRUE)
  draw_outer_title(title)
  invisible(NULL)
}

render_slcptac_plot <- function(plot, title = NULL) {
  if (is.null(plot)) {
    draw_outer_title(title)
    return(invisible(NULL))
  }

  if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) {
    draw_in_title_viewport(title, function() {
      print(plot, newpage = FALSE)
    })
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

  actual_data_reason <- get_actual_data_skip_reason(task)
  if (nzchar(actual_data_reason)) {
    return(actual_data_reason)
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
  status %in% c("success", "cached")
}

is_error_status <- function(status) {
  !(status %in% c("success", "cached", "skipped"))
}

make_runtime_summary <- function(task_summary, start_time) {
  total_runtime_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  task_status <- if (nrow(task_summary) > 0L) task_summary$Status else character(0)

  data.frame(
    Section = c("slcptac_tasks", "total"),
    Tasks = c(nrow(task_summary), nrow(task_summary)),
    Completed = c(sum(is_success_status(task_status)), sum(is_success_status(task_status))),
    Skipped = c(sum(task_status == "skipped"), sum(task_status == "skipped")),
    Failed = c(sum(is_error_status(task_status)), sum(is_error_status(task_status))),
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

  unlink(RESULT_ROOT, recursive = TRUE, force = TRUE)
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

if (SMOKE_TEST_MODE) {
  if (QUICKANALYSIS_VERBOSE) {
    cat("\nSLCPTAC smoke test mode is enabled. Generating temporary CPTAC-like bulk data...\n")
  }
  create_slcptac_smoke_bulk_data(SLCPTAC_BULK_DATA_ROOT)
}

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

if (QUICKANALYSIS_VERBOSE) {
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
}

prepare_linkedomics_bulk_data(slcptac_tasks)

if (QUICKANALYSIS_VERBOSE) {
  cat("\nWriting SLCPTAC analysis catalog and data audit...\n")
}
write_analysis_catalog(selected_analyses)
run_data_summary(slcptac_tasks)

if (QUICKANALYSIS_VERBOSE) {
  cat("\nSLCPTAC plotting/statistical tasks: ", length(slcptac_tasks), "\n", sep = "")
}

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
  failed_tasks <- task_summary[is_error_status(task_summary$Status) | task_summary$Status == "skipped", , drop = FALSE]
  write_table(failed_tasks, TABLE_ROOT, "000_slcptac_failed_tasks")
}

runtime_summary <- make_runtime_summary(task_summary = task_summary, start_time = SCRIPT_START_TIME)
write_table(runtime_summary, TABLE_ROOT, "000_slcptac_runtime_summary")

if (QUICKANALYSIS_VERBOSE) {
  cat("\nSLCPTAC quick analysis finished.\n")
  if (nrow(task_summary) > 0L) {
    cat("Completed tasks: ", sum(is_success_status(task_summary$Status)), "\n", sep = "")
    cat("Skipped tasks: ", sum(task_summary$Status == "skipped"), "\n", sep = "")
    cat("Failed tasks: ", sum(is_error_status(task_summary$Status)), "\n", sep = "")
    cat("Task summary: ", file.path(TABLE_ROOT, "000_slcptac_task_summary.csv"), "\n", sep = "")
    cat("Failed task summary: ", file.path(TABLE_ROOT, "000_slcptac_failed_tasks.csv"), "\n", sep = "")
  }
  cat("Analysis catalog: ", file.path(TABLE_ROOT, "000_slcptac_analysis_catalog.csv"), "\n", sep = "")
  cat("Task design: ", file.path(TABLE_ROOT, "000_slcptac_task_design.csv"), "\n", sep = "")
  cat("Required file audit: ", file.path(TABLE_ROOT, "000_slcptac_required_files.csv"), "\n", sep = "")
  cat("Runtime summary: ", file.path(TABLE_ROOT, "000_slcptac_runtime_summary.csv"), "\n", sep = "")
}

cat("\n02 SLCPTAC quick analysis finished: ", RESULT_ROOT, "\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
