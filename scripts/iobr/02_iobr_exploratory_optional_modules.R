# IOBR探索性扩展分析
#
# 本脚本放置和当前实验目的相关但条件依赖更强、或不如01号脚本直接的IOBR模块：
# 1. TCGA/GTEx：基于01号脚本输出的TME/signature表，围绕目标基因High/Low分组做ROC、
#    聚类组成等探索；TCGA额外尝试生存分析；
# 2. GSE114012：不做目标基因High/Low扩展，仅保留LRC vs BULK设计相关的IOBR结果；
# 3. 基于GSE114012现有NGS差异分析表调用IOBR::sig_gsea，作为LRC vs BULK差异分析的旁证；
# 4. 离群样本检测作为QC模块保存到QC/all_samples层级，不参与功能分组解释。
# 5. 每个脚本只清理自己生成的OUTPUT_ROOT；01号core与02号exploratory互不清理。
#
# 运行前建议先运行：
# Rscript scripts/iobr/01_iobr_core_target_gene_and_gse114012.R


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


# 1. 分析设计区 ----------------------------------------------------------------

TARGET_GENES <- iobr_parse_env_vector("IOBR_TARGET_GENES", c("ATF3"))
CORE_OUTPUT_ROOT <- file.path(PROJECT_ROOT, "results", "iobr", "core")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "results", "iobr", "exploratory")
TEMP_ROOT <- file.path(PROJECT_ROOT, "temporary", "iobr", "exploratory")
DATA_ROOT <- file.path(PROJECT_ROOT, "data", "iobr")
IOBR_CACHE_DIR <- file.path(DATA_ROOT, "IOBR_cache")

RUN_SURVIVAL <- iobr_parse_env_logical("IOBR_EXPLORE_RUN_SURVIVAL", TRUE)
RUN_ROC <- iobr_parse_env_logical("IOBR_EXPLORE_RUN_ROC", TRUE)
RUN_CLUSTER <- iobr_parse_env_logical("IOBR_EXPLORE_RUN_CLUSTER", TRUE)
RUN_OUTLIER <- iobr_parse_env_logical("IOBR_EXPLORE_RUN_OUTLIER", TRUE)
RUN_SIG_GSEA <- iobr_parse_env_logical("IOBR_EXPLORE_RUN_SIG_GSEA", TRUE)

TOP_FEATURE_N <- iobr_parse_env_integer("IOBR_EXPLORE_TOP_FEATURE_N", 12L)
GSE114012_DEG_ROOT <- file.path(PROJECT_ROOT, "results", "ngs", "GSE114012", "tables")
SIG_GSEA_CATEGORIES <- iobr_parse_env_vector("IOBR_EXPLORE_SIG_GSEA_CATEGORIES", c("H", "C2", "C5"))
SIG_GSEA_MAX_DESIGNS <- iobr_parse_env_integer("IOBR_EXPLORE_SIG_GSEA_MAX_DESIGNS", Inf)

CLEAR_PREVIOUS_RUN_OUTPUTS <- iobr_parse_env_logical("IOBR_EXPLORE_CLEAR_PREVIOUS", TRUE)
IOBR_VERBOSE <- iobr_parse_env_logical("IOBR_EXPLORE_VERBOSE", FALSE)
DISABLE_FORK_PARALLEL <- iobr_parse_env_logical("IOBR_EXPLORE_DISABLE_FORK", interactive())
PARALLEL_BACKEND <- Sys.getenv(
  "IOBR_EXPLORE_PARALLEL_BACKEND",
  unset = Sys.getenv("PARALLEL_RUNTIME_BACKEND", unset = "auto")
)
MAX_PARALLEL_WORKERS <- iobr_parse_env_integer(
  "IOBR_EXPLORE_PARALLEL_WORKERS",
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
  "IOBR", "ggplot2", "dplyr", "survival", "pROC", "qs2", "Cairo"
))

iobr_prepare_output_tree(
  output_root = OUTPUT_ROOT,
  temporary_root = TEMP_ROOT,
  clear_previous = CLEAR_PREVIOUS_RUN_OUTPUTS
)

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
})

SCRIPT_START_TIME <- start_runtime_timer()
PARALLEL_WORKERS <- if (is.na(MAX_PARALLEL_WORKERS)) {
  get_available_worker_count()
} else {
  max(1L, MAX_PARALLEL_WORKERS)
}

explore_config <- data.frame(
  Key = c(
    "Target_genes", "Run_survival", "Run_roc", "Run_cluster",
    "Run_outlier", "Run_sig_gsea", "Top_feature_n",
    "Sig_gsea_categories", "Sig_gsea_max_designs",
    "Clear_previous_run_outputs", "Parallel_backend", "Parallel_workers",
    "Core_output_root", "Output_root", "Temporary_root"
  ),
  Value = c(
    paste(TARGET_GENES, collapse = ", "),
    RUN_SURVIVAL,
    RUN_ROC,
    RUN_CLUSTER,
    RUN_OUTLIER,
    RUN_SIG_GSEA,
    TOP_FEATURE_N,
    paste(SIG_GSEA_CATEGORIES, collapse = ", "),
    SIG_GSEA_MAX_DESIGNS,
    CLEAR_PREVIOUS_RUN_OUTPUTS,
    PARALLEL_BACKEND,
    PARALLEL_WORKERS,
    CORE_OUTPUT_ROOT,
    OUTPUT_ROOT,
    TEMP_ROOT
  ),
  stringsAsFactors = FALSE
)
iobr_write_module_csv(explore_config, OUTPUT_ROOT, "run_summary", "000_iobr_exploratory_run_configuration")


# 2. 小工具 --------------------------------------------------------------------

copy_iobr_generated_outputs <- function(from_dir, module, file_prefix) {
  if (!dir.exists(from_dir)) {
    return(invisible(character(0)))
  }

  files <- list.files(from_dir, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0L) {
    return(invisible(character(0)))
  }

  copied <- character(0)
  table_files <- files[tolower(tools::file_ext(files)) %in% c("csv", "txt", "tsv")]
  for (file in table_files) {
    target <- file.path(
      iobr_table_dir(OUTPUT_ROOT, module, file_prefix),
      paste0(iobr_sanitize(file_prefix), "_", iobr_sanitize(basename(file)))
    )
    file.copy(file, target, overwrite = TRUE)
    copied <- c(copied, target)
  }

  pdf_files <- files[tolower(tools::file_ext(files)) == "pdf"]
  for (file in pdf_files) {
    target <- file.path(
      iobr_plot_dir(OUTPUT_ROOT, module, "pdf", file_prefix),
      paste0(iobr_sanitize(file_prefix), "_", iobr_sanitize(basename(file)))
    )
    file.copy(file, target, overwrite = TRUE)
    copied <- c(copied, target)

    if (requireNamespace("pdftools", quietly = TRUE)) {
      png_target <- file.path(
        iobr_plot_dir(OUTPUT_ROOT, module, "png", file_prefix),
        paste0(tools::file_path_sans_ext(basename(target)), ".png")
      )
      try(pdftools::pdf_convert(target, format = "png", dpi = 300, filenames = png_target), silent = TRUE)
    }
  }

  png_files <- files[tolower(tools::file_ext(files)) == "png"]
  for (file in png_files) {
    target <- file.path(
      iobr_plot_dir(OUTPUT_ROOT, module, "png", file_prefix),
      paste0(iobr_sanitize(file_prefix), "_", iobr_sanitize(basename(file)))
    )
    file.copy(file, target, overwrite = TRUE)
    copied <- c(copied, target)
  }

  invisible(copied)
}

read_core_feature_files <- function() {
  files <- list.files(
    CORE_OUTPUT_ROOT,
    pattern = "_scores[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  files <- files[!grepl("/tables/(run_summary|capability_catalog)/", files)]
  out <- data.frame(
    Feature_File = files,
    Module_Source = ifelse(
      grepl("/signature_score/", files),
      "signature_score",
      ifelse(grepl("/deconvolution/", files), "deconvolution", basename(dirname(files)))
    ),
    File_Stem = tools::file_path_sans_ext(basename(files)),
    stringsAsFactors = FALSE
  )
  if (nrow(out) == 0L) {
    return(out)
  }

  contexts <- lapply(seq_len(nrow(out)), function(i) {
    iobr_parse_output_context(out$File_Stem[i], out$Module_Source[i])
  })
  out$Dataset_ID <- vapply(contexts, `[[`, character(1), "dataset_id")
  out$Dataset_Family <- vapply(contexts, `[[`, character(1), "dataset_family")
  out$Target_Block <- vapply(contexts, `[[`, character(1), "target_block")
  out$Analysis_Design <- vapply(contexts, `[[`, character(1), "design")
  out$IOBR_Method <- vapply(contexts, `[[`, character(1), "method")
  out
}

select_numeric_features <- function(dat, top_n = TOP_FEATURE_N) {
  features <- iobr_feature_columns(dat)
  metadata_patterns <- c(
    "^days_to_", "^age_", "^OS_", "^DSS_", "^DFI_", "^PFI_",
    "^vital_status$", "^gender$", "^sex$", "^stage", "^grade",
    "_clinical$", "^Target_", "^Sample_", "^Dataset_", "^ProjectID$",
    "^time$", "^status$"
  )
  metadata_regex <- paste(metadata_patterns, collapse = "|")
  features <- features[!grepl(metadata_regex, features, ignore.case = TRUE)]
  features <- features[vapply(features, function(feature) {
    values <- suppressWarnings(as.numeric(dat[[feature]]))
    sum(is.finite(values)) >= 6L && stats::sd(values, na.rm = TRUE) > 0
  }, logical(1))]
  head(features, min(length(features), top_n))
}

extract_sig_gsea_table <- function(x) {
  if (is.data.frame(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (is.list(x)) {
    for (item in x) {
      found <- extract_sig_gsea_table(item)
      if (is.data.frame(found) && nrow(found) > 0L) {
        return(found)
      }
    }
  }
  data.frame()
}

run_survival_task <- function(task) {
  module <- "survival"
  dat <- read_report_csv(task$Feature_File)
  if (!all(c("time", "status") %in% colnames(dat))) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "No time/status columns."))
  }
  dat$time <- suppressWarnings(as.numeric(dat$time))
  dat$status <- suppressWarnings(as.numeric(dat$status))
  dat <- dat[is.finite(dat$time) & !is.na(dat$status) & dat$time > 0, , drop = FALSE]
  if (nrow(dat) < 20L || length(unique(dat$status)) < 2L) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "Too few survival events."))
  }

  features <- select_numeric_features(dat)
  if (length(features) == 0L) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "No numeric features."))
  }

  res <- suppressWarnings(suppressMessages(IOBR::batch_surv(
    pdata = dat,
    variable = features,
    time = "time",
    status = "status",
    best_cutoff = FALSE
  )))
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  out_file <- iobr_write_module_csv(res, OUTPUT_ROOT, module, paste0(task$File_Stem, "_batch_survival"))

  p <- suppressWarnings(suppressMessages(try(
    IOBR::sig_forest(
      res,
      signature = intersect(c("ID", "signature", "Variable", "variable"), colnames(res))[1]
    ),
    silent = TRUE
  )))
  if (!inherits(p, "try-error") && inherits(p, "ggplot")) {
    iobr_save_module_plot(
      iobr_append_title(p, paste0(task$File_Stem, " survival forest")),
      OUTPUT_ROOT,
      module,
      paste0(task$File_Stem, "_survival_forest"),
      plot_type = "forest",
      n_items = nrow(res),
      max_label_chars = max(nchar(apply(res, 1, paste, collapse = " ")), na.rm = TRUE),
      title = paste0(task$File_Stem, " survival forest")
    )
  }

  iobr_status_row(task$Dataset_ID, module, task$File_Stem, "success", output_file = out_file, n_rows = nrow(res), n_cols = ncol(res))
}

run_roc_task <- function(task) {
  module <- "roc"
  dat <- read_report_csv(task$Feature_File)
  if (!"Target_Group" %in% colnames(dat)) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "No Target_Group column."))
  }
  if (length(unique(dat$Target_Group[!is.na(dat$Target_Group)])) != 2L) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "Target_Group is not binary."))
  }

  features <- select_numeric_features(dat)
  if (length(features) == 0L || !requireNamespace("pROC", quietly = TRUE)) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "No numeric features or pROC missing."))
  }

  roc_rows <- lapply(features, function(feature) {
    ok <- is.finite(dat[[feature]]) & !is.na(dat$Target_Group)
    if (sum(ok) < 6L) {
      return(NULL)
    }
    roc_obj <- suppressWarnings(pROC::roc(dat$Target_Group[ok], dat[[feature]][ok], quiet = TRUE))
    data.frame(
      Feature = feature,
      AUC = as.numeric(pROC::auc(roc_obj)),
      Direction = roc_obj$direction,
      stringsAsFactors = FALSE
    )
  })
  roc_table <- iobr_bind_rows(roc_rows)
  if (nrow(roc_table) == 0L) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "No valid ROC results."))
  }
  roc_table <- roc_table[order(roc_table$AUC, decreasing = TRUE), , drop = FALSE]
  out_file <- iobr_write_module_csv(roc_table, OUTPUT_ROOT, module, paste0(task$File_Stem, "_target_group_roc"))

  p <- iobr_make_top_barplot(
    roc_table,
    feature_col = "Feature",
    value_col = "AUC",
    title = paste0(task$File_Stem, " Target_Group ROC"),
    top_n = min(TOP_FEATURE_N, nrow(roc_table))
  )
  iobr_save_module_plot(
    p,
    OUTPUT_ROOT,
    module,
    paste0(task$File_Stem, "_target_group_roc_auc"),
    plot_type = "bar",
    n_items = min(TOP_FEATURE_N, nrow(roc_table)),
    max_label_chars = max(nchar(as.character(roc_table$Feature)), na.rm = TRUE),
    title = paste0(task$File_Stem, " Target_Group ROC")
  )

  iobr_status_row(task$Dataset_ID, module, task$File_Stem, "success", output_file = out_file, n_rows = nrow(roc_table), n_cols = ncol(roc_table))
}

run_cluster_task <- function(task) {
  module <- "cluster"
  dat <- read_report_csv(task$Feature_File)
  features <- select_numeric_features(dat)
  if (length(features) < 3L || !requireNamespace("NbClust", quietly = TRUE)) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "Need >=3 features and NbClust."))
  }

  clustered <- suppressWarnings(suppressMessages(IOBR::tme_cluster(
    input = dat[, c("ID", features), drop = FALSE],
    features = features,
    id = "ID",
    scale = TRUE,
    method = "kmeans",
    min_nc = 2,
    max.nc = 6
  )))
  clustered <- as.data.frame(clustered, stringsAsFactors = FALSE)
  out_file <- iobr_write_module_csv(clustered, OUTPUT_ROOT, module, paste0(task$File_Stem, "_tme_cluster"))

  if ("cluster" %in% colnames(clustered)) {
    plot_data <- merge(dat[, intersect(c("ID", "Target_Group"), colnames(dat)), drop = FALSE], clustered[, c("ID", "cluster"), drop = FALSE], by = "ID")
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cluster, fill = Target_Group)) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(x = "Cluster", y = "Proportion", title = paste0(task$File_Stem, " cluster composition")) +
      iobr_common_ggplot_theme(
        base_size = BASE_FONT_SIZE,
        axis_text_size = 10,
        axis_title_size = 12.5,
        title_size = 13
      )
    iobr_save_module_plot(
      p,
      OUTPUT_ROOT,
      module,
      paste0(task$File_Stem, "_cluster_composition"),
      plot_type = "composition",
      n_groups = length(unique(plot_data$cluster)),
      max_label_chars = max(nchar(c(as.character(plot_data$cluster), as.character(plot_data$Target_Group))), na.rm = TRUE),
      title = paste0(task$File_Stem, " cluster composition")
    )
  }

  iobr_status_row(task$Dataset_ID, module, task$File_Stem, "success", output_file = out_file, n_rows = nrow(clustered), n_cols = ncol(clustered))
}

run_outlier_task <- function(task) {
  module <- "outlier"
  input <- iobr_cache_read(task$Input_File)
  res <- suppressWarnings(suppressMessages(try(IOBR::find_outlier_samples(
    eset = input$score_expr,
    project = task$Dataset_ID,
    plot_hculst = FALSE,
    show_plot = FALSE,
    save = FALSE
  ), silent = TRUE)))
  if (inherits(res, "try-error") || is.null(res)) {
    return(iobr_status_row(task$Dataset_ID, module, task$Dataset_ID, "failed", as.character(res)))
  }
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  out_file <- iobr_write_module_csv(res, OUTPUT_ROOT, module, paste0(task$Dataset_ID, "_outlier_samples"))
  iobr_status_row(task$Dataset_ID, module, task$Dataset_ID, "success", output_file = out_file, n_rows = nrow(res), n_cols = ncol(res))
}

run_sig_gsea_task <- function(task) {
  module <- "sig_gsea"
  deg <- iobr_read_deg_for_sig_gsea(task$DEG_File)
  if (nrow(deg) < 20L) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "skipped", "Too few DEG rows."))
  }

  temp_dir <- file.path(TEMP_ROOT, module, iobr_sanitize(task$File_Stem))
  unlink(temp_dir, recursive = TRUE, force = TRUE)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  res <- suppressWarnings(suppressMessages(try(IOBR::sig_gsea(
    deg = deg,
    path = temp_dir,
    gene_symbol = "symbol",
    logfc = "log2FoldChange",
    org = "hsa",
    msigdb = TRUE,
    category = task$Category,
    project = task$File_Stem,
    show_plot = FALSE,
    print_bar = FALSE,
    plot_single_sig = FALSE,
    fig.type = "pdf",
    verbose = FALSE,
    seed = 123
  ), silent = TRUE)))
  if (inherits(res, "try-error")) {
    return(iobr_status_row(task$Dataset_ID, module, task$File_Stem, "failed", as.character(res)))
  }

  sig_gsea_table <- extract_sig_gsea_table(res)
  if (nrow(sig_gsea_table) > 0L) {
    out_file <- iobr_write_module_csv(sig_gsea_table, OUTPUT_ROOT, module, paste0(task$File_Stem, "_", task$Category, "_sig_gsea_result"))

    term_col <- intersect(c("Description", "ID", "pathway", "term", "Term"), colnames(sig_gsea_table))[1]
    value_col <- intersect(c("NES", "enrichmentScore", "score"), colnames(sig_gsea_table))[1]
    p_col <- intersect(c("p.adjust", "pvalue", "p.value", "qvalue"), colnames(sig_gsea_table))[1]
    if (!is.na(term_col)) {
      plot_table <- sig_gsea_table
      if (!is.na(value_col)) {
        plot_table$Plot_Value <- suppressWarnings(as.numeric(plot_table[[value_col]]))
        y_label <- value_col
      } else if (!is.na(p_col)) {
        plot_table$Plot_Value <- -log10(pmax(suppressWarnings(as.numeric(plot_table[[p_col]])), .Machine$double.xmin))
        y_label <- paste0("-log10(", p_col, ")")
      } else {
        plot_table$Plot_Value <- seq_len(nrow(plot_table))
        y_label <- "Rank"
      }
      plot_table <- plot_table[is.finite(plot_table$Plot_Value), , drop = FALSE]
      plot_table <- head(plot_table[order(abs(plot_table$Plot_Value), decreasing = TRUE), , drop = FALSE], TOP_FEATURE_N)
      if (nrow(plot_table) > 0L) {
        if (!is.na(p_col)) {
          plot_table$P_Label <- iobr_make_p_label(plot_table[[p_col]])
        } else {
          plot_table$P_Label <- ""
        }
        p_gsea <- iobr_make_top_barplot(
          plot_table,
          feature_col = term_col,
          value_col = "Plot_Value",
          title = paste(task$File_Stem, task$Category, "IOBR sig_gsea"),
          xlab = NULL,
          top_n = nrow(plot_table),
          label_col = "P_Label",
          ylab = y_label
        )
        iobr_save_module_plot(
          p_gsea,
          OUTPUT_ROOT,
          module,
          paste0(task$File_Stem, "_", task$Category, "_sig_gsea_top", TOP_FEATURE_N),
          plot_type = "bar",
          n_items = nrow(plot_table),
          max_label_chars = max(nchar(as.character(plot_table[[term_col]])), na.rm = TRUE),
          title = paste(task$File_Stem, task$Category, "IOBR sig_gsea")
        )
      }
    }
  } else {
    out_file <- iobr_write_module_csv(
      data.frame(Result_Class = class(res), stringsAsFactors = FALSE),
      OUTPUT_ROOT,
      module,
      paste0(task$File_Stem, "_", task$Category, "_sig_gsea_result_class")
    )
  }
  copy_iobr_generated_outputs(temp_dir, module, paste0(task$File_Stem, "_", task$Category))
  iobr_status_row(task$Dataset_ID, module, task$File_Stem, "success", output_file = out_file)
}


# 3. 构建任务表 ----------------------------------------------------------------

task_list <- list()

core_feature_files <- read_core_feature_files()
if (nrow(core_feature_files) > 0L) {
  core_feature_files$Task_ID_Base <- seq_len(nrow(core_feature_files))

  for (i in seq_len(nrow(core_feature_files))) {
    is_target_gene_table <- grepl("^GENE_", core_feature_files$Target_Block[i])
    is_tcga <- identical(core_feature_files$Dataset_Family[i], "TCGA")
    is_gtex <- identical(core_feature_files$Dataset_Family[i], "GTEx")

    if (RUN_SURVIVAL && is_tcga && is_target_gene_table) {
      task_list[[length(task_list) + 1L]] <- data.frame(Task_Type = "survival", core_feature_files[i, ], stringsAsFactors = FALSE)
    }
    if (RUN_ROC && is_target_gene_table && (is_tcga || is_gtex)) {
      task_list[[length(task_list) + 1L]] <- data.frame(Task_Type = "roc", core_feature_files[i, ], stringsAsFactors = FALSE)
    }
    if (RUN_CLUSTER && is_target_gene_table && (is_tcga || is_gtex)) {
      task_list[[length(task_list) + 1L]] <- data.frame(Task_Type = "cluster", core_feature_files[i, ], stringsAsFactors = FALSE)
    }
  }
}

manifest_file <- file.path(CORE_OUTPUT_ROOT, "tables", "run_summary", "002_prepared_input_manifest.csv")
if (RUN_OUTLIER && file.exists(manifest_file)) {
  manifest <- read_report_csv(manifest_file)
  for (i in seq_len(nrow(manifest))) {
    task_list[[length(task_list) + 1L]] <- data.frame(
      Task_Type = "outlier",
      Dataset_ID = manifest$Dataset_ID[i],
      Input_File = manifest$Cache_File[i],
      Feature_File = NA_character_,
      Module_Source = "prepared_input",
      File_Stem = manifest$Dataset_ID[i],
      stringsAsFactors = FALSE
    )
  }
}

if (RUN_SIG_GSEA && dir.exists(GSE114012_DEG_ROOT)) {
  deg_files <- list.files(GSE114012_DEG_ROOT, pattern = "all_genes\\.csv$", recursive = TRUE, full.names = TRUE)
  if (is.finite(SIG_GSEA_MAX_DESIGNS)) {
    deg_files <- head(deg_files, SIG_GSEA_MAX_DESIGNS)
  }
  for (deg_file in deg_files) {
    analysis_name <- basename(dirname(dirname(deg_file)))
    for (category in SIG_GSEA_CATEGORIES) {
      task_list[[length(task_list) + 1L]] <- data.frame(
        Task_Type = "sig_gsea",
        Dataset_ID = "GSE114012",
        Feature_File = NA_character_,
        Module_Source = "GSE114012_DEG",
        File_Stem = paste("GSE114012", analysis_name, category, sep = "_"),
        DEG_File = deg_file,
        Category = category,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(task_list) == 0L) {
  iobr_write_module_csv(
    data.frame(Message = "No exploratory IOBR tasks were prepared.", stringsAsFactors = FALSE),
    OUTPUT_ROOT,
    "run_summary",
    "001_no_tasks"
  )
  cat("\n02 IOBR exploratory analysis finished with no tasks: ", OUTPUT_ROOT, "\n", sep = "")
  quit(save = "no", status = 0)
}

explore_task_table <- dplyr::bind_rows(task_list)
explore_task_table$Task_ID <- seq_len(nrow(explore_task_table))
iobr_write_module_csv(explore_task_table, OUTPUT_ROOT, "run_summary", "001_iobr_exploratory_task_table")


# 4. 并行运行探索任务 ----------------------------------------------------------

run_one_iobr_exploratory_task <- function(task_id) {
  task <- explore_task_table[task_id, , drop = FALSE]
  tryCatch(suppressWarnings(suppressMessages({
    if (identical(task$Task_Type, "survival")) {
      return(run_survival_task(task))
    }
    if (identical(task$Task_Type, "roc")) {
      return(run_roc_task(task))
    }
    if (identical(task$Task_Type, "cluster")) {
      return(run_cluster_task(task))
    }
    if (identical(task$Task_Type, "outlier")) {
      return(run_outlier_task(task))
    }
    if (identical(task$Task_Type, "sig_gsea")) {
      return(run_sig_gsea_task(task))
    }
    iobr_status_row(task$Dataset_ID, task$Task_Type, task$File_Stem, "failed", paste("Unknown task type:", task$Task_Type))
  })), error = function(error) {
    iobr_status_row(task$Dataset_ID, task$Task_Type, task$File_Stem, "failed", conditionMessage(error))
  })
}

parallel_strategy <- setup_parallel_strategy(
  total_tasks = nrow(explore_task_table),
  max_workers = PARALLEL_WORKERS,
  inner_label = "IOBR exploratory inner workers",
  nested_label = "IOBR nested workers"
)

task_results <- run_parallel_tasks_with_progress(
  task_ids = explore_task_table$Task_ID,
  task_function = run_one_iobr_exploratory_task,
  workers = parallel_strategy$task_workers,
  progress_label = "IOBR explore"
)

status_table <- iobr_bind_rows(task_results)
iobr_write_module_csv(status_table, OUTPUT_ROOT, "run_summary", "002_iobr_exploratory_task_status")
failed <- status_table[status_table$Status %in% c("failed"), , drop = FALSE]
if (nrow(failed) > 0L) {
  iobr_write_module_csv(failed, OUTPUT_ROOT, "run_summary", "003_iobr_exploratory_failed_tasks")
}

runtime_seconds <- print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
runtime_table <- data.frame(
  Total_Runtime_Seconds = runtime_seconds,
  Total_Runtime = format_runtime_seconds(runtime_seconds),
  Finished_At = as.character(Sys.time()),
  Output_Root = OUTPUT_ROOT,
  stringsAsFactors = FALSE
)
iobr_write_module_csv(runtime_table, OUTPUT_ROOT, "run_summary", "004_iobr_exploratory_runtime")

cat("\n02 IOBR exploratory analysis finished: ", OUTPUT_ROOT, "\n", sep = "")
