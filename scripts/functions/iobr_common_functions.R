# IOBR分析公共函数
#
# 本文件服务于scripts/iobr目录下的IOBR快速分析脚本：
# - 统一项目路径、IOBR缓存、表格和图片保存逻辑；
# - 统一读取GSE114012 SummarizedExperiment与差异分析结果；
# - 统一调用parallel_runtime_functions.R的多进程进度条；
# - 每个IOBR模块独立记录success/skipped/failed，避免单个方法失败中断全流程。


# 0. 环境变量与路径 ------------------------------------------------------------

iobr_parse_env_vector <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  parts <- trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
  parts[nzchar(parts)]
}

iobr_parse_env_logical <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(default)
  }

  value %in% c("1", "true", "t", "yes", "y")
}

iobr_parse_env_integer <- function(name, default) {
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

iobr_find_project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(current, "scripts", "functions", "iobr_common_functions.R")
    if (file.exists(marker)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot locate project root from: ", getwd())
    }
    current <- parent
  }
}

iobr_source_project_functions <- function(project_root) {
  source(file.path(project_root, "scripts", "functions", "plotting_common_functions.R"))
  source(file.path(project_root, "scripts", "functions", "result_table_io_functions.R"))
  source(file.path(project_root, "scripts", "functions", "parallel_runtime_functions.R"))
  source(file.path(project_root, "scripts", "functions", "limma_de_functions.R"))
  invisible(TRUE)
}

iobr_sanitize <- function(x, default = "analysis") {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- default
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == ""] <- default
  x
}


# 1. 包安装、能力目录与输出 ----------------------------------------------------

iobr_install_if_needed <- function(auto_install = TRUE) {
  if (requireNamespace("IOBR", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  if (!auto_install) {
    stop(
      "IOBR is not installed. Install with: ",
      "remotes::install_github('IOBR/IOBR', dependencies = TRUE)"
    )
  }

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("IOBR/IOBR", dependencies = TRUE, upgrade = "never")

  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("IOBR installation was attempted but the package is still unavailable.")
  }

  invisible(TRUE)
}

iobr_setup_runtime <- function(
    project_root,
    iobr_cache_dir,
    auto_install = TRUE,
    parallel_backend = "auto",
    quiet_strategy = TRUE) {
  dir.create(iobr_cache_dir, recursive = TRUE, showWarnings = FALSE)
  options(IOBR.cache_dir = iobr_cache_dir)
  options(parallel_runtime_backend = parallel_backend)
  options(parallel_runtime_force_single_line_progress = TRUE)
  options(parallel_runtime_quiet_strategy = quiet_strategy)

  iobr_install_if_needed(auto_install = auto_install)
  suppressPackageStartupMessages(library(IOBR))

  if (exists("set_iobr_cache_dir", envir = asNamespace("IOBR"), mode = "function")) {
    try(IOBR::set_iobr_cache_dir(iobr_cache_dir), silent = TRUE)
  }

  invisible(TRUE)
}

iobr_prepare_output_tree <- function(output_root, temporary_root, clear_previous = TRUE) {
  if (clear_previous) {
    unlink(output_root, recursive = TRUE, force = TRUE)
    unlink(temporary_root, recursive = TRUE, force = TRUE)
  }

  dirs <- c(
    output_root,
    file.path(output_root, "tables"),
    file.path(output_root, "plots"),
    file.path(output_root, "logs"),
    temporary_root
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

iobr_write_csv <- function(dat, file, n_rows = 21, na = "NA") {
  write_csv_with_report_previews(dat, file, n_rows = n_rows, na = na)
}

iobr_save_ggplot <- function(plot, plot_root, module, file_stem, width = 7.2, height = 5.4) {
  pdf_dir <- file.path(plot_root, module, "pdf")
  png_dir <- file.path(plot_root, module, "png")
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

  file_stem <- iobr_sanitize(file_stem)
  pdf_file <- file.path(pdf_dir, paste0(file_stem, ".pdf"))
  save_ggplot_pdf_png(plot, pdf_file = pdf_file, width = width, height = height)
}

iobr_empty_result <- function(module, task, status, message = "", output_file = "") {
  data.frame(
    Module = module,
    Task = task,
    Status = status,
    Message = as.character(message),
    Output_File = as.character(output_file),
    stringsAsFactors = FALSE
  )
}

iobr_safe_task <- function(expr, module, task) {
  tryCatch(
    expr,
    error = function(error) {
      iobr_empty_result(module, task, "failed", conditionMessage(error))
    }
  )
}

iobr_get_exported_functions <- function() {
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    return(character(0))
  }
  sort(getNamespaceExports("IOBR"))
}

iobr_get_tme_methods <- function() {
  if (requireNamespace("IOBR", quietly = TRUE) &&
      exists("tme_deconvolution_methods", envir = asNamespace("IOBR"))) {
    methods <- get("tme_deconvolution_methods", envir = asNamespace("IOBR"))
    return(unique(unname(as.character(methods))))
  }

  c(
    "mcpcounter", "epic", "xcell", "cibersort", "cibersort_abs",
    "ips", "estimate", "timer", "quantiseq", "svr", "lsei"
  )
}

iobr_get_signature_methods <- function() {
  if (requireNamespace("IOBR", quietly = TRUE) &&
      exists("signature_score_calculation_methods", envir = asNamespace("IOBR"))) {
    methods <- get("signature_score_calculation_methods", envir = asNamespace("IOBR"))
    return(unique(unname(as.character(methods))))
  }
  c("pca", "zscore", "ssgsea", "integration")
}

iobr_write_capability_catalog <- function(output_root, scope = "core") {
  exported <- iobr_get_exported_functions()
  core_functions <- c(
    "deconvo_tme", "deconvo_mcpcounter", "deconvo_epic", "deconvo_xcell",
    "deconvo_cibersort", "deconvo_estimate", "deconvo_ips",
    "deconvo_timer", "deconvo_quantiseq", "calculate_sig_score",
    "sigScore", "batch_wilcoxon", "batch_kruskal", "batch_cor",
    "sig_box", "sig_box_batch", "sig_gsea", "iobr_pca", "tme_cluster"
  )
  optional_functions <- c(
    "LR_cal", "sig_roc", "batch_surv", "sig_surv_plot", "PrognosticModel",
    "PrognosticResult", "RegressionResult", "lasso_select", "add_riskscore",
    "find_mutations", "make_mut_matrix", "RemoveBatchEffect", "remove_batcheffect",
    "iobr_deg", "limma.dif", "feature_select", "find_outlier_samples"
  )

  catalog <- data.frame(
    Function = sort(unique(c(core_functions, optional_functions, exported))),
    Exported_By_IOBR = sort(unique(c(core_functions, optional_functions, exported))) %in% exported,
    Planned_Scope = ifelse(
      sort(unique(c(core_functions, optional_functions, exported))) %in% core_functions,
      "core_directly_relevant",
      ifelse(
        sort(unique(c(core_functions, optional_functions, exported))) %in% optional_functions,
        "optional_or_input_dependent",
        "available_not_prioritized"
      )
    ),
    Current_Script = scope,
    stringsAsFactors = FALSE
  )

  iobr_write_csv(
    catalog,
    file.path(output_root, "tables", "capability_catalog", "000_iobr_function_catalog.csv")
  )
}


# 2. 数据读取与整理 ------------------------------------------------------------

iobr_load_gse114012_inputs <- function(
    se_file,
    clinical_file,
    assay_name = "tpm",
    symbol_column = "Symbol",
    biotype_filter = "protein_coding") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package SummarizedExperiment is required.")
  }
  se <- readRDS(se_file)
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(assay_name %in% names(SummarizedExperiment::assays(se)))
  stopifnot("Sample_ID" %in% colnames(clinical))

  expr <- as.matrix(SummarizedExperiment::assay(se, assay_name))
  row_data <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)

  missing_samples <- setdiff(colnames(expr), clinical$Sample_ID)
  if (length(missing_samples) > 0L) {
    stop("Samples missing from clinical table: ", paste(missing_samples, collapse = ", "))
  }
  clinical <- clinical[match(colnames(expr), clinical$Sample_ID), , drop = FALSE]

  expr_symbol <- iobr_collapse_expression_by_symbol(
    expr = expr,
    row_data = row_data,
    symbol_column = symbol_column,
    biotype_filter = biotype_filter
  )

  list(
    se = se,
    clinical = clinical,
    expression = expr_symbol,
    row_data = row_data,
    se_file = se_file,
    clinical_file = clinical_file,
    assay_name = assay_name
  )
}

iobr_collapse_expression_by_symbol <- function(
    expr,
    row_data,
    symbol_column = "Symbol",
    biotype_filter = "protein_coding") {
  stopifnot(nrow(expr) == nrow(row_data))
  stopifnot(symbol_column %in% colnames(row_data))

  keep <- !is.na(row_data[[symbol_column]]) & trimws(row_data[[symbol_column]]) != ""
  if (!is.null(biotype_filter) &&
      "Biotype" %in% colnames(row_data) &&
      length(biotype_filter) > 0L) {
    keep <- keep & row_data$Biotype %in% biotype_filter
  }

  expr <- expr[keep, , drop = FALSE]
  symbols <- trimws(as.character(row_data[[symbol_column]][keep]))
  symbols <- make.unique(symbols, sep = "__dup__")
  base_symbols <- sub("__dup__[0-9]+$", "", symbols)

  storage.mode(expr) <- "double"
  collapsed <- rowsum(expr, group = base_symbols, reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(factor(base_symbols, levels = rownames(collapsed))))
  collapsed <- sweep(collapsed, 1, counts, "/")
  collapsed <- collapsed[order(rownames(collapsed)), , drop = FALSE]
  collapsed
}

iobr_detect_sample_id_column <- function(col_data) {
  candidates <- c("Sample_ID", "sample_id", "ID", "Run", "Title")
  detected <- candidates[candidates %in% colnames(col_data)][1]
  if (is.na(detected)) {
    return(NULL)
  }
  detected
}

iobr_filter_samples <- function(
    col_data,
    sample_id,
    detail_values = NULL,
    detail_column = "group_detail",
    barcode_pattern = NULL) {
  keep <- rep(TRUE, nrow(col_data))
  detail_filter_applied <- FALSE

  if (!is.null(detail_values) &&
      length(detail_values) > 0L &&
      detail_column %in% colnames(col_data)) {
    detail <- tolower(trimws(as.character(col_data[[detail_column]])))
    keep <- keep & detail %in% tolower(detail_values)
    detail_filter_applied <- TRUE
  }

  if (!is.null(barcode_pattern) &&
      nzchar(barcode_pattern) &&
      (!detail_filter_applied || !any(keep)) &&
      length(sample_id) == nrow(col_data)) {
    keep <- grepl(barcode_pattern, sample_id, ignore.case = TRUE)
  }

  keep
}

iobr_prepare_local_se_input <- function(
    dataset_id,
    se_file,
    target_gene,
    project_id = dataset_id,
    assay_name = "tpm",
    deconvolution_assay = assay_name,
    score_assay = assay_name,
    score_log2 = TRUE,
    sample_detail_values = NULL,
    sample_detail_column = "group_detail",
    sample_barcode_pattern = NULL,
    symbol_column = "Symbol",
    biotype_filter = "protein_coding",
    high_group_label = "High",
    low_group_label = "Low",
    tumor = TRUE,
    timer_indication = NULL,
    clinical_file = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package SummarizedExperiment is required.")
  }
  if (!file.exists(se_file)) {
    stop("SE file not found: ", se_file)
  }

  se <- readRDS(se_file)
  stopifnot(inherits(se, "SummarizedExperiment"))
  assay_names <- names(SummarizedExperiment::assays(se))
  if (!deconvolution_assay %in% assay_names) {
    stop("Assay not found for deconvolution: ", deconvolution_assay)
  }
  if (!score_assay %in% assay_names) {
    stop("Assay not found for signature scoring: ", score_assay)
  }

  col_data <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  sample_id_column <- iobr_detect_sample_id_column(col_data)
  sample_id <- if (is.null(sample_id_column)) {
    colnames(se)
  } else {
    as.character(col_data[[sample_id_column]])
  }
  sample_id[is.na(sample_id) | sample_id == ""] <- colnames(se)[is.na(sample_id) | sample_id == ""]
  stopifnot(!any(duplicated(sample_id)))

  keep_sample <- iobr_filter_samples(
    col_data = col_data,
    sample_id = sample_id,
    detail_values = sample_detail_values,
    detail_column = sample_detail_column,
    barcode_pattern = sample_barcode_pattern
  )
  if (!any(keep_sample)) {
    stop("No samples remain after filtering for dataset: ", dataset_id)
  }

  deconvolution_expr <- as.matrix(SummarizedExperiment::assay(se, deconvolution_assay))
  score_expr <- as.matrix(SummarizedExperiment::assay(se, score_assay))
  row_data <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)

  deconvolution_expr <- deconvolution_expr[, keep_sample, drop = FALSE]
  score_expr <- score_expr[, keep_sample, drop = FALSE]
  col_data <- col_data[keep_sample, , drop = FALSE]
  sample_id <- sample_id[keep_sample]
  colnames(deconvolution_expr) <- sample_id
  colnames(score_expr) <- sample_id

  deconvolution_expr <- iobr_collapse_expression_by_symbol(
    expr = deconvolution_expr,
    row_data = row_data,
    symbol_column = symbol_column,
    biotype_filter = biotype_filter
  )
  score_expr <- iobr_collapse_expression_by_symbol(
    expr = score_expr,
    row_data = row_data,
    symbol_column = symbol_column,
    biotype_filter = biotype_filter
  )

  if (score_log2) {
    score_expr <- log2(score_expr + 1)
  }

  target_gene <- trimws(as.character(target_gene)[1])
  if (!target_gene %in% rownames(score_expr)) {
    stop("Target gene not found in ", dataset_id, ": ", target_gene)
  }

  target_expression <- as.numeric(score_expr[target_gene, ])
  target_median <- median(target_expression, na.rm = TRUE)
  target_group <- ifelse(
    target_expression >= target_median,
    high_group_label,
    low_group_label
  )
  target_group <- factor(target_group, levels = c(low_group_label, high_group_label))

  pdata <- cbind(
    data.frame(
      ID = sample_id,
      Sample_ID = sample_id,
      Dataset_ID = dataset_id,
      ProjectID = project_id,
      Target_Gene = target_gene,
      Target_Expression = target_expression,
      Target_Median = target_median,
      Target_Group = as.character(target_group),
      stringsAsFactors = FALSE
    ),
    col_data,
    stringsAsFactors = FALSE
  )
  pdata <- pdata[, !duplicated(colnames(pdata)), drop = FALSE]

  if (!is.null(clinical_file) && file.exists(clinical_file)) {
    clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
    if ("Sample_ID" %in% colnames(clinical)) {
      pdata <- merge(pdata, clinical, by = "Sample_ID", all.x = TRUE, suffixes = c("", "_clinical"))
      pdata <- pdata[match(sample_id, pdata$Sample_ID), , drop = FALSE]
    }
  }

  list(
    dataset_id = dataset_id,
    project_id = project_id,
    se_file = se_file,
    assay_name = assay_name,
    deconvolution_assay = deconvolution_assay,
    score_assay = score_assay,
    score_log2 = score_log2,
    deconvolution_expr = deconvolution_expr,
    score_expr = score_expr,
    pdata = pdata,
    target_gene = target_gene,
    tumor = tumor,
    timer_indication = timer_indication,
    sample_count = ncol(score_expr),
    gene_count = nrow(score_expr)
  )
}

iobr_add_survival_columns <- function(pdata) {
  pdata <- as.data.frame(pdata, stringsAsFactors = FALSE, check.names = FALSE)
  if ("days_to_death" %in% colnames(pdata) || "days_to_last_follow_up" %in% colnames(pdata)) {
    death <- if ("days_to_death" %in% colnames(pdata)) {
      suppressWarnings(as.numeric(pdata$days_to_death))
    } else {
      rep(NA_real_, nrow(pdata))
    }
    follow <- if ("days_to_last_follow_up" %in% colnames(pdata)) {
      suppressWarnings(as.numeric(pdata$days_to_last_follow_up))
    } else {
      rep(NA_real_, nrow(pdata))
    }
    pdata$time <- ifelse(is.na(death), follow, death)
  }

  if ("OS_status" %in% colnames(pdata)) {
    pdata$status <- suppressWarnings(as.numeric(pdata$OS_status))
  } else if ("vital_status" %in% colnames(pdata)) {
    vital <- tolower(trimws(as.character(pdata$vital_status)))
    pdata$status <- ifelse(vital %in% c("dead", "deceased", "1"), 1, 0)
  }

  pdata
}

iobr_bind_rows <- function(items) {
  items <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, items)
  if (length(items) == 0L) {
    return(data.frame())
  }
  dplyr::bind_rows(items)
}

iobr_feature_columns <- function(dat, exclude = character(0)) {
  default_exclude <- c(
    "ID", "Sample_ID", "Dataset_ID", "ProjectID", "Target_Gene",
    "Target_Expression", "Target_Median", "Target_Group",
    "time", "status"
  )
  iobr_numeric_feature_columns(dat, exclude = unique(c(default_exclude, exclude)))
}

iobr_status_row <- function(
    dataset_id,
    module,
    task,
    status,
    message = "",
    output_file = "",
    n_rows = NA_integer_,
    n_cols = NA_integer_) {
  data.frame(
    Dataset_ID = dataset_id,
    Module = module,
    Task = task,
    Status = status,
    Message = as.character(message),
    Rows = n_rows,
    Columns = n_cols,
    Output_File = as.character(output_file),
    stringsAsFactors = FALSE
  )
}

iobr_prepare_signature_collection <- function(selected_signatures = NULL) {
  data(signature_collection, package = "IOBR", envir = environment())
  signatures <- signature_collection
  if (!is.null(selected_signatures) && length(selected_signatures) > 0L) {
    selected_signatures <- intersect(selected_signatures, names(signatures))
    signatures <- signatures[selected_signatures]
  }
  signatures
}

iobr_load_lm22_reference <- function() {
  ref <- IOBR::load_data("lm22")
  as.data.frame(ref, check.names = FALSE)
}

iobr_make_sample_pdata <- function(clinical) {
  pdata <- clinical
  pdata$ID <- pdata$Sample_ID
  pdata
}

iobr_get_analysis_design_table <- function(clinical) {
  get_analysis_designs(clinical)
}

iobr_make_design_group_table <- function(clinical, design_row) {
  design_samples <- prepare_design_samples(
    sample_info = clinical,
    group_column_index = design_row$Column_Index,
    experiment_group = design_row$Experiment_Group
  )
  data.frame(
    ID = design_samples$sample_info$Sample_ID,
    Analysis_Name = design_row$Analysis_Name,
    Group = as.character(design_samples$group_list),
    Control_Group = design_samples$control_group,
    Experiment_Group = design_row$Experiment_Group,
    stringsAsFactors = FALSE
  )
}

iobr_get_deg_file_info_safe <- function(table_root) {
  get_deg_file_info(table_root)
}

iobr_read_deg_for_sig_gsea <- function(file) {
  dat <- read_report_csv(file)
  symbol_col <- intersect(c("Symbol", "symbol", "Gene", "gene"), colnames(dat))[1]
  logfc_col <- intersect(c("logFC", "log2FoldChange", "avg_log2FC"), colnames(dat))[1]
  if (is.na(symbol_col) || is.na(logfc_col)) {
    stop("DEG file lacks Symbol/logFC columns: ", file)
  }
  dat <- dat[!is.na(dat[[symbol_col]]) & dat[[symbol_col]] != "", , drop = FALSE]
  data.frame(
    symbol = as.character(dat[[symbol_col]]),
    log2FoldChange = as.numeric(dat[[logfc_col]]),
    stringsAsFactors = FALSE
  )
}


# 3. 结果矩阵、统计和绘图 ------------------------------------------------------

iobr_normalize_id_column <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  id_candidates <- intersect(c("ID", "Sample", "sample", "Sample_ID", "Mixture"), colnames(dat))
  if (length(id_candidates) == 0L) {
    dat$ID <- rownames(dat)
  } else if (!identical(id_candidates[1], "ID")) {
    colnames(dat)[match(id_candidates[1], colnames(dat))] <- "ID"
  }
  dat$ID <- as.character(dat$ID)
  dat
}

iobr_numeric_feature_columns <- function(dat, exclude = character(0)) {
  dat <- as.data.frame(dat)
  numeric_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  setdiff(numeric_cols, exclude)
}

iobr_merge_feature_tables <- function(tables) {
  tables <- Filter(function(x) is.data.frame(x) && "ID" %in% colnames(x), tables)
  if (length(tables) == 0L) {
    return(data.frame())
  }
  Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), tables)
}

iobr_add_atf3_expression <- function(feature_table, expr, target_gene) {
  if (!target_gene %in% rownames(expr)) {
    stop("Target gene not found in expression matrix: ", target_gene)
  }
  atf3 <- data.frame(
    ID = colnames(expr),
    Target_Gene = target_gene,
    Target_Expression = as.numeric(expr[target_gene, ]),
    stringsAsFactors = FALSE
  )
  merge(atf3, feature_table, by = "ID", all.x = TRUE)
}

iobr_batch_cor <- function(dat, target, features, method = "spearman") {
  if (requireNamespace("IOBR", quietly = TRUE)) {
    res <- try(IOBR::batch_cor(data = dat, target = target, feature = features, method = method), silent = TRUE)
    if (!inherits(res, "try-error")) {
      res <- as.data.frame(res, stringsAsFactors = FALSE)
      res$Method_Source <- "IOBR::batch_cor"
      return(res)
    }
  }

  rows <- lapply(features, function(feature) {
    ok <- is.finite(dat[[target]]) & is.finite(dat[[feature]])
    if (sum(ok) < 3L) {
      return(NULL)
    }
    test <- suppressWarnings(cor.test(dat[[target]][ok], dat[[feature]][ok], method = method))
    data.frame(
      sig_names = feature,
      p.value = test$p.value,
      statistic = unname(test$estimate),
      Method_Source = "stats::cor.test",
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, rows)
  if (is.null(res)) {
    return(data.frame())
  }
  res$p.adj <- p.adjust(res$p.value, method = "BH")
  res[order(res$p.value), , drop = FALSE]
}

iobr_group_test <- function(dat, group_col, features) {
  groups <- unique(dat[[group_col]][!is.na(dat[[group_col]]) & dat[[group_col]] != ""])
  if (length(groups) < 2L) {
    return(data.frame())
  }

  if (requireNamespace("IOBR", quietly = TRUE)) {
    res <- tryCatch({
      if (length(groups) == 2L) {
        IOBR::batch_wilcoxon(data = dat, target = group_col, feature = features)
      } else {
        IOBR::batch_kruskal(data = dat, group = group_col, feature = features)
      }
    }, error = function(error) error)
    if (!inherits(res, "error")) {
      res <- as.data.frame(res, stringsAsFactors = FALSE)
      res$Method_Source <- if (length(groups) == 2L) "IOBR::batch_wilcoxon" else "IOBR::batch_kruskal"
      return(res)
    }
  }

  rows <- lapply(features, function(feature) {
    formula <- stats::as.formula(paste0("`", feature, "` ~ `", group_col, "`"))
    if (length(groups) == 2L) {
      test <- suppressWarnings(wilcox.test(formula, data = dat))
      stat <- unname(test$statistic)
      source <- "stats::wilcox.test"
    } else {
      test <- suppressWarnings(kruskal.test(formula, data = dat))
      stat <- unname(test$statistic)
      source <- "stats::kruskal.test"
    }
    data.frame(sig_names = feature, p.value = test$p.value, statistic = stat, Method_Source = source)
  })
  res <- do.call(rbind, rows)
  res$p.adj <- p.adjust(res$p.value, method = "BH")
  res[order(res$p.value), , drop = FALSE]
}

iobr_make_top_barplot <- function(dat, feature_col, value_col, title, xlab = NULL, top_n = 30) {
  dat <- dat[is.finite(dat[[value_col]]), , drop = FALSE]
  dat <- dat[order(abs(dat[[value_col]]), decreasing = TRUE), , drop = FALSE]
  dat <- head(dat, top_n)
  dat[[feature_col]] <- factor(dat[[feature_col]], levels = rev(dat[[feature_col]]))

  ggplot2::ggplot(dat, ggplot2::aes(x = .data[[feature_col]], y = .data[[value_col]])) +
    ggplot2::geom_col(fill = "#2166AC", alpha = 0.86) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = xlab, y = value_col, title = title) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = TEXT_FONT_FACE),
      axis.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE)
    )
}

iobr_make_feature_boxplot <- function(dat, feature, group_col, title) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[[group_col]], y = .data[[feature]], fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.82) +
    ggplot2::geom_jitter(width = 0.12, size = 1.8, alpha = 0.78) +
    ggplot2::labs(x = NULL, y = feature, title = title) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = TEXT_FONT_FACE),
      axis.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE)
    )
}

iobr_make_pca_plot <- function(dat, features, color_col, title) {
  features <- features[vapply(features, function(feature) {
    values <- suppressWarnings(as.numeric(dat[[feature]]))
    sum(is.finite(values)) >= 3L && stats::sd(values, na.rm = TRUE) > 0
  }, logical(1))]
  if (length(features) < 2L) {
    stop("Too few non-constant features for PCA.")
  }

  complete <- stats::complete.cases(dat[, features, drop = FALSE])
  pca_data <- dat[complete, c("ID", color_col, features), drop = FALSE]
  if (nrow(pca_data) < 3L || length(features) < 2L) {
    stop("Too few complete rows/features for PCA.")
  }
  matrix_data <- scale(as.matrix(pca_data[, features, drop = FALSE]))
  pca <- prcomp(matrix_data, center = FALSE, scale. = FALSE)
  plot_data <- data.frame(
    ID = pca_data$ID,
    Group = pca_data[[color_col]],
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

  ggplot2::ggplot(plot_data, ggplot2::aes(PC1, PC2, color = Group)) +
    ggplot2::geom_point(size = 3.2, alpha = 0.86) +
    ggplot2::labs(
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)"),
      title = title
    ) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = TEXT_FONT_FACE),
      axis.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE)
    )
}

iobr_make_heatmap_plot <- function(dat, features, sample_order, title) {
  features <- features[features %in% colnames(dat)]
  dat <- dat[match(sample_order, dat$ID), c("ID", features), drop = FALSE]
  long <- reshape(
    dat,
    varying = features,
    v.names = "Value",
    timevar = "Feature",
    times = features,
    direction = "long"
  )
  long$ID <- factor(long$ID, levels = sample_order)
  long$Feature <- factor(long$Feature, levels = rev(features))

  ggplot2::ggplot(long, ggplot2::aes(x = ID, y = Feature, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D73027") +
    ggplot2::labs(x = NULL, y = NULL, title = title, fill = "Score") +
    ggplot2::theme_bw(base_size = 9, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = TEXT_FONT_FACE),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = ggplot2::element_text(size = 7, face = TEXT_FONT_FACE),
      panel.grid = ggplot2::element_blank()
    )
}

iobr_feature_table_from_score <- function(score_table, prefix = NULL) {
  score_table <- iobr_normalize_id_column(score_table)
  numeric_cols <- iobr_numeric_feature_columns(score_table)
  if (!is.null(prefix) && nzchar(prefix)) {
    rename_cols <- setdiff(colnames(score_table), "ID")
    colnames(score_table)[match(rename_cols, colnames(score_table))] <- paste0(prefix, "__", rename_cols)
  }
  score_table
}


# 4. IOBR快速分析缓存与通用任务函数 ------------------------------------------

iobr_cache_write <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qs_save(object, file)
  } else {
    saveRDS(object, file)
  }
  invisible(file)
}

iobr_cache_read <- function(file) {
  if (requireNamespace("qs2", quietly = TRUE)) {
    return(qs2::qs_read(file))
  }
  readRDS(file)
}

iobr_make_cache_file <- function(cache_root, file_stem) {
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  extension <- if (requireNamespace("qs2", quietly = TRUE)) ".qs2" else ".rds"
  file.path(cache_root, paste0(iobr_sanitize(file_stem), extension))
}

iobr_table_dir <- function(output_root, module) {
  path <- file.path(output_root, "tables", iobr_sanitize(module))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

iobr_plot_dir <- function(output_root, module, format = c("pdf", "png")) {
  format <- match.arg(format)
  path <- file.path(output_root, "plots", iobr_sanitize(module), format)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

iobr_write_module_csv <- function(dat, output_root, module, file_stem, n_rows = 21) {
  iobr_write_csv(
    dat,
    file.path(iobr_table_dir(output_root, module), paste0(iobr_sanitize(file_stem), ".csv")),
    n_rows = n_rows,
    na = "NA"
  )
}

iobr_save_module_plot <- function(plot, output_root, module, file_stem, width = 8, height = 5.5) {
  pdf_file <- file.path(
    iobr_plot_dir(output_root, module, "pdf"),
    paste0(iobr_sanitize(file_stem), ".pdf")
  )
  save_ggplot_pdf_png(plot, pdf_file = pdf_file, width = width, height = height)
}

iobr_append_title <- function(plot, title) {
  if (inherits(plot, "ggplot")) {
    return(plot + ggplot2::ggtitle(title) + ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = TEXT_FONT_FACE)
    ))
  }
  plot
}

iobr_empty_feature_result <- function() {
  data.frame(
    sig_names = character(0),
    p.value = numeric(0),
    statistic = numeric(0),
    p.adj = numeric(0),
    Method_Source = character(0),
    stringsAsFactors = FALSE
  )
}

iobr_standardize_test_result <- function(dat, feature_col = NULL) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(dat) == 0L) {
    return(iobr_empty_feature_result())
  }

  if (is.null(feature_col)) {
    feature_col <- intersect(
      c("sig_names", "feature", "Feature", "signature", "Signature", "variables", "Variable"),
      colnames(dat)
    )[1]
  }
  if (is.na(feature_col) || is.null(feature_col)) {
    feature_col <- colnames(dat)[1]
  }
  if (!identical(feature_col, "sig_names")) {
    colnames(dat)[match(feature_col, colnames(dat))] <- "sig_names"
  }

  p_col <- intersect(c("p.value", "p.value.Wilcox", "pvalue", "p", "P.Value"), colnames(dat))[1]
  if (is.na(p_col)) {
    dat$p.value <- NA_real_
  } else if (!identical(p_col, "p.value")) {
    colnames(dat)[match(p_col, colnames(dat))] <- "p.value"
  }

  stat_col <- intersect(c("statistic", "estimate", "cor", "Correlation", "rho", "R"), colnames(dat))[1]
  if (is.na(stat_col)) {
    dat$statistic <- NA_real_
  } else if (!identical(stat_col, "statistic")) {
    colnames(dat)[match(stat_col, colnames(dat))] <- "statistic"
  }

  if (!"p.adj" %in% colnames(dat)) {
    padj_col <- intersect(c("p.adj", "p.adjust", "p.adjusted", "FDR", "adj.P.Val"), colnames(dat))[1]
    if (!is.na(padj_col)) {
      colnames(dat)[match(padj_col, colnames(dat))] <- "p.adj"
    } else {
      dat$p.adj <- p.adjust(suppressWarnings(as.numeric(dat$p.value)), method = "BH")
    }
  }

  dat$sig_names <- as.character(dat$sig_names)
  dat$p.value <- suppressWarnings(as.numeric(dat$p.value))
  dat$statistic <- suppressWarnings(as.numeric(dat$statistic))
  dat$p.adj <- suppressWarnings(as.numeric(dat$p.adj))
  dat[order(dat$p.value, na.last = TRUE), , drop = FALSE]
}

iobr_get_signature_groups <- function(selected_groups = NULL) {
  e <- new.env(parent = emptyenv())
  data(sig_group, package = "IOBR", envir = e)
  groups <- e$sig_group
  if (!is.null(selected_groups) && length(selected_groups) > 0L) {
    selected_groups <- intersect(selected_groups, names(groups))
    groups <- groups[selected_groups]
  }
  groups
}

iobr_get_signature_collection_list <- function(selected_signatures = NULL) {
  e <- new.env(parent = emptyenv())
  data(signature_collection, package = "IOBR", envir = e)
  signatures <- e$signature_collection
  if (!is.null(selected_signatures) && length(selected_signatures) > 0L) {
    selected_signatures <- intersect(selected_signatures, names(signatures))
    signatures <- signatures[selected_signatures]
  }
  signatures
}

iobr_select_signatures <- function(
    selected_groups = NULL,
    selected_signatures = NULL,
    max_signatures = Inf) {
  signatures <- iobr_get_signature_collection_list(selected_signatures)

  if (!is.null(selected_groups) && length(selected_groups) > 0L) {
    groups <- iobr_get_signature_groups(selected_groups)
    wanted <- unique(unlist(groups, use.names = FALSE))
    signatures <- signatures[intersect(wanted, names(signatures))]
  }

  if (is.finite(max_signatures) && length(signatures) > max_signatures) {
    signatures <- signatures[seq_len(max_signatures)]
  }

  signatures
}

iobr_write_signature_catalog <- function(output_root) {
  signatures <- iobr_get_signature_collection_list()
  groups <- iobr_get_signature_groups()

  signature_table <- data.frame(
    Signature = names(signatures),
    Gene_Count = vapply(signatures, length, integer(1)),
    stringsAsFactors = FALSE
  )
  group_table <- do.call(rbind, lapply(names(groups), function(group_name) {
    data.frame(
      Signature_Group = group_name,
      Signature = groups[[group_name]],
      stringsAsFactors = FALSE
    )
  }))

  iobr_write_module_csv(signature_table, output_root, "capability_catalog", "001_iobr_signature_collection")
  iobr_write_module_csv(group_table, output_root, "capability_catalog", "002_iobr_signature_groups")
}

iobr_prepare_and_cache_local_input <- function(
    project_root,
    cache_root,
    dataset_id,
    se_file,
    target_gene,
    project_id = dataset_id,
    assay_name = "tpm",
    deconvolution_assay = assay_name,
    score_assay = assay_name,
    score_log2 = TRUE,
    sample_detail_values = NULL,
    sample_detail_column = "group_detail",
    sample_barcode_pattern = NULL,
    clinical_file = NULL,
    tumor = TRUE,
    timer_indication = NULL) {
  input <- iobr_prepare_local_se_input(
    dataset_id = dataset_id,
    se_file = se_file,
    target_gene = target_gene,
    project_id = project_id,
    assay_name = assay_name,
    deconvolution_assay = deconvolution_assay,
    score_assay = score_assay,
    score_log2 = score_log2,
    sample_detail_values = sample_detail_values,
    sample_detail_column = sample_detail_column,
    sample_barcode_pattern = sample_barcode_pattern,
    clinical_file = clinical_file,
    tumor = tumor,
    timer_indication = timer_indication
  )
  input$pdata <- iobr_add_survival_columns(input$pdata)

  cache_file <- iobr_make_cache_file(
    cache_root,
    paste(dataset_id, target_gene, "iobr_input", sep = "_")
  )
  iobr_cache_write(input, cache_file)

  data.frame(
    Dataset_ID = dataset_id,
    ProjectID = project_id,
    Target_Gene = target_gene,
    Sample_Count = input$sample_count,
    Gene_Count = input$gene_count,
    SE_File = normalizePath(se_file, winslash = "/", mustWork = TRUE),
    Cache_File = cache_file,
    Tumor = tumor,
    TIMER_Indication = ifelse(is.null(timer_indication), NA_character_, timer_indication),
    stringsAsFactors = FALSE
  )
}

iobr_get_group_columns <- function(input, include_target_group = TRUE, gse_designs = TRUE) {
  pdata <- input$pdata
  group_columns <- character(0)

  if (include_target_group && "Target_Group" %in% colnames(pdata)) {
    group_columns <- c(group_columns, "Target_Group")
  }
  if (gse_designs) {
    design_cols <- grep("^analysis_", colnames(pdata), value = TRUE)
    group_columns <- c(group_columns, design_cols)
  }

  unique(group_columns[vapply(group_columns, function(column) {
    values <- pdata[[column]]
    values <- values[!is.na(values) & trimws(as.character(values)) != ""]
    length(unique(values)) >= 2L
  }, logical(1))])
}

iobr_make_group_test_table <- function(feature_table, group_columns, features) {
  if (length(group_columns) == 0L || length(features) == 0L) {
    return(data.frame())
  }

  results <- lapply(group_columns, function(group_col) {
    res <- iobr_standardize_test_result(
      iobr_group_test(feature_table, group_col = group_col, features = features)
    )
    if (nrow(res) == 0L) {
      return(NULL)
    }
    res$Group_Column <- group_col
    res
  })

  iobr_bind_rows(results)
}

iobr_make_correlation_table <- function(feature_table, target = "Target_Expression", features, method = "spearman") {
  if (!target %in% colnames(feature_table) || length(features) == 0L) {
    return(data.frame())
  }
  iobr_standardize_test_result(
    iobr_batch_cor(feature_table, target = target, features = features, method = method)
  )
}

iobr_save_feature_summary_plots <- function(
    feature_table,
    output_root,
    module,
    file_prefix,
    group_columns = character(0),
    correlation_table = data.frame(),
    group_test_table = data.frame(),
    top_n = 15L) {
  features <- iobr_feature_columns(feature_table)
  status <- list()

  try({
  if (nrow(correlation_table) > 0L && "statistic" %in% colnames(correlation_table)) {
    plot_data <- head(
      correlation_table[order(abs(correlation_table$statistic), decreasing = TRUE), , drop = FALSE],
      top_n
    )
    if (nrow(plot_data) > 0L) {
      p <- iobr_make_top_barplot(
        plot_data,
        feature_col = "sig_names",
        value_col = "statistic",
        title = paste0(file_prefix, " correlation with target gene"),
        xlab = NULL,
        top_n = top_n
      )
      files <- iobr_save_module_plot(
        p,
        output_root = output_root,
        module = module,
        file_stem = paste0(file_prefix, "_target_correlation_top", top_n),
        width = 8.5,
        height = max(5.5, 0.24 * nrow(plot_data) + 2.4)
      )
      status[[length(status) + 1L]] <- files
    }
  }
  }, silent = TRUE)

  try({
  if (nrow(group_test_table) > 0L) {
    for (group_col in unique(group_test_table$Group_Column)) {
      current <- group_test_table[group_test_table$Group_Column == group_col, , drop = FALSE]
      current$minus_log10_p <- -log10(pmax(current$p.value, .Machine$double.xmin))
      current <- head(current[order(current$p.value, na.last = TRUE), , drop = FALSE], top_n)
      if (nrow(current) == 0L) {
        next
      }

      p <- iobr_make_top_barplot(
        current,
        feature_col = "sig_names",
        value_col = "minus_log10_p",
        title = paste0(file_prefix, " group differences: ", group_col),
        xlab = NULL,
        top_n = top_n
      )
      iobr_save_module_plot(
        p,
        output_root = output_root,
        module = module,
        file_stem = paste0(file_prefix, "_", group_col, "_group_difference_top", top_n),
        width = 8.5,
        height = max(5.5, 0.24 * nrow(current) + 2.4)
      )

      box_features <- intersect(current$sig_names, features)
      if (length(box_features) > 0L) {
        for (feature in head(box_features, min(6L, length(box_features)))) {
          p_box <- try(
            IOBR::sig_box(
              data = feature_table,
              signature = feature,
              variable = group_col,
              jitter = TRUE,
              point_size = 1.8,
              size_of_font = 10,
              size_of_pvalue = 4.5,
              show_pairwise_p = TRUE,
              show_overall_p = FALSE
            ),
            silent = TRUE
          )
          if (inherits(p_box, "try-error") || !inherits(p_box, "ggplot")) {
            p_box <- iobr_make_feature_boxplot(
              feature_table,
              feature = feature,
              group_col = group_col,
              title = paste(file_prefix, feature, "by", group_col)
            )
          } else {
            p_box <- iobr_append_title(p_box, paste(file_prefix, feature, "by", group_col))
          }
          iobr_save_module_plot(
            p_box,
            output_root = output_root,
            module = module,
            file_stem = paste0(file_prefix, "_", group_col, "_", feature, "_boxplot"),
            width = 5.8,
            height = 5.2
          )
        }
      }
    }
  }
  }, silent = TRUE)

  pca_features <- character(0)
  if (nrow(correlation_table) > 0L) {
    pca_features <- correlation_table$sig_names[order(abs(correlation_table$statistic), decreasing = TRUE)]
  }
  if (length(pca_features) < 2L && nrow(group_test_table) > 0L) {
    pca_features <- c(pca_features, group_test_table$sig_names[order(group_test_table$p.value, na.last = TRUE)])
  }
  pca_features <- unique(intersect(pca_features, features))
  pca_features <- head(pca_features, min(top_n, length(pca_features)))
  pca_features <- pca_features[vapply(pca_features, function(feature) {
    values <- suppressWarnings(as.numeric(feature_table[[feature]]))
    sum(is.finite(values)) >= 3L && stats::sd(values, na.rm = TRUE) > 0
  }, logical(1))]

  try({
  if (length(pca_features) >= 2L && length(group_columns) > 0L) {
    group_col <- group_columns[1]
    p_pca <- try(
      IOBR::iobr_pca(
        data = t(as.matrix(feature_table[, pca_features, drop = FALSE])),
        is.matrix = TRUE,
        scale = TRUE,
        is.log = FALSE,
        pdata = feature_table[, c("ID", group_col), drop = FALSE],
        id_pdata = "ID",
        group = group_col,
        repel = FALSE,
        addEllipses = FALSE
      ),
      silent = TRUE
    )
    if (inherits(p_pca, "try-error") || !inherits(p_pca, "ggplot")) {
      p_pca <- iobr_make_pca_plot(
        feature_table,
        features = pca_features,
        color_col = group_col,
        title = paste0(file_prefix, " PCA")
      )
    } else {
      p_pca <- iobr_append_title(p_pca, paste0(file_prefix, " PCA"))
    }
    iobr_save_module_plot(
      p_pca,
      output_root = output_root,
      module = module,
      file_stem = paste0(file_prefix, "_", group_col, "_pca"),
      width = 6.5,
      height = 5.8
    )

    sample_order <- feature_table$ID[order(feature_table[[group_col]], feature_table$Target_Expression)]
    p_heat <- iobr_make_heatmap_plot(
      feature_table,
      features = pca_features,
      sample_order = sample_order,
      title = paste0(file_prefix, " top feature heatmap")
    )
    iobr_save_module_plot(
      p_heat,
      output_root = output_root,
      module = module,
      file_stem = paste0(file_prefix, "_", group_col, "_heatmap"),
      width = max(7, min(18, 0.065 * length(sample_order) + 4.5)),
      height = max(5.5, min(12, 0.26 * length(pca_features) + 2.8))
    )
  }
  }, silent = TRUE)

  invisible(status)
}

iobr_run_deconvolution_task <- function(task, output_root, correlation_method = "spearman", top_n = 15L) {
  input <- iobr_cache_read(task$Input_File)
  method <- task$Method
  module <- "deconvolution"
  file_prefix <- iobr_sanitize(paste(task$Dataset_ID, task$Target_Gene, method, sep = "_"))

  result <- iobr_safe_task({
    reference <- NULL
    if (method %in% c("svr", "lsei")) {
      reference <- iobr_load_lm22_reference()
    }
    group_list <- NULL
    if (identical(method, "timer") && !is.null(input$timer_indication)) {
      group_list <- rep(input$timer_indication, ncol(input$deconvolution_expr))
    }

    deconv <- suppressWarnings(suppressMessages(IOBR::deconvo_tme(
      eset = input$deconvolution_expr,
      project = input$project_id,
      method = method,
      tumor = input$tumor,
      perm = task$Perm,
      reference = reference,
      group_list = group_list,
      arrays = FALSE,
      plot = FALSE
    )))
    deconv <- iobr_normalize_id_column(deconv)
    feature_table <- merge(input$pdata, deconv, by = "ID", all.x = TRUE)
    features <- iobr_feature_columns(feature_table, exclude = setdiff(colnames(input$pdata), "ID"))
    group_columns <- iobr_get_group_columns(input)
    cor_table <- iobr_make_correlation_table(
      feature_table,
      target = "Target_Expression",
      features = features,
      method = correlation_method
    )
    group_table <- iobr_make_group_test_table(
      feature_table,
      group_columns = group_columns,
      features = features
    )

    feature_file <- iobr_write_module_csv(feature_table, output_root, module, paste0(file_prefix, "_scores"))
    cor_file <- iobr_write_module_csv(cor_table, output_root, module, paste0(file_prefix, "_target_correlation"))
    group_file <- iobr_write_module_csv(group_table, output_root, module, paste0(file_prefix, "_group_tests"))

    iobr_save_feature_summary_plots(
      feature_table = feature_table,
      output_root = output_root,
      module = module,
      file_prefix = file_prefix,
      group_columns = group_columns,
      correlation_table = cor_table,
      group_test_table = group_table,
      top_n = top_n
    )

    if (length(features) > 1L && nrow(feature_table) > 1L) {
      cell_plot <- try(
        IOBR::cell_bar_plot(
          input = feature_table[, c("ID", head(features, min(length(features), 20L))), drop = FALSE],
          id = "ID",
          title = paste0(file_prefix, " cell fractions"),
          features = head(features, min(length(features), 20L)),
          coord_flip = TRUE
        ),
        silent = TRUE
      )
      if (!inherits(cell_plot, "try-error") && inherits(cell_plot, "ggplot")) {
        iobr_save_module_plot(
          iobr_append_title(cell_plot, paste0(file_prefix, " cell fractions")),
          output_root = output_root,
          module = module,
          file_stem = paste0(file_prefix, "_cell_barplot"),
          width = 9,
          height = 6
        )
      }
    }

    iobr_status_row(
      dataset_id = task$Dataset_ID,
      module = module,
      task = paste("deconvo_tme", method, sep = "::"),
      status = "success",
      output_file = feature_file,
      n_rows = nrow(feature_table),
      n_cols = ncol(feature_table)
    )
  }, module = module, task = paste(task$Dataset_ID, method, sep = "_"))

  result
}

iobr_run_signature_score_task <- function(
    task,
    output_root,
    signatures,
    correlation_method = "spearman",
    top_n = 15L) {
  input <- iobr_cache_read(task$Input_File)
  method <- task$Method
  module <- "signature_score"
  file_prefix <- iobr_sanitize(paste(task$Dataset_ID, task$Target_Gene, method, sep = "_"))

  result <- iobr_safe_task({
    score <- suppressWarnings(suppressMessages(IOBR::calculate_sig_score(
      pdata = input$pdata[, c("ID"), drop = FALSE],
      eset = input$score_expr,
      signature = signatures,
      method = method,
      mini_gene_count = task$Mini_Gene_Count,
      column_of_sample = "ID",
      print_gene_proportion = FALSE,
      print_filtered_signatures = FALSE,
      adjust_eset = FALSE,
      parallel.size = task$Inner_Workers
    )))
    score <- iobr_normalize_id_column(score)
    feature_table <- merge(input$pdata, score, by = "ID", all.x = TRUE)
    features <- iobr_feature_columns(feature_table, exclude = setdiff(colnames(input$pdata), "ID"))
    group_columns <- iobr_get_group_columns(input)
    cor_table <- iobr_make_correlation_table(
      feature_table,
      target = "Target_Expression",
      features = features,
      method = correlation_method
    )
    group_table <- iobr_make_group_test_table(
      feature_table,
      group_columns = group_columns,
      features = features
    )

    feature_file <- iobr_write_module_csv(feature_table, output_root, module, paste0(file_prefix, "_scores"))
    cor_file <- iobr_write_module_csv(cor_table, output_root, module, paste0(file_prefix, "_target_correlation"))
    group_file <- iobr_write_module_csv(group_table, output_root, module, paste0(file_prefix, "_group_tests"))

    iobr_save_feature_summary_plots(
      feature_table = feature_table,
      output_root = output_root,
      module = module,
      file_prefix = file_prefix,
      group_columns = group_columns,
      correlation_table = cor_table,
      group_test_table = group_table,
      top_n = top_n
    )

    iobr_status_row(
      dataset_id = task$Dataset_ID,
      module = module,
      task = paste("calculate_sig_score", method, sep = "::"),
      status = "success",
      output_file = feature_file,
      n_rows = nrow(feature_table),
      n_cols = ncol(feature_table)
    )
  }, module = module, task = paste(task$Dataset_ID, method, sep = "_"))

  result
}
