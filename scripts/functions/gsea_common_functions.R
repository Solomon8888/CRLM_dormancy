# GSEA公共配置声明和函数
#
# 本文件由06号GSEA脚本中的内部声明与常用函数整理而来。
# 06号脚本负责GSEA运算，07号脚本负责GSEA绘图；两者共同调用这里的函数。

GENE_ID_COLUMN_BY_TYPE <- c(
  ENTREZ = "Entrez",
  SYMBOL = "Symbol",
  ENSEMBL = "Ensembl"
)

READABLE_KEY_TYPE_BY_GENE_ID_TYPE <- c(
  ENTREZ = "ENTREZID",
  SYMBOL = "SYMBOL",
  ENSEMBL = "ENSEMBL"
)

ORGDB_PACKAGE_BY_SPECIES <- c(
  human = "org.Hs.eg.db",
  `homo sapiens` = "org.Hs.eg.db",
  mouse = "org.Mm.eg.db",
  `mus musculus` = "org.Mm.eg.db",
  rat = "org.Rn.eg.db",
  `rattus norvegicus` = "org.Rn.eg.db"
)

MSIGDB_GENE_COLUMN_BY_TYPE <- c(
  ENTREZ = "ncbi_gene",
  SYMBOL = "gene_symbol",
  ENSEMBL = "ensembl_gene"
)

BREAK_RANK_TIES <- TRUE
RANK_TIE_EPSILON <- 1e-12
GSEA_RANDOM_SEED <- 20260604

GSEA_WARNING_PATTERNS_TO_HIDE <- c(
  "P-values are less than",
  "P-values were not calculated properly",
  "Invalid p-values detected",
  "NA values detected in gene set IDs",
  "Duplicate gene set IDs detected",
  "aes_string() was deprecated",
  "aes_() was deprecated",
  "Using `size` aesthetic for lines was deprecated"
)

GSEA_RESULT_COLUMNS <- c(
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)


# 2. 常用函数 -----------------------------------------------------------------

get_msigdb_db_species <- function(species) {
  # msigdbr区分MSigDB数据库物种(db_species)和目标物种(species)。
  # 为减少头部配置，人类/大鼠默认使用人类MSigDB库，小鼠使用小鼠MSigDB库。
  species_lower <- tolower(trimws(species))

  if (species_lower %in% c("mouse", "mus musculus", "mm", "mmusculus")) {
    return("MM")
  }

  "HS"
}

get_orgdb_package <- function(species) {
  species_key <- tolower(trimws(species))
  orgdb_package <- ORGDB_PACKAGE_BY_SPECIES[[species_key]]

  if (is.null(orgdb_package)) {
    stop("No OrgDb package is configured for SPECIES = ", species)
  }

  orgdb_package
}

get_orgdb_object <- function(species) {
  orgdb_package <- get_orgdb_package(species)
  if (!requireNamespace(orgdb_package, quietly = TRUE)) {
    stop(
      "READABLE_GENE_SYMBOLS requires package ",
      orgdb_package,
      ". Please install it before running GSEA."
    )
  }

  get(orgdb_package, envir = asNamespace(orgdb_package))
}

get_readable_key_type <- function(gene_id_type) {
  key_type <- READABLE_KEY_TYPE_BY_GENE_ID_TYPE[[gene_id_type]]
  if (is.null(key_type)) {
    stop("No readable keyType is configured for GENE_ID_TYPE = ", gene_id_type)
  }

  key_type
}

make_msigdb_key <- function(collection, subcollection) {
  if (is.na(subcollection) || subcollection == "") {
    return(collection)
  }

  paste(collection, subcollection, sep = ":")
}

make_msigdb_output_name <- function(collection, subcollection) {
  if (is.na(subcollection) || subcollection == "") {
    if (collection == "H") {
      return("hallmark")
    }
    return(collection)
  }

  # 有子类别时优先用子类别命名输出目录，例如GO:BP保存为GO_BP。
  # 目录结构已经位于GSEA下，不再重复写入C5等大类前缀。
  gsub("[:/ -]+", "_", subcollection)
}

build_msigdb_geneset_catalog <- function() {
  catalog_table <- load_msigdb_catalog_table()

  catalog_list <- lapply(seq_len(nrow(catalog_table)), function(i) {
    collection <- as.character(catalog_table$gs_collection[i])
    subcollection <- as.character(catalog_table$gs_subcollection[i])
    if (is.na(subcollection)) {
      subcollection <- ""
    }

    list(
      key = make_msigdb_key(collection, subcollection),
      collection = collection,
      subcollection = if (subcollection == "") NULL else subcollection,
      output_name = make_msigdb_output_name(collection, subcollection),
      description = as.character(catalog_table$gs_collection_name[i])
    )
  })

  names(catalog_list) <- vapply(catalog_list, `[[`, character(1), "key")
  catalog_list
}

select_msigdb_genesets <- function(catalog, genesets_to_run) {
  genesets_to_run <- unique(trimws(as.character(genesets_to_run)))

  if (length(genesets_to_run) == 1 && tolower(genesets_to_run) == "all") {
    selected_catalog <- catalog
  } else {
    missing_genesets <- setdiff(genesets_to_run, names(catalog))
    if (length(missing_genesets) > 0) {
      stop(
        "Undefined MSigDB gene set keys in GSEA_GENESETS_TO_RUN: ",
        paste(missing_genesets, collapse = ", ")
      )
    }

    selected_catalog <- catalog[genesets_to_run]
  }

  names(selected_catalog) <- vapply(
    selected_catalog,
    `[[`,
    character(1),
    "output_name"
  )

  if (any(duplicated(names(selected_catalog)))) {
    names(selected_catalog) <- make.unique(names(selected_catalog), sep = "_")
  }

  selected_catalog
}

get_runtime_genesets_to_run <- function() {
  # 日常运行使用脚本头部GSEA_GENESETS_TO_RUN。
  # 仅测试脚本时，可在命令行临时设置环境变量：
  # GSEA_TEST_GENESETS="H,C6" Rscript scripts/GSE114012/06_gsea_analysis.R
  # 这样无需修改脚本即可只跑少数基因集。
  test_genesets <- Sys.getenv("GSEA_TEST_GENESETS", unset = "")
  test_genesets <- trimws(test_genesets)

  if (test_genesets == "") {
    return(GSEA_GENESETS_TO_RUN)
  }

  unique(trimws(strsplit(test_genesets, ",", fixed = TRUE)[[1]]))
}

clean_previous_gsea_outputs <- function(selected_analyses) {
  # 每次重跑都清空所选分析设计的GSEA表格和图片结果，避免旧结果残留。
  # ANALYSES_TO_RUN为"all"时，selected_analyses即全部差异分析设计。
  if (!CLEAN_GSEA_OUTPUT_DIR) {
    return(invisible(FALSE))
  }

  table_dirs <- file.path(TABLE_OUTPUT_ROOT, selected_analyses, "GSEA")
  plot_dirs <- file.path(PLOT_ROOT, sanitize_file_name(selected_analyses))

  unlink(table_dirs, recursive = TRUE, force = TRUE)
  unlink(plot_dirs, recursive = TRUE, force = TRUE)

  # 若当前运行覆盖全部分析，则直接清理GSEA图片根目录下的陈旧散落文件/目录。
  if (identical(ANALYSES_TO_RUN, "all") && dir.exists(PLOT_ROOT)) {
    unlink(PLOT_ROOT, recursive = TRUE, force = TRUE)
  }

  invisible(TRUE)
}

with_gsea_warnings_suppressed <- function(expr) {
  # 定向静默批量运行中反复出现、但不影响结果的第三方包提示。
  withCallingHandlers(
    expr,
    warning = function(w) {
      warning_text <- conditionMessage(w)
      if (any(vapply(
        GSEA_WARNING_PATTERNS_TO_HIDE,
        grepl,
        logical(1),
        x = warning_text,
        fixed = TRUE
      ))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

object_signature <- function(x) {
  paste(capture.output(dput(x)), collapse = "")
}

cache_hash <- function(...) {
  # 使用md5短哈希避免缓存文件名过长，同时让配置变化能自动生成新缓存。
  text <- paste(..., collapse = "|")
  tmp_file <- tempfile("gsea_cache_hash_")
  writeLines(text, tmp_file, useBytes = TRUE)
  hash_value <- unname(tools::md5sum(tmp_file))
  unlink(tmp_file)
  hash_value
}

get_file_stamp <- function(file) {
  file_info <- file.info(file)
  paste(file_info$size, as.integer(file_info$mtime), sep = "_")
}

read_qs2_cache <- function(cache_file) {
  if (!USE_QS2_CACHE || REFRESH_QS2_CACHE || !file.exists(cache_file)) {
    return(NULL)
  }

  qs2::qs_read(cache_file)
}

write_qs2_cache <- function(object, cache_file) {
  if (!USE_QS2_CACHE) {
    return(invisible(FALSE))
  }

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(object, cache_file)
  invisible(TRUE)
}

is_reference_cache_fresh <- function(cache_file) {
  if (!file.exists(cache_file)) {
    return(FALSE)
  }

  cache_age_days <- as.numeric(
    difftime(Sys.time(), file.info(cache_file)$mtime, units = "days")
  )
  cache_age_days <= MSIGDB_REFERENCE_MAX_AGE_DAYS
}

read_reference_qs2_cache <- function(cache_file) {
  if (!is_reference_cache_fresh(cache_file)) {
    return(NULL)
  }

  qs2::qs_read(cache_file)
}

write_reference_qs2_cache <- function(object, cache_file) {
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(object, cache_file)
  invisible(TRUE)
}

get_msigdb_catalog_cache_file <- function() {
  file.path(
    MSIGDB_REFERENCE_DIR,
    paste0("catalog_", get_msigdb_db_species(SPECIES), ".qs2")
  )
}

load_msigdb_catalog_table <- function() {
  cache_file <- get_msigdb_catalog_cache_file()
  cached_catalog <- read_reference_qs2_cache(cache_file)
  if (!is.null(cached_catalog)) {
    return(cached_catalog)
  }

  catalog_table <- msigdbr::msigdbr_collections(
    db_species = get_msigdb_db_species(SPECIES)
  )
  write_reference_qs2_cache(catalog_table, cache_file)
  catalog_table
}

get_msigdb_data_cache_file <- function(config) {
  subcollection <- ifelse(
    is.null(config$subcollection),
    "all",
    sanitize_file_name(config$subcollection)
  )

  file.path(
    MSIGDB_REFERENCE_DIR,
    get_msigdb_db_species(SPECIES),
    sanitize_file_name(SPECIES),
    paste0(
      sanitize_file_name(config$collection),
      "_",
      subcollection,
      ".qs2"
    )
  )
}

get_gsea_cache_file <- function(analysis_name, deg_file, geneset_name, config) {
  hash_value <- cache_hash(
    "gsea",
    SPECIES,
    GENE_ID_TYPE,
    RANK_METRIC_COLUMN,
    READABLE_GENE_SYMBOLS,
    if (READABLE_GENE_SYMBOLS) get_orgdb_package(SPECIES) else "not_readable",
    get_file_stamp(deg_file),
    geneset_name,
    object_signature(config),
    object_signature(GSEA_PARAMS)
  )

  file.path(
    QS2_CACHE_DIR,
    "gsea_result",
    sanitize_file_name(analysis_name),
    paste0(sanitize_file_name(geneset_name), "_", hash_value, ".qs2")
  )
}

get_msigdb_collection <- function(config) {
  # 根据脚本头部配置获取一套msigdbr基因集。
  cache_file <- get_msigdb_data_cache_file(config)
  cached_msigdb_data <- read_reference_qs2_cache(cache_file)
  if (!is.null(cached_msigdb_data)) {
    return(list(
      data = cached_msigdb_data,
      cache_source = "reference_cache",
      cache_file = cache_file
    ))
  }

  args <- list(
    db_species = get_msigdb_db_species(SPECIES),
    species = SPECIES,
    collection = config$collection
  )

  if (!is.null(config$subcollection)) {
    args$subcollection <- config$subcollection
  }

  msigdb_data <- do.call(msigdbr::msigdbr, args)
  write_reference_qs2_cache(msigdb_data, cache_file)

  list(
    data = msigdb_data,
    cache_source = "computed",
    cache_file = cache_file
  )
}

prepare_msigdb_terms <- function(msigdb_data, gene_id_type) {
  # 将msigdbr结果整理为clusterProfiler::GSEA需要的TERM2GENE/TERM2NAME格式。
  gene_column <- MSIGDB_GENE_COLUMN_BY_TYPE[[gene_id_type]]
  stopifnot(!is.null(gene_column))
  stopifnot(gene_column %in% colnames(msigdb_data))
  stopifnot("gs_name" %in% colnames(msigdb_data))

  term2gene <- msigdb_data[, c("gs_name", gene_column), drop = FALSE]
  colnames(term2gene) <- c("term", "gene")

  term2gene$term <- as.character(term2gene$term)
  term2gene$gene <- as.character(term2gene$gene)
  term2gene <- term2gene[
    !is.na(term2gene$term) & term2gene$term != "" &
      !is.na(term2gene$gene) & term2gene$gene != "",
    ,
    drop = FALSE
  ]
  term2gene <- unique(term2gene)

  term2name <- unique(data.frame(
    term = as.character(msigdb_data$gs_name),
    name = as.character(msigdb_data$gs_name),
    stringsAsFactors = FALSE
  ))

  list(
    term2gene = term2gene,
    term2name = term2name
  )
}

load_msigdb_terms <- function(geneset_name, config) {
  msigdb_cache <- get_msigdb_collection(config)
  terms <- prepare_msigdb_terms(msigdb_cache$data, GENE_ID_TYPE)
  terms$config <- config
  terms$Cache_Source <- msigdb_cache$cache_source
  terms$Reference_Cache_File <- msigdb_cache$cache_file
  terms
}

count_nes_direction <- function(result_table, direction = c("positive", "negative")) {
  # 终端summary用；不参与官方GSEA结果表保存。
  direction <- match.arg(direction)
  if (!"NES" %in% colnames(result_table)) {
    return(0L)
  }

  nes <- as.numeric(result_table$NES)
  if (direction == "positive") {
    return(sum(nes > 0, na.rm = TRUE))
  }

  sum(nes < 0, na.rm = TRUE)
}

prepare_gsea_table_for_text_output <- function(dat) {
  # md/tex属于同一官方结果表的不同文本格式；NA按R输出习惯保留为"NA"。
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  dat <- data.frame(
    lapply(dat, function(column) {
      column <- as.character(column)
      column[is.na(column)] <- "NA"
      column
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  dat
}

write_gsea_markdown_table <- function(dat, md_file) {
  dat <- prepare_gsea_table_for_text_output(dat)

  header <- paste(escape_markdown_vector(colnames(dat)), collapse = " | ")
  separator <- paste(rep("---", ncol(dat)), collapse = " | ")
  rows <- if (nrow(dat) > 0) {
    apply(dat, 1, function(row) {
      paste(escape_markdown_vector(row), collapse = " | ")
    })
  } else {
    character(0)
  }

  writeLines(
    c(
      paste0("| ", header, " |"),
      paste0("| ", separator, " |"),
      paste0("| ", rows, " |")
    ),
    md_file,
    useBytes = TRUE
  )
}

write_gsea_latex_table <- function(dat, tex_file) {
  dat <- prepare_gsea_table_for_text_output(dat)

  col_spec <- paste(rep("l", ncol(dat)), collapse = "")
  header <- paste(escape_latex_vector(colnames(dat)), collapse = " & ")
  rows <- if (nrow(dat) > 0) {
    apply(dat, 1, function(row) {
      paste(escape_latex_vector(row), collapse = " & ")
    })
  } else {
    character(0)
  }

  writeLines(
    c(
      "\\begingroup",
      "\\tiny",
      "\\setlength{\\tabcolsep}{2pt}",
      "\\renewcommand{\\arraystretch}{1.18}",
      sprintf("\\begin{tabular}{@{}%s@{}}", col_spec),
      "\\toprule",
      paste0(header, " \\\\"),
      "\\midrule",
      if (length(rows) > 0) paste0(rows, " \\\\"),
      "\\bottomrule",
      "\\end{tabular}",
      "\\endgroup"
    ),
    tex_file,
    useBytes = TRUE
  )
}

prepare_gene_list <- function(deg_table, gene_id_type, rank_metric_column) {
  # 从all_genes.csv构建GSEA排序向量：names为基因ID，values为排序统计量。
  gene_column <- GENE_ID_COLUMN_BY_TYPE[[gene_id_type]]
  stopifnot(!is.null(gene_column))
  stopifnot(gene_column %in% colnames(deg_table))
  stopifnot(rank_metric_column %in% colnames(deg_table))

  dat <- data.frame(
    Gene_ID = trimws(as.character(deg_table[[gene_column]])),
    Rank_Metric = as.numeric(deg_table[[rank_metric_column]]),
    stringsAsFactors = FALSE
  )

  dat <- dat[
    !is.na(dat$Gene_ID) & dat$Gene_ID != "" &
      is.finite(dat$Rank_Metric),
    ,
    drop = FALSE
  ]

  stopifnot(nrow(dat) > 0)

  # 同一基因ID出现多次时，保留绝对统计量最大的记录。
  dat <- dat[order(-abs(dat$Rank_Metric), dat$Gene_ID), , drop = FALSE]
  dat <- dat[!duplicated(dat$Gene_ID), , drop = FALSE]

  if (BREAK_RANK_TIES && any(duplicated(dat$Rank_Metric))) {
    tie_breaker <- rank(dat$Gene_ID, ties.method = "first")
    dat$Rank_Metric <- dat$Rank_Metric + tie_breaker * RANK_TIE_EPSILON
  }

  gene_list <- dat$Rank_Metric
  names(gene_list) <- dat$Gene_ID
  gene_list <- sort(gene_list, decreasing = TRUE)

  stopifnot(length(gene_list) > 1)
  stopifnot(!any(duplicated(names(gene_list))))
  gene_list
}

run_clusterprofiler_gsea <- function(gene_list, term2gene, term2name, params) {
  # 只把官方GSEA支持的参数传入clusterProfiler::GSEA，保持结果字段为官方格式。
  # nproc经由...传递给fgsea/multilevel底层函数。
  # 当前脚本优先做任务级并行；若任务数少，再自动把更多核心分配给单个GSEA。
  set.seed(GSEA_RANDOM_SEED)
  params_for_run <- params
  if (!"nproc" %in% names(params_for_run)) {
    params_for_run$nproc <- GSEA_INNER_NPROC
  }

  gsea_call <- function() {
    do.call(
      clusterProfiler::GSEA,
      c(
        list(
          geneList = gene_list,
          TERM2GENE = term2gene,
          TERM2NAME = term2name
        ),
        params_for_run
      )
    )
  }

  with_gsea_warnings_suppressed(gsea_call())
}

make_gsea_result_readable <- function(gsea_result) {
  # 将core_enrichment里的基因ID转换为Symbol，方便直接阅读富集到的基因。
  if (!READABLE_GENE_SYMBOLS || GENE_ID_TYPE == "SYMBOL") {
    return(gsea_result)
  }
  if (!inherits(gsea_result, c("gseaResult", "enrichResult", "compareClusterResult"))) {
    return(gsea_result)
  }

  clusterProfiler::setReadable(
    x = gsea_result,
    OrgDb = get_orgdb_object(SPECIES),
    keyType = get_readable_key_type(GENE_ID_TYPE),
    toType = "SYMBOL"
  )
}

load_or_run_gsea <- function(
    analysis_name,
    deg_file,
    geneset_name,
    config,
    gene_list,
    term2gene,
    term2name) {
  cache_file <- get_gsea_cache_file(
    analysis_name = analysis_name,
    deg_file = deg_file,
    geneset_name = geneset_name,
    config = config
  )

  cached_result <- read_qs2_cache(cache_file)
  if (!is.null(cached_result)) {
    return(list(result = cached_result, source = "cache"))
  }

  gsea_result <- run_clusterprofiler_gsea(
    gene_list = gene_list,
    term2gene = term2gene,
    term2name = term2name,
    params = GSEA_PARAMS
  )
  gsea_result <- make_gsea_result_readable(gsea_result)
  write_qs2_cache(gsea_result, cache_file)

  list(result = gsea_result, source = "computed")
}

is_gsea_result_object <- function(gsea_result) {
  inherits(gsea_result, c("gseaResult", "enrichResult", "compareClusterResult"))
}

make_empty_gsea_result_table <- function() {
  empty_table <- as.data.frame(
    matrix(nrow = 0, ncol = length(GSEA_RESULT_COLUMNS)),
    stringsAsFactors = FALSE
  )
  colnames(empty_table) <- GSEA_RESULT_COLUMNS
  empty_table
}

get_gsea_result_table <- function(gsea_result) {
  # as.data.frame(gseaResult)保留clusterProfiler官方结果字段。
  if (is_gsea_result_object(gsea_result)) {
    return(as.data.frame(gsea_result))
  }

  if (is.data.frame(gsea_result)) {
    return(gsea_result)
  }

  make_empty_gsea_result_table()
}

write_gsea_result_tables <- function(gsea_result, csv_file) {
  result_table <- get_gsea_result_table(gsea_result)
  write.csv(result_table, csv_file, row.names = FALSE, na = "NA")
  write_gsea_markdown_table(result_table, sub("[.]csv$", ".md", csv_file))
  write_gsea_latex_table(result_table, sub("[.]csv$", ".tex", csv_file))
  result_table
}

format_gsea_description_for_plot <- function(description) {
  # 只优化绘图标签，不改变官方GSEA结果表。
  description <- as.character(description)

  if (SIMPLIFY_PATHWAY_PREFIX_IN_PLOT) {
    description <- gsub("^HALLMARKS?[_ -]+", "", description)
    description <- gsub("^GO[_ -]*BP[_ -]+", "", description)
    description <- gsub("^GOBP[_ -]+", "", description)
    description <- gsub("^GO[_ -]*CC[_ -]+", "", description)
    description <- gsub("^GOCC[_ -]+", "", description)
    description <- gsub("^GO[_ -]*MF[_ -]+", "", description)
    description <- gsub("^GOMF[_ -]+", "", description)
    description <- gsub("^KEGG[_ -]+", "", description)
    description <- gsub("^REACTOME[_ -]+", "", description)
    description <- gsub("^BIOCARTA[_ -]+", "", description)
    description <- gsub("^WP[_ -]+", "", description)
    description <- gsub("^C[0-9]+[_ -]+", "", description)
  }

  if (REPLACE_UNDERSCORE_WITH_SPACE_IN_PLOT) {
    description <- gsub("_", " ", description)
  }

  description <- gsub("\\s+", " ", description)
  trimws(description)
}

prepare_gsea_result_for_plot <- function(gsea_result) {
  # dotplot无法稳定处理ID/NES/P值为NA的通路；这里只过滤绘图对象，不改变CSV结果。
  plot_result <- gsea_result
  result_table <- as.data.frame(plot_result)

  required_columns <- c("ID", "Description", "NES", "pvalue", "p.adjust")
  missing_columns <- setdiff(required_columns, colnames(result_table))
  if (length(missing_columns) > 0) {
    return(plot_result)
  }

  valid_index <- !is.na(result_table$ID) &
    result_table$ID != "" &
    !is.na(result_table$Description) &
    result_table$Description != "" &
    is.finite(as.numeric(result_table$NES)) &
    is.finite(as.numeric(result_table$pvalue)) &
    is.finite(as.numeric(result_table$p.adjust))

  plot_result@result <- result_table[valid_index, , drop = FALSE]
  plot_result@result$Description <- format_gsea_description_for_plot(
    plot_result@result$Description
  )
  plot_result
}

get_longest_plot_label_line <- function(plot_labels) {
  if (length(plot_labels) == 0) {
    return(20L)
  }

  label_lines <- unlist(strsplit(as.character(plot_labels), "\n", fixed = TRUE))
  max(nchar(label_lines), na.rm = TRUE)
}

get_plot_label_line_counts <- function(plot_labels, shown_terms) {
  if (length(plot_labels) == 0) {
    return(rep(1L, max(as.integer(shown_terms), 1L)))
  }

  vapply(strsplit(as.character(plot_labels), "\n", fixed = TRUE), length, integer(1))
}

get_gsea_dotplot_size <- function(shown_terms, plot_labels) {
  shown_terms <- max(as.integer(shown_terms), 1L)
  label_line_counts <- get_plot_label_line_counts(plot_labels, shown_terms)
  total_label_lines <- sum(label_line_counts, na.rm = TRUE)
  label_height <- total_label_lines * DOTPLOT_LABEL_LINE_HEIGHT +
    shown_terms * DOTPLOT_TERM_GAP_HEIGHT

  body_size <- max(DOTPLOT_BODY_BASE_SIZE, label_height)
  body_size <- min(max(body_size, DOTPLOT_BODY_MIN_SIZE), DOTPLOT_BODY_MAX_SIZE)

  longest_label <- get_longest_plot_label_line(plot_labels)
  label_width <- DOTPLOT_LABEL_BASE_WIDTH +
    longest_label * DOTPLOT_LABEL_WIDTH_PER_CHARACTER
  label_width <- min(
    max(label_width, DOTPLOT_LABEL_MIN_WIDTH),
    DOTPLOT_LABEL_MAX_WIDTH
  )

  width <- label_width + body_size + DOTPLOT_LEGEND_WIDTH
  height <- body_size + DOTPLOT_VERTICAL_PADDING

  list(
    width = width,
    height = height
  )
}

make_empty_gsea_plot <- function() {
  # 没有可展示通路时仍输出空白占位图，保证批量结果目录完整。
  ggplot2::ggplot() +
    ggplot2::geom_text(
      ggplot2::aes(x = 0, y = 0, label = "No GSEA terms to display"),
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      size = 4.2
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void(base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

apply_common_gsea_theme <- function(plot) {
  suppressMessages(plot + ggplot2::scale_size(range = GSEAVIS_POINT_SIZE_RANGE)) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = TEXT_FONT_FAMILY,
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      plot.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      axis.text.y = ggplot2::element_text(
        family = TEXT_FONT_FAMILY,
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      axis.title = ggplot2::element_text(
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      legend.title = ggplot2::element_text(
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      legend.text = ggplot2::element_text(
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      panel.border = ggplot2::element_rect(
        color = TEXT_COLOR,
        linewidth = AXIS_LINE_WIDTH,
        fill = NA
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey90",
        color = TEXT_COLOR,
        linewidth = AXIS_LINE_WIDTH
      ),
      strip.text = ggplot2::element_text(
        family = TEXT_FONT_FAMILY,
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      axis.ticks = ggplot2::element_line(
        color = TEXT_COLOR,
        linewidth = AXIS_LINE_WIDTH
      ),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

has_dotplot_display_terms <- function(result_table) {
  # GseaVis::dotplotGsea会先按p值阈值筛选通路。
  # 若筛选后没有任何通路，ggplot分面在导出PDF/PNG时会报错；
  # 这里提前判断并改为输出空白占位图，不改变GSEA官方结果表。
  if (nrow(result_table) == 0) {
    return(FALSE)
  }

  display_table <- result_table

  if (!is.null(GSEAVIS_DOTPLOT_PARAMS$pval) && "pvalue" %in% colnames(display_table)) {
    display_table <- display_table[
      !is.na(display_table$pvalue) &
        display_table$pvalue <= GSEAVIS_DOTPLOT_PARAMS$pval,
      ,
      drop = FALSE
    ]
  }

  if (!is.null(GSEAVIS_DOTPLOT_PARAMS$pajust) && "p.adjust" %in% colnames(display_table)) {
    display_table <- display_table[
      !is.na(display_table$p.adjust) &
        display_table$p.adjust <= GSEAVIS_DOTPLOT_PARAMS$pajust,
      ,
      drop = FALSE
    ]
  }

  nrow(display_table) > 0
}

make_gsea_dotplot <- function(gsea_result, result_table, analysis_name, geneset_name) {
  if (!has_dotplot_display_terms(result_table)) {
    return(list(
      plot = make_empty_gsea_plot(),
      shown_terms = 0L,
      plot_labels = character(0)
    ))
  }

  dotplot_object <- with_gsea_warnings_suppressed(
    do.call(
      GseaVis::dotplotGsea,
      c(list(data = gsea_result), GSEAVIS_DOTPLOT_PARAMS)
    )
  )

  if (is.list(dotplot_object) && "plot" %in% names(dotplot_object)) {
    plot_object <- dotplot_object$plot
    shown_terms <- if ("df" %in% names(dotplot_object)) {
      nrow(dotplot_object$df)
    } else {
      min(nrow(result_table), GSEAVIS_DOTPLOT_PARAMS$topn * 2L)
    }
  } else {
    plot_object <- dotplot_object
    shown_terms <- min(nrow(result_table), GSEAVIS_DOTPLOT_PARAMS$topn * 2L)
  }

  if (is.list(dotplot_object) && "df" %in% names(dotplot_object) &&
      nrow(dotplot_object$df) == 0) {
    return(list(
      plot = make_empty_gsea_plot(),
      shown_terms = 0L,
      plot_labels = character(0)
    ))
  }

  plot_labels <- if (is.list(dotplot_object) && "df" %in% names(dotplot_object)) {
    levels(dotplot_object$df$Description)
  } else {
    character(0)
  }

  list(
    plot = apply_common_gsea_theme(plot_object + ggplot2::labs(title = NULL)),
    shown_terms = shown_terms,
    plot_labels = plot_labels
  )
}

gsea_scores_compat <- function(gene_list, gene_set, exponent = 1, fortify = TRUE) {
  # 兼容GseaVis旧接口需要的DOSE::gseaScores。
  # 该函数按经典加权GSEA公式计算running enrichment score，只用于绘图，不改变GSEA结果。
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_set <- intersect(as.character(gene_set), names(gene_list))

  hit_index <- names(gene_list) %in% gene_set
  hit_count <- sum(hit_index)
  gene_count <- length(gene_list)

  if (hit_count == 0 || hit_count == gene_count) {
    stop("Invalid gene set for GseaVis plotting.")
  }

  rank_weight <- abs(gene_list)^exponent
  hit_score <- cumsum(ifelse(
    hit_index,
    rank_weight / sum(rank_weight[hit_index]),
    0
  ))
  miss_score <- cumsum(ifelse(!hit_index, 1 / (gene_count - hit_count), 0))
  running_score <- hit_score - miss_score

  data.frame(
    x = seq_along(gene_list),
    runningScore = running_score,
    position = as.integer(hit_index)
  )
}

patch_gseavis_gsinfo_if_needed <- function() {
  # 当前GseaVis版本的gsInfo会调用旧版DOSE内部函数gseaScores；
  # 若本机DOSE已移除该函数，则在当前R会话内替换GseaVis::gsInfo以保证gseaNb可用。
  if (exists("gseaScores", asNamespace("DOSE"), inherits = FALSE)) {
    return(invisible(FALSE))
  }

  gsinfo_compat <- function(object, geneSetID) {
    geneList <- object@geneList
    if (is.numeric(geneSetID)) {
      geneSetID <- object@result[geneSetID, "ID"]
    }

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gsea_scores_compat(geneList, geneSet, exponent, fortify = TRUE)
    df$ymin <- 0
    df$ymax <- 0

    pos <- df$position == 1
    h <- diff(range(df$runningScore)) / 20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    df$Description <- object@result[geneSetID, "Description"]
    df
  }

  gseavis_namespace <- asNamespace("GseaVis")
  was_locked <- bindingIsLocked("gsInfo", gseavis_namespace)
  if (was_locked) {
    unlockBinding("gsInfo", gseavis_namespace)
  }
  assign("gsInfo", gsinfo_compat, envir = gseavis_namespace)
  if (was_locked) {
    lockBinding("gsInfo", gseavis_namespace)
  }
  invisible(TRUE)
}

get_analysis_design_row <- function(analysis_designs, analysis_name) {
  design_index <- match(analysis_name, analysis_designs$Analysis_Name)
  if (is.na(design_index)) {
    stop("Cannot find analysis design for ", analysis_name)
  }

  analysis_designs[design_index, , drop = FALSE]
}

get_single_pathway_sample_info <- function(
    analysis_designs,
    sample_info_all,
    analysis_name) {
  design_row <- get_analysis_design_row(analysis_designs, analysis_name)
  design_samples <- prepare_design_samples(
    sample_info = sample_info_all,
    group_column_index = design_row$Column_Index,
    experiment_group = design_row$Experiment_Group
  )

  sample_info <- design_samples$sample_info
  group_value <- as.character(design_samples$group_list)

  # 表达热图中优先展示实验组(LRC)，再展示对照组(BULK)，便于解读方向。
  sample_order <- order(
    group_value != design_row$Experiment_Group,
    group_value,
    sample_info$Sample_ID
  )
  sample_info <- sample_info[sample_order, , drop = FALSE]

  label_column <- SINGLE_PATHWAY_SAMPLE_LABEL_COLUMN
  if (!label_column %in% colnames(sample_info)) {
    label_column <- "Sample_ID"
  }

  sample_labels <- trimws(as.character(sample_info[[label_column]]))
  empty_label <- is.na(sample_labels) | sample_labels == ""
  sample_labels[empty_label] <- sample_info$Sample_ID[empty_label]
  sample_labels <- make.unique(sample_labels, sep = "_")

  list(
    sample_info = sample_info,
    sample_labels = sample_labels
  )
}

get_expression_cache_file <- function(analysis_name, sample_ids, sample_labels) {
  hash_value <- cache_hash(
    "single_pathway_expression",
    analysis_name,
    get_file_stamp(SE_RDS_FILE),
    TPM_ASSAY_NAME,
    SINGLE_PATHWAY_SAMPLE_LABEL_COLUMN,
    paste(sample_ids, collapse = ","),
    paste(sample_labels, collapse = ",")
  )

  file.path(
    QS2_CACHE_DIR,
    "single_pathway_expression",
    paste0(sanitize_file_name(analysis_name), "_", hash_value, ".qs2")
  )
}

prepare_single_pathway_expression_table <- function(
    analysis_name,
    expression_matrix_all,
    gene_annotation_all,
    analysis_designs,
    sample_info_all) {
  sample_data <- get_single_pathway_sample_info(
    analysis_designs = analysis_designs,
    sample_info_all = sample_info_all,
    analysis_name = analysis_name
  )

  sample_ids <- sample_data$sample_info$Sample_ID
  sample_labels <- sample_data$sample_labels
  cache_file <- get_expression_cache_file(analysis_name, sample_ids, sample_labels)
  cached_expression <- read_qs2_cache(cache_file)
  if (!is.null(cached_expression)) {
    return(cached_expression)
  }

  stopifnot(all(sample_ids %in% colnames(expression_matrix_all)))
  stopifnot("Symbol" %in% colnames(gene_annotation_all))

  expr_matrix <- expression_matrix_all[, sample_ids, drop = FALSE]
  colnames(expr_matrix) <- sample_labels

  gene_symbol <- trimws(as.character(gene_annotation_all$Symbol))
  keep_gene <- !is.na(gene_symbol) & gene_symbol != ""
  expr_matrix <- expr_matrix[keep_gene, , drop = FALSE]
  gene_symbol <- gene_symbol[keep_gene]

  # 同名Symbol只保留平均表达最高的一条记录，避免GseaVis热图中重复基因名。
  mean_expr <- rowMeans(expr_matrix, na.rm = TRUE)
  keep_order <- order(-mean_expr, gene_symbol)
  expr_matrix <- expr_matrix[keep_order, , drop = FALSE]
  gene_symbol <- gene_symbol[keep_order]
  unique_index <- !duplicated(gene_symbol)
  expr_matrix <- expr_matrix[unique_index, , drop = FALSE]
  gene_symbol <- gene_symbol[unique_index]

  expression_table <- data.frame(
    gene_name = gene_symbol,
    expr_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expression_data <- list(
    expression_table = expression_table,
    sample_order = sample_labels,
    sample_info = sample_data$sample_info
  )

  write_qs2_cache(expression_data, cache_file)
  expression_data
}

select_single_pathway_ids <- function(result_table) {
  if (!DRAW_SINGLE_PATHWAY_GSEA || nrow(result_table) == 0) {
    return(character(0))
  }

  stopifnot(SINGLE_PATHWAY_PVALUE_COLUMN %in% colnames(result_table))
  plot_table <- result_table[
    is.finite(as.numeric(result_table[[SINGLE_PATHWAY_PVALUE_COLUMN]])) &
      as.numeric(result_table[[SINGLE_PATHWAY_PVALUE_COLUMN]]) <=
        SINGLE_PATHWAY_PVALUE_CUTOFF,
    ,
    drop = FALSE
  ]

  if (nrow(plot_table) == 0) {
    return(character(0))
  }

  if (length(SINGLE_PATHWAY_KEYWORDS) > 0) {
    search_text <- paste(plot_table$ID, plot_table$Description, sep = " ")
    keyword_pattern <- paste(SINGLE_PATHWAY_KEYWORDS, collapse = "|")
    matched_index <- grepl(keyword_pattern, search_text, ignore.case = TRUE)
    plot_table <- plot_table[matched_index, , drop = FALSE]

    if (nrow(plot_table) == 0) {
      return(character(0))
    }

    plot_table <- plot_table[order(
      as.numeric(plot_table[[SINGLE_PATHWAY_PVALUE_COLUMN]]),
      -abs(as.numeric(plot_table$NES))
    ), , drop = FALSE]
    return(head(unique(plot_table$ID), SINGLE_PATHWAY_MAX_KEYWORD_TERMS))
  }

  plot_table <- plot_table[order(
    as.numeric(plot_table[[SINGLE_PATHWAY_PVALUE_COLUMN]]),
    -abs(as.numeric(plot_table$NES))
  ), , drop = FALSE]
  head(unique(plot_table$ID), SINGLE_PATHWAY_TOP_N)
}

get_core_gene_count <- function(result_table, term_id) {
  term_index <- match(term_id, result_table$ID)
  if (is.na(term_index) || !"core_enrichment" %in% colnames(result_table)) {
    return(20L)
  }

  core_genes <- unique(unlist(strsplit(result_table$core_enrichment[term_index], "/")))
  length(core_genes[core_genes != "" & !is.na(core_genes)])
}

get_single_pathway_plot_size <- function(result_table, term_id, sample_count) {
  core_gene_count <- get_core_gene_count(result_table, term_id)
  term_index <- match(term_id, result_table$ID)
  title_text <- if (is.na(term_index)) term_id else result_table$Description[term_index]
  title_text <- format_gsea_description_for_plot(title_text)
  title_lines <- length(strsplit(
    paste(strwrap(title_text, width = GSEAVIS_SINGLE_PATHWAY_PARAMS$termWidth), collapse = "\n"),
    "\n",
    fixed = TRUE
  )[[1]])

  plot_size <- SINGLE_PATHWAY_PLOT_BASE_SIZE +
    core_gene_count * SINGLE_PATHWAY_PLOT_SIZE_PER_GENE +
    sample_count * SINGLE_PATHWAY_PLOT_SIZE_PER_SAMPLE +
    title_lines * SINGLE_PATHWAY_PLOT_SIZE_PER_TITLE_LINE

  # 单通路图要求整体保持正方形。这里额外根据热图样本名和核心基因名数量
  # 估算最低画布大小，避免样本名或基因名在GseaVis热图区域发生重叠。
  heatmap_ratio <- GSEAVIS_SINGLE_PATHWAY_PARAMS$ght.relHight
  sample_required_size <- sample_count *
    SINGLE_PATHWAY_SAMPLE_TEXT_MIN_SIZE / 72 /
    max(heatmap_ratio * SINGLE_PATHWAY_SAMPLE_TEXT_HEIGHT_FACTOR, 0.05)
  gene_required_size <- core_gene_count *
    SINGLE_PATHWAY_GENE_TEXT_MIN_SIZE / 72 /
    max(SINGLE_PATHWAY_GENE_TEXT_WIDTH_FACTOR, 0.05)

  plot_size <- max(plot_size, sample_required_size, gene_required_size)
  plot_size <- min(
    max(plot_size, SINGLE_PATHWAY_PLOT_MIN_SIZE),
    SINGLE_PATHWAY_PLOT_MAX_SIZE
  )

  # 主图区域按正方形导出，避免上方GSEA曲线和下方表达热图比例失衡。
  list(width = plot_size, height = plot_size)
}

get_single_pathway_gene_text_size <- function(plot_width, core_gene_count) {
  # GseaVis表达热图中基因名为90度竖排；
  # 这里按图宽和核心基因数动态缩小字号，尽量保证相邻基因名不重叠。
  core_gene_count <- max(as.integer(core_gene_count), 1L)
  available_points <- plot_width * 72 * SINGLE_PATHWAY_GENE_TEXT_WIDTH_FACTOR /
    core_gene_count

  min(
    max(available_points, SINGLE_PATHWAY_GENE_TEXT_MIN_SIZE),
    SINGLE_PATHWAY_GENE_TEXT_MAX_SIZE
  )
}

get_single_pathway_sample_text_size <- function(plot_height, sample_count) {
  # GseaVis没有单独暴露表达热图样本名字号参数；
  # 这里通过patchwork全局主题设置axis.text.y，并按热图高度动态计算字号。
  sample_count <- max(as.integer(sample_count), 1L)
  heatmap_height <- plot_height * GSEAVIS_SINGLE_PATHWAY_PARAMS$ght.relHight
  available_points <- heatmap_height * 72 *
    SINGLE_PATHWAY_SAMPLE_TEXT_HEIGHT_FACTOR / sample_count

  min(
    max(available_points, SINGLE_PATHWAY_SAMPLE_TEXT_MIN_SIZE),
    SINGLE_PATHWAY_SAMPLE_TEXT_MAX_SIZE
  )
}

apply_common_single_gsea_theme <- function(plot, sample_text_size) {
  common_theme <- ggplot2::theme(
    text = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    plot.title = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR,
      hjust = 0.5
    ),
    axis.text = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    axis.text.y = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR,
      size = sample_text_size
    ),
    axis.title = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    legend.title = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    legend.text = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    panel.border = ggplot2::element_rect(
      color = TEXT_COLOR,
      linewidth = AXIS_LINE_WIDTH,
      fill = NA
    )
  )

  # gseaNb返回patchwork对象时，&可将主题应用到所有子图；若失败则回退到普通ggplot叠加。
  tryCatch(
    plot & common_theme,
    error = function(e) plot + common_theme
  )
}

make_single_pathway_gsea_plot <- function(
    gsea_result,
    term_id,
    expression_data,
    result_table,
    plot_size) {
  params <- GSEAVIS_SINGLE_PATHWAY_PARAMS
  params$object <- gsea_result
  params$geneSetID <- term_id
  params$ght.geneText.size <- get_single_pathway_gene_text_size(
    plot_width = plot_size$width,
    core_gene_count = get_core_gene_count(result_table, term_id)
  )
  sample_text_size <- get_single_pathway_sample_text_size(
    plot_height = plot_size$height,
    sample_count = length(expression_data$sample_order)
  )

  if (isTRUE(params$add.geneExpHt)) {
    params$exp <- expression_data$expression_table
    params$sample.order <- expression_data$sample_order
  }

  # setReadable后core_enrichment是Symbol，但geneList仍保留原始Entrez。
  # GseaVis在kegg=TRUE时会使用gseaResult@gene2Symbol匹配表达热图，适合当前readable结果。
  params$kegg <- READABLE_GENE_SYMBOLS && GENE_ID_TYPE != "SYMBOL"

  plot <- with_gsea_warnings_suppressed(
    do.call(GseaVis::gseaNb, params)
  )

  apply_common_single_gsea_theme(plot, sample_text_size = sample_text_size)
}

save_single_pathway_gsea_plots <- function(
    gsea_result,
    result_table,
    analysis_name,
    geneset_name,
    plot_output_dir,
    expression_data) {
  selected_term_ids <- select_single_pathway_ids(result_table)
  if (length(selected_term_ids) == 0) {
    return(0L)
  }

  patch_gseavis_gsinfo_if_needed()
  single_pathway_dir <- file.path(plot_output_dir, "single_pathway")
  dir.create(single_pathway_dir, recursive = TRUE, showWarnings = FALSE)

  save_one_single_pathway_plot <- function(i) {
    patch_gseavis_gsinfo_if_needed()
    term_id <- selected_term_ids[i]
    term_dir <- file.path(
      single_pathway_dir,
      paste0(sprintf("%02d_", i), sanitize_file_name(term_id))
    )
    dir.create(term_dir, recursive = TRUE, showWarnings = FALSE)

    plot_size <- get_single_pathway_plot_size(
      result_table = result_table,
      term_id = term_id,
      sample_count = length(expression_data$sample_order)
    )
    single_plot <- make_single_pathway_gsea_plot(
      gsea_result = gsea_result,
      term_id = term_id,
      expression_data = expression_data,
      result_table = result_table,
      plot_size = plot_size
    )

    save_ggplot_pdf_png(
      plot = single_plot,
      pdf_file = file.path(term_dir, "gsea_plot.pdf"),
      width = plot_size$width,
      height = plot_size$height
    )

    TRUE
  }

  plot_index <- seq_along(selected_term_ids)
  if (.Platform$OS.type != "windows" &&
      SINGLE_PATHWAY_PLOT_WORKERS > 1 &&
      length(plot_index) > 1) {
    saved_flags <- parallel::mclapply(
      plot_index,
      save_one_single_pathway_plot,
      mc.cores = min(SINGLE_PATHWAY_PLOT_WORKERS, length(plot_index)),
      mc.preschedule = FALSE
    )
  } else {
    saved_flags <- lapply(plot_index, save_one_single_pathway_plot)
  }

  sum(vapply(saved_flags, isTRUE, logical(1)))
}
