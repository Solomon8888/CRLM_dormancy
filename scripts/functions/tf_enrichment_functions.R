# 转录因子富集分析公共函数
#
# 本文件只保存08号脚本需要复用的输入解析、TF方法运行和结果保存函数。
# 主脚本负责配置与并行调度；这里的函数尽量保持方法学含义清晰，避免把主流程写得过长。


sanitize_tf_path_name <- function(x) {
  # 将分析名或交集方案名转成稳定目录名。
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

is_all_keyword <- function(x) {
  length(x) == 1L && tolower(as.character(x)) == "all"
}

get_runtime_vector <- function(env_name, default_value) {
  # 支持临时环境变量覆盖配置，主要用于小规模测试，不影响常规交互式运行。
  env_value <- Sys.getenv(env_name, unset = "")
  if (env_value == "") {
    return(default_value)
  }

  trimws(strsplit(env_value, ",", fixed = TRUE)[[1]])
}

read_result_table <- function(file_name) {
  read.csv(
    file_name,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

get_analysis_name_from_significant_file <- function(file_name) {
  # significant_genes.csv位于tables/<analysis_name>/DEG/。
  basename(dirname(dirname(file_name)))
}

get_tf_deg_inputs <- function(table_root, analyses_to_run = "all") {
  significant_files <- list.files(
    table_root,
    pattern = "significant_genes[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  significant_files <- significant_files[
    basename(dirname(significant_files)) == "DEG"
  ]

  if (length(significant_files) == 0L) {
    stop("No significant_genes.csv files were found under: ", table_root)
  }

  analysis_names <- vapply(
    significant_files,
    get_analysis_name_from_significant_file,
    character(1)
  )

  keep <- rep(TRUE, length(significant_files))
  if (!is_all_keyword(analyses_to_run)) {
    keep <- analysis_names %in% analyses_to_run
  }

  significant_files <- significant_files[keep]
  analysis_names <- analysis_names[keep]

  if (length(significant_files) == 0L) {
    stop("No DEG inputs matched ANALYSES_TO_RUN.")
  }

  inputs <- Map(function(analysis_name, significant_file) {
    list(
      input_type = "DEG",
      input_name = analysis_name,
      gene_file = significant_file,
      rank_files = setNames(significant_file, analysis_name),
      file_prefix = "deg",
      output_label = analysis_name
    )
  }, analysis_names, significant_files)

  names(inputs) <- analysis_names
  inputs[order(names(inputs))]
}

get_tf_intersection_inputs <- function(intersect_root, schemes_to_run = "all") {
  if (!dir.exists(intersect_root)) {
    return(list())
  }

  gene_list_files <- list.files(
    intersect_root,
    pattern = "gene_list[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  gene_list_files <- gene_list_files[
    basename(dirname(gene_list_files)) != "DEG"
  ]

  if (length(gene_list_files) == 0L) {
    return(list())
  }

  scheme_names <- basename(dirname(gene_list_files))

  keep <- rep(TRUE, length(gene_list_files))
  if (!is_all_keyword(schemes_to_run)) {
    keep <- scheme_names %in% schemes_to_run
  }

  gene_list_files <- gene_list_files[keep]
  scheme_names <- scheme_names[keep]

  inputs <- Map(function(scheme_name, gene_file) {
    scheme_dir <- dirname(gene_file)
    rank_files <- list.files(
      scheme_dir,
      pattern = "deg_results[.]csv$",
      recursive = TRUE,
      full.names = TRUE
    )

    rank_names <- basename(dirname(rank_files))
    if (length(rank_files) > 0L) {
      rank_files <- setNames(rank_files, rank_names)
    }

    list(
      input_type = "INTERSECT",
      input_name = scheme_name,
      gene_file = gene_file,
      rank_files = rank_files,
      file_prefix = "intersect",
      output_label = scheme_name
    )
  }, scheme_names, gene_list_files)

  names(inputs) <- scheme_names
  inputs[order(names(inputs))]
}

extract_gene_symbols <- function(dat, symbol_column = "Symbol") {
  if (!symbol_column %in% colnames(dat)) {
    stop("Missing symbol column: ", symbol_column)
  }

  symbols <- toupper(trimws(as.character(dat[[symbol_column]])))
  symbols <- symbols[!is.na(symbols) & symbols != "" & symbols != "NA"]
  sort(unique(symbols))
}

extract_ranked_signature <- function(
    dat,
    symbol_column = "Symbol",
    rank_metric_column = "t",
    fallback_rank_columns = c("t", "logFC")) {
  rank_column <- rank_metric_column
  if (!rank_column %in% colnames(dat)) {
    rank_column <- fallback_rank_columns[fallback_rank_columns %in% colnames(dat)][1]
  }
  if (is.na(rank_column) || !rank_column %in% colnames(dat)) {
    stop(
      "No usable rank metric column was found. Tried: ",
      paste(unique(c(rank_metric_column, fallback_rank_columns)), collapse = ", ")
    )
  }
  if (!symbol_column %in% colnames(dat)) {
    stop("Missing symbol column: ", symbol_column)
  }

  symbols <- toupper(trimws(as.character(dat[[symbol_column]])))
  values <- suppressWarnings(as.numeric(dat[[rank_column]]))
  keep <- !is.na(symbols) & symbols != "" & symbols != "NA" & !is.na(values)
  symbols <- symbols[keep]
  values <- values[keep]

  if (length(values) == 0L) {
    return(numeric(0))
  }

  # 同一Symbol可能由多个转录本/Ensembl ID对应；保留绝对统计量最大的记录。
  split_values <- split(values, symbols)
  collapsed <- vapply(split_values, function(x) {
    x[which.max(abs(x))]
  }, numeric(1))

  collapsed[order(abs(collapsed), decreasing = TRUE)]
}

make_binary_gene_matrix <- function(
    genes,
    network_targets,
    condition_name = "gene_set") {
  genes <- sort(unique(toupper(genes)))
  background <- sort(unique(toupper(c(network_targets, genes))))
  mat <- matrix(
    0,
    nrow = length(background),
    ncol = 1,
    dimnames = list(background, condition_name)
  )
  mat[intersect(genes, background), 1] <- 1
  mat
}

make_rank_matrix <- function(
    rank_files,
    symbol_column = "Symbol",
    rank_metric_column = "t",
    fallback_rank_columns = c("t", "logFC")) {
  if (length(rank_files) == 0L) {
    return(matrix(nrow = 0, ncol = 0))
  }

  ranked_vectors <- lapply(rank_files, function(file_name) {
    extract_ranked_signature(
      dat = read_result_table(file_name),
      symbol_column = symbol_column,
      rank_metric_column = rank_metric_column,
      fallback_rank_columns = fallback_rank_columns
    )
  })

  ranked_vectors <- ranked_vectors[vapply(ranked_vectors, length, integer(1)) > 0L]
  if (length(ranked_vectors) == 0L) {
    return(matrix(nrow = 0, ncol = 0))
  }

  all_genes <- sort(unique(unlist(lapply(ranked_vectors, names), use.names = FALSE)))
  mat <- matrix(
    0,
    nrow = length(all_genes),
    ncol = length(ranked_vectors),
    dimnames = list(all_genes, names(ranked_vectors))
  )

  for (condition_name in names(ranked_vectors)) {
    x <- ranked_vectors[[condition_name]]
    mat[names(x), condition_name] <- as.numeric(x)
  }

  mat
}

load_dorothea_regulon <- function(
    species = "human",
    confidence_levels = c("A", "B", "C")) {
  # DoRothEA官方Bioconductor包目前直接提供human/mouse regulon数据表。
  species_key <- tolower(species)
  if (species_key %in% c("human", "homo sapiens", "hs")) {
    data("dorothea_hs", package = "dorothea")
    network <- dorothea_hs
  } else if (species_key %in% c("mouse", "mus musculus", "mm")) {
    data("dorothea_mm", package = "dorothea")
    network <- dorothea_mm
  } else {
    stop(
      "DoRothEA package currently supports human/mouse regulons in this script. ",
      "Please use SPECIES <- \"human\" or SPECIES <- \"mouse\"."
    )
  }

  network <- network[network$confidence %in% confidence_levels, , drop = FALSE]
  network$tf <- toupper(as.character(network$tf))
  network$target <- toupper(as.character(network$target))
  network
}

get_ncbi_tax_id <- function(species = "human") {
  species_key <- tolower(species)
  if (species_key %in% c("human", "homo sapiens", "hs", "9606")) {
    return(9606L)
  }
  if (species_key %in% c("mouse", "mus musculus", "mm", "10090")) {
    return(10090L)
  }
  if (species_key %in% c("rat", "rattus norvegicus", "rn", "10116")) {
    return(10116L)
  }
  stop("Unsupported species: ", species)
}

format_network_for_decoupler <- function(
    dat,
    source_column,
    target_column,
    mor_column = NULL,
    likelihood_column = NULL) {
  network <- data.frame(
    source = toupper(trimws(as.character(dat[[source_column]]))),
    target = toupper(trimws(as.character(dat[[target_column]]))),
    stringsAsFactors = FALSE
  )

  if (!is.null(mor_column) && mor_column %in% colnames(dat)) {
    network$mor <- suppressWarnings(as.numeric(dat[[mor_column]]))
  } else {
    network$mor <- 1
  }

  network$mor[is.na(network$mor)] <- 0
  network$mor[network$mor > 0] <- 1
  network$mor[network$mor < 0] <- -1

  if (!is.null(likelihood_column) && likelihood_column %in% colnames(dat)) {
    network$likelihood <- suppressWarnings(as.numeric(dat[[likelihood_column]]))
  } else {
    network$likelihood <- 1
  }
  network$likelihood[is.na(network$likelihood)] <- 1

  network <- network[
    !is.na(network$source) & !is.na(network$target) &
      network$source != "" & network$target != "",
    ,
    drop = FALSE
  ]
  collapse_network_edges(network)
}

collapse_network_edges <- function(network) {
  # decoupleR要求source-target边唯一；多个证据支持同一边时，折叠为一条边。
  if (nrow(network) == 0L) {
    return(network)
  }

  edge_key <- paste(network$source, network$target, sep = "\t")
  split_rows <- split(seq_len(nrow(network)), edge_key)

  collapsed <- lapply(split_rows, function(idx) {
    edge <- network[idx[1], , drop = FALSE]
    mor_values <- network$mor[idx]
    mor_values <- mor_values[!is.na(mor_values)]
    if (length(mor_values) > 0L) {
      non_zero_mor <- mor_values[mor_values != 0]
      edge$mor <- if (length(non_zero_mor) > 0L) {
        sign(sum(non_zero_mor))
      } else {
        0
      }
    }

    likelihood_values <- network$likelihood[idx]
    likelihood_values <- likelihood_values[!is.na(likelihood_values)]
    if (length(likelihood_values) > 0L) {
      edge$likelihood <- max(likelihood_values)
    }

    edge
  })

  do.call(rbind, collapsed)
}

load_trrust_network <- function(species = "human") {
  # TRRUST v2官方资源通过OmnipathR::trrust_download获取。
  species_key <- tolower(species)
  organism <- if (species_key %in% c("human", "homo sapiens", "hs")) {
    "human"
  } else if (species_key %in% c("mouse", "mus musculus", "mm")) {
    "mouse"
  } else {
    stop("TRRUST currently supports human/mouse in OmnipathR.")
  }

  trrust <- OmnipathR::trrust_download(organism = organism)
  network <- format_network_for_decoupler(
    dat = trrust,
    source_column = "source_genesymbol",
    target_column = "target_genesymbol",
    mor_column = "effect"
  )
  network$reference <- trrust$reference[
    match(
      paste(network$source, network$target),
      paste(
        toupper(as.character(trrust$source_genesymbol)),
        toupper(as.character(trrust$target_genesymbol))
      )
    )
  ]
  network
}

read_collectri_static_table <- function(species = "human") {
  # 当前OmnipathR版本在解析CollecTRI evidences时可能报错；
  # 因此这里使用OmnipathR官方static_tables索引定位官方静态TSV，再由R直接解析。
  tax_id <- as.character(get_ncbi_tax_id(species))
  static_tables <- OmnipathR::static_tables()
  collectri_row <- static_tables[
    static_tables$query == "interactions" &
      static_tables$resource == "collectri" &
      static_tables$organism == tax_id,
    ,
    drop = FALSE
  ]

  if (nrow(collectri_row) == 0L) {
    stop("No CollecTRI static table was found for organism tax_id: ", tax_id)
  }

  tmp_file <- tempfile(fileext = ".tsv.gz")
  utils::download.file(
    collectri_row$url[1],
    tmp_file,
    quiet = TRUE,
    mode = "wb"
  )
  read.delim(
    gzfile(tmp_file),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

load_collectri_network <- function(species = "human", split_complexes = FALSE) {
  # 优先使用decoupleR官方接口；若当前OmnipathR解析失败，则使用官方静态表兜底。
  collectri <- try(
    decoupleR::get_collectri(
      organism = species,
      split_complexes = split_complexes,
      load_meta = FALSE
    ),
    silent = TRUE
  )

  if (inherits(collectri, "try-error")) {
    collectri <- read_collectri_static_table(species = species)
  }

  if (all(c("source", "target", "mor") %in% colnames(collectri))) {
    network <- format_network_for_decoupler(
      dat = collectri,
      source_column = "source",
      target_column = "target",
      mor_column = "mor"
    )
  } else {
    mor_column <- if ("consensus_direction" %in% colnames(collectri)) {
      "consensus_direction"
    } else {
      NULL
    }
    likelihood_column <- if ("curation_effort" %in% colnames(collectri)) {
      "curation_effort"
    } else {
      NULL
    }
    network <- format_network_for_decoupler(
      dat = collectri,
      source_column = "source_genesymbol",
      target_column = "target_genesymbol",
      mor_column = mor_column,
      likelihood_column = likelihood_column
    )
  }

  network
}

run_dorothea_ora_enrichment <- function(
    input,
    network,
    symbol_column = "Symbol",
    params = list()) {
  dat <- read_result_table(input$gene_file)
  genes <- extract_gene_symbols(dat, symbol_column = symbol_column)
  if (length(genes) == 0L) {
    return(data.frame())
  }

  mat <- make_binary_gene_matrix(
    genes = genes,
    network_targets = network$target,
    condition_name = input$input_name
  )

  args <- c(
    list(
      mat = mat,
      network = network,
      .source = quote(tf),
      .target = quote(target),
      n_up = length(genes),
      n_bottom = 0
    ),
    params
  )

  as.data.frame(do.call(decoupleR::run_ora, args))
}

run_viper_activity_enrichment <- function(
    input,
    network,
    symbol_column = "Symbol",
    rank_metric_column = "t",
    fallback_rank_columns = c("t", "logFC"),
    params = list()) {
  if (length(input$rank_files) == 0L) {
    return(data.frame())
  }

  regulon <- dorothea::df2regulon(network)
  rank_conditions <- names(input$rank_files)
  result_matrices <- lapply(rank_conditions, function(condition_name) {
    signature <- extract_ranked_signature(
      dat = read_result_table(input$rank_files[[condition_name]]),
      symbol_column = symbol_column,
      rank_metric_column = rank_metric_column,
      fallback_rank_columns = fallback_rank_columns
    )

    if (length(signature) == 0L) {
      return(NULL)
    }

    args <- c(
      list(
        eset = signature,
        regulon = regulon
      ),
      params
    )
    result <- as.matrix(do.call(viper::viper, args))
    colnames(result) <- condition_name
    result
  })
  names(result_matrices) <- rank_conditions

  result_matrices <- result_matrices[!vapply(result_matrices, is.null, logical(1))]
  if (length(result_matrices) == 0L) {
    return(data.frame())
  }

  condition_names <- vapply(
    result_matrices,
    function(x) colnames(x)[1],
    character(1)
  )

  all_tfs <- sort(unique(unlist(lapply(result_matrices, rownames), use.names = FALSE)))
  output_matrix <- matrix(
    NA_real_,
    nrow = length(all_tfs),
    ncol = length(result_matrices),
    dimnames = list(all_tfs, condition_names)
  )

  for (condition_index in seq_along(result_matrices)) {
    condition_name <- colnames(result_matrices[[condition_index]])[1]
    output_matrix[rownames(result_matrices[[condition_index]]), condition_name] <-
      result_matrices[[condition_index]][, 1]
  }

  data.frame(TF = rownames(output_matrix), output_matrix, check.names = FALSE)
}

run_chea3_enrichment <- function(
    input,
    symbol_column = "Symbol",
    min_genes = 5,
    api_url = "https://maayanlab.cloud/chea3/api/enrich/") {
  dat <- read_result_table(input$gene_file)
  genes <- extract_gene_symbols(dat, symbol_column = symbol_column)
  if (length(genes) < min_genes) {
    return(list())
  }

  rChEA3::queryChEA3(
    genes = genes,
    query_name = paste(input$input_type, input$input_name, sep = "_"),
    verbose = FALSE,
    url = api_url
  )
}

run_enrichr_enrichment <- function(
    input,
    databases,
    symbol_column = "Symbol",
    min_genes = 5,
    background = NULL,
    include_overlap = TRUE,
    sleep_time = 1,
    enrichr_site = "Enrichr") {
  dat <- read_result_table(input$gene_file)
  genes <- extract_gene_symbols(dat, symbol_column = symbol_column)
  if (length(genes) < min_genes) {
    return(list())
  }

  # enrichR的站点配置依赖包加载时初始化的options；显式attach可避免namespace调用时缺少option。
  suppressPackageStartupMessages(
    require("enrichR", character.only = TRUE)
  )
  if (!is.null(enrichr_site) && enrichr_site != "" && enrichr_site != "Enrichr") {
    enrichR::setEnrichrSite(enrichr_site)
  }
  enrichR::enrichr(
    genes = genes,
    databases = databases,
    background = background,
    include_overlap = include_overlap,
    sleepTime = sleep_time
  )
}

run_network_ora_enrichment <- function(
    input,
    network,
    symbol_column = "Symbol",
    params = list()) {
  dat <- read_result_table(input$gene_file)
  genes <- extract_gene_symbols(dat, symbol_column = symbol_column)
  if (length(genes) == 0L) {
    return(data.frame())
  }

  mat <- make_binary_gene_matrix(
    genes = genes,
    network_targets = network$target,
    condition_name = input$input_name
  )

  args <- c(
    list(
      mat = mat,
      network = network,
      .source = quote(source),
      .target = quote(target),
      n_up = length(genes),
      n_bottom = 0
    ),
    params
  )

  as.data.frame(do.call(decoupleR::run_ora, args))
}

run_network_viper_enrichment <- function(
    input,
    network,
    symbol_column = "Symbol",
    rank_metric_column = "t",
    fallback_rank_columns = c("t", "logFC"),
    params = list()) {
  rank_matrix <- make_rank_matrix(
    rank_files = input$rank_files,
    symbol_column = symbol_column,
    rank_metric_column = rank_metric_column,
    fallback_rank_columns = fallback_rank_columns
  )

  if (nrow(rank_matrix) == 0L || ncol(rank_matrix) == 0L) {
    return(data.frame())
  }

  args <- c(
    list(
      mat = rank_matrix,
      network = network,
      .source = quote(source),
      .target = quote(target),
      .mor = quote(mor)
    ),
    params
  )

  as.data.frame(do.call(decoupleR::run_viper, args))
}

save_tf_table <- function(dat, output_dir, file_stem, preview_rows = 21) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(dat) || nrow(dat) == 0L || ncol(dat) == 0L) {
    dat <- data.frame(Message = "No results returned for this input.")
  }

  csv_file <- file.path(output_dir, paste0(file_stem, ".csv"))
  write_csv_with_report_previews(
    dat = dat,
    csv_file = csv_file,
    n_rows = preview_rows
  )
  csv_file
}

save_chea3_tables <- function(result_list, output_dir, file_prefix, preview_rows = 21) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(result_list) == 0L) {
    return(save_tf_table(
      dat = data.frame(Message = "No ChEA3 results returned for this input."),
      output_dir = output_dir,
      file_stem = paste0(file_prefix, "_no_results"),
      preview_rows = preview_rows
    ))
  }

  output_files <- character(0)
  for (library_name in names(result_list)) {
    dat <- as.data.frame(result_list[[library_name]])
    file_stem <- paste0(file_prefix, "_", sanitize_tf_path_name(library_name))
    output_files <- c(
      output_files,
      save_tf_table(
        dat = dat,
        output_dir = output_dir,
        file_stem = file_stem,
        preview_rows = preview_rows
      )
    )
  }

  output_files
}

save_named_result_tables <- function(result_list, output_dir, file_prefix, preview_rows = 21) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(result_list) == 0L) {
    return(save_tf_table(
      dat = data.frame(Message = "No results returned for this input."),
      output_dir = output_dir,
      file_stem = paste0(file_prefix, "_no_results"),
      preview_rows = preview_rows
    ))
  }

  output_files <- character(0)
  for (result_name in names(result_list)) {
    dat <- as.data.frame(result_list[[result_name]])
    file_stem <- paste0(file_prefix, "_", sanitize_tf_path_name(result_name))
    output_files <- c(
      output_files,
      save_tf_table(
        dat = dat,
        output_dir = output_dir,
        file_stem = file_stem,
        preview_rows = preview_rows
      )
    )
  }

  output_files
}

run_tf_enrichment_task <- function(
    task,
    tf_root,
    dorothea_network,
    trrust_network,
    collectri_network,
    symbol_column,
    rank_metric_column,
    fallback_rank_columns,
    dorothea_ora_params,
    trrust_ora_params,
    viper_params,
    collectri_viper_params,
    chea3_min_genes,
    chea3_api_url,
    enrichr_databases,
    enrichr_min_genes,
    enrichr_background,
    enrichr_include_overlap,
    enrichr_sleep_time,
    enrichr_site,
    preview_rows) {
  method <- task$method
  input <- task$input
  output_dir <- file.path(
    tf_root,
    method,
    sanitize_tf_path_name(input$output_label)
  )

  if (method == "dorothea") {
    result <- run_dorothea_ora_enrichment(
      input = input,
      network = dorothea_network,
      symbol_column = symbol_column,
      params = dorothea_ora_params
    )
    output_files <- save_tf_table(
      dat = result,
      output_dir = output_dir,
      file_stem = paste0("dorothea_", input$file_prefix, "_tf_enrichment"),
      preview_rows = preview_rows
    )
  } else if (method == "viper") {
    result <- run_viper_activity_enrichment(
      input = input,
      network = dorothea_network,
      symbol_column = symbol_column,
      rank_metric_column = rank_metric_column,
      fallback_rank_columns = fallback_rank_columns,
      params = viper_params
    )
    output_files <- save_tf_table(
      dat = result,
      output_dir = output_dir,
      file_stem = paste0("viper_", input$file_prefix, "_tf_activity"),
      preview_rows = preview_rows
    )
  } else if (method == "chea3") {
    result <- run_chea3_enrichment(
      input = input,
      symbol_column = symbol_column,
      min_genes = chea3_min_genes,
      api_url = chea3_api_url
    )
    output_files <- save_chea3_tables(
      result_list = result,
      output_dir = output_dir,
      file_prefix = paste0("chea3_", input$file_prefix),
      preview_rows = preview_rows
    )
  } else if (method == "enrichr") {
    result <- run_enrichr_enrichment(
      input = input,
      databases = enrichr_databases,
      symbol_column = symbol_column,
      min_genes = enrichr_min_genes,
      background = enrichr_background,
      include_overlap = enrichr_include_overlap,
      sleep_time = enrichr_sleep_time,
      enrichr_site = enrichr_site
    )
    output_files <- save_named_result_tables(
      result_list = result,
      output_dir = output_dir,
      file_prefix = paste0("enrichr_", input$file_prefix),
      preview_rows = preview_rows
    )
  } else if (method == "trrust") {
    result <- run_network_ora_enrichment(
      input = input,
      network = trrust_network,
      symbol_column = symbol_column,
      params = trrust_ora_params
    )
    output_files <- save_tf_table(
      dat = result,
      output_dir = output_dir,
      file_stem = paste0("trrust_", input$file_prefix, "_tf_enrichment"),
      preview_rows = preview_rows
    )
  } else if (method == "collectri") {
    result <- run_network_viper_enrichment(
      input = input,
      network = collectri_network,
      symbol_column = symbol_column,
      rank_metric_column = rank_metric_column,
      fallback_rank_columns = fallback_rank_columns,
      params = collectri_viper_params
    )
    output_files <- save_tf_table(
      dat = result,
      output_dir = output_dir,
      file_stem = paste0("collectri_", input$file_prefix, "_tf_activity"),
      preview_rows = preview_rows
    )
  } else {
    stop("Unsupported TF enrichment method: ", method)
  }

  gene_count <- length(extract_gene_symbols(
    read_result_table(input$gene_file),
    symbol_column = symbol_column
  ))

  data.frame(
    Method = method,
    Input_Type = input$input_type,
    Input_Name = input$input_name,
    Gene_Count = gene_count,
    Rank_Conditions = paste(names(input$rank_files), collapse = ";"),
    Output_Directory = output_dir,
    Result_Tables = length(output_files),
    stringsAsFactors = FALSE
  )
}
