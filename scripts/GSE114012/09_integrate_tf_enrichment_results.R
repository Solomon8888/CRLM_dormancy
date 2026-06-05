# GSE114012转录因子富集结果整合与交集分析
#
# 读取08号脚本输出的六类TF富集/活性推断结果：
# DoRothEA、ChEA3、VIPER、ENRICHR、TRRUST、CollecTRI。
#
# 本脚本分两步：
# 1. 在每个DEG或intersect方案内，为每种TF分析方法整理一版method_final结果；
# 2. 按脚本头部配置的多套方法组合，对method_final结果做全量TF交集，
#    输出Top候选TF表和summary表。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TF_RESULT_ROOT <- file.path(RESULT_ROOT, "TF")
TF_SUMMARY_ROOT <- file.path(RESULT_ROOT, "TF_summary")
OBSOLETE_TF_VENN_PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "TF_intersection")

# 需要整合哪些输入类型。
# DEG       = results/ngs/GSE114012/TF/DEG/<analysis_name>/
# INTERSECT = results/ngs/GSE114012/TF/intersect/<intersection_scheme>/
INPUT_TYPES_TO_RUN <- c("DEG", "INTERSECT")

# 需要整合哪些具体方案。设为"all"时自动读取全部。
DEG_ANALYSES_TO_RUN <- "all"
INTERSECTION_ANALYSES_TO_RUN <- "all"

# 每套交集最终输出Top多少个候选TF。
INTERSECTION_TOP_N_TO_REPORT <- 10

# 重跑时清空09号整合结果，避免旧结果残留。
CLEAN_TF_INTEGRATION_OUTPUT <- TRUE

# 六种方法的固定顺序。后续summary和交集表均按此顺序展示。
TF_METHODS <- c("dorothea", "chea3", "viper", "enrichr", "trrust", "collectri")

TF_METHOD_LABELS <- c(
  dorothea = "DoRothEA",
  chea3 = "ChEA3",
  viper = "VIPER",
  enrichr = "ENRICHR",
  trrust = "TRRUST",
  collectri = "CollecTRI"
)

# 多套TF方法交集方案。
# 说明：
# - 前四组是按当前研究设计优先关注的组合；
# - 后续组合按“证据来源/是否有方向性/网络覆盖度”补充，用于比较不同证据逻辑下的候选TF稳定性。
TF_INTERSECTION_SCHEMES <- list(
  ALL_6_METHODS = c("dorothea", "chea3", "viper", "enrichr", "trrust", "collectri"),
  WITHOUT_CHEA3 = c("dorothea", "viper", "enrichr", "trrust", "collectri"),
  DOROTHEA_CHEA3_VIPER = c("dorothea", "chea3", "viper"),
  ENRICHR_TRRUST_COLLECTRI = c("enrichr", "trrust", "collectri"),
  ACTIVITY_METHODS_VIPER_COLLECTRI = c("viper", "collectri"),
  ORA_METHODS_DOROTHEA_TRRUST = c("dorothea", "trrust"),
  API_METHODS_CHEA3_ENRICHR = c("chea3", "enrichr"),
  SIGNED_NETWORK_DOROTHEA_VIPER_COLLECTRI = c("dorothea", "viper", "collectri"),
  LIST_BASED_ORA_API_DOROTHEA_CHEA3_ENRICHR_TRRUST = c("dorothea", "chea3", "enrichr", "trrust"),
  CHIP_LITERATURE_EVIDENCE_CHEA3_ENRICHR_TRRUST = c("chea3", "enrichr", "trrust"),
  BROAD_DATABASE_EVIDENCE_CHEA3_ENRICHR_COLLECTRI = c("chea3", "enrichr", "collectri"),
  CURATED_REGULON_NETWORK_DOROTHEA_TRRUST_COLLECTRI = c("dorothea", "trrust", "collectri")
)

# 运行哪些交集方案。默认全部运行。
INTERSECTION_SCHEMES_TO_RUN <- names(TF_INTERSECTION_SCHEMES)

# ChEA3中用于判断library证据的官方library结果。
# Integrated结果单独记录rank，不计入library证据数量。
CHEA3_EVIDENCE_LIBRARIES <- c(
  "ARCHS4--Coexpression",
  "GTEx--Coexpression",
  "ENCODE--ChIP-seq",
  "ReMap--ChIP-seq",
  "Literature--ChIP-seq",
  "Enrichr--Queries"
)

# 兼容历史参数名；当前R脚本只保存完整CSV。
TABLE_PREVIEW_ROWS <- 21

options(width = 200)


# 1. 加载公共函数 --------------------------------------------------------------

source(REPORT_TABLE_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()


# 2. 基础工具函数 --------------------------------------------------------------

standardize_tf <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "" | is.na(x) | x == "NA"] <- NA_character_
  x
}

sanitize_file_name <- function(x, default = "analysis") {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(x == "" | is.na(x), default, x)
}

to_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

first_existing_column <- function(dat, column_names) {
  column_names[column_names %in% colnames(dat)][1]
}

read_tf_csv <- function(file_name) {
  read.csv(
    resolve_report_csv_file(file_name),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

resolve_tf_result_file <- function(file_name) {
  resolve_report_csv_file(file_name)
}

clean_output_dirs <- function() {
  if (!CLEAN_TF_INTEGRATION_OUTPUT) {
    return(invisible(FALSE))
  }

  unlink(TF_SUMMARY_ROOT, recursive = TRUE, force = TRUE)
  unlink(OBSOLETE_TF_VENN_PLOT_ROOT, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

get_method_prefix <- function(input_type) {
  if (toupper(input_type) == "DEG") {
    return("deg")
  }

  "intersect"
}

get_tf_input_dirs <- function(tf_root, input_types) {
  inputs <- list()

  for (input_type in input_types) {
    dir_name <- if (toupper(input_type) == "DEG") "DEG" else "intersect"
    input_root <- file.path(tf_root, dir_name)
    if (!dir.exists(input_root)) {
      next
    }

    analysis_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
    analysis_names <- basename(analysis_dirs)

    if (toupper(input_type) == "DEG" && !identical(DEG_ANALYSES_TO_RUN, "all")) {
      keep <- analysis_names %in% DEG_ANALYSES_TO_RUN
      analysis_dirs <- analysis_dirs[keep]
      analysis_names <- analysis_names[keep]
    }

    if (toupper(input_type) == "INTERSECT" && !identical(INTERSECTION_ANALYSES_TO_RUN, "all")) {
      keep <- analysis_names %in% INTERSECTION_ANALYSES_TO_RUN
      analysis_dirs <- analysis_dirs[keep]
      analysis_names <- analysis_names[keep]
    }

    for (i in seq_along(analysis_dirs)) {
      inputs[[length(inputs) + 1L]] <- list(
        Input_Type = toupper(input_type),
        Analysis_Name = analysis_names[i],
        Input_Dir = analysis_dirs[i],
        Method_Prefix = get_method_prefix(input_type)
      )
    }
  }

  inputs
}

get_method_csv_files <- function(input_dir, method_name) {
  method_dir <- file.path(input_dir, method_name)
  if (!dir.exists(method_dir)) {
    return(character(0))
  }

  current_files <- list.files(
    method_dir,
    pattern = "[.]csv$",
    full.names = TRUE,
    recursive = FALSE
  )

  legacy_dir <- file.path(method_dir, "csv")
  legacy_files <- if (dir.exists(legacy_dir)) {
    list.files(legacy_dir, pattern = "[.]csv$", full.names = TRUE)
  } else {
    character(0)
  }

  out <- c(current_files, legacy_files)
  prefer_report_csv_files(out, function(x) basename(normalize_report_csv_path(x)))
}

get_library_name_from_file <- function(file_name, method_name, prefix) {
  stem <- tools::file_path_sans_ext(basename(file_name))
  sub(
    paste0("^", method_name, "_", prefix, "_"),
    "",
    stem
  )
}

add_missing_columns <- function(dat, columns) {
  for (column in columns) {
    if (!column %in% colnames(dat)) {
      dat[[column]] <- NA
    }
  }

  dat
}

sort_by_rank_columns <- function(dat) {
  dat <- add_missing_columns(dat, c("Method_Rank", "Method_Score", "Method_P_Value"))
  dat$Method_Rank <- to_number(dat$Method_Rank)
  dat$Method_Score <- to_number(dat$Method_Score)
  dat$Method_P_Value <- to_number(dat$Method_P_Value)

  dat[order(
    dat$Method_Rank,
    dat$Method_P_Value,
    -abs(dat$Method_Score),
    dat$TF,
    na.last = TRUE
  ), , drop = FALSE]
}

deduplicate_tf_rows <- function(dat) {
  dat <- dat[!is.na(dat$TF) & dat$TF != "", , drop = FALSE]
  dat <- sort_by_rank_columns(dat)
  dat[!duplicated(dat$TF), , drop = FALSE]
}


# 3. ChEA3证据整理 ------------------------------------------------------------

build_chea3_evidence <- function(input_dir, prefix) {
  chea3_files <- get_method_csv_files(input_dir, "chea3")
  if (length(chea3_files) == 0L) {
    return(data.frame(TF = character(0)))
  }

  evidence_records <- list()
  integrated_top <- data.frame()

  for (file_name in chea3_files) {
    library_name <- get_library_name_from_file(file_name, "chea3", prefix)
    dat <- read_tf_csv(file_name)
    if (!"TF" %in% colnames(dat)) {
      next
    }

    dat$TF <- standardize_tf(dat$TF)
    dat <- dat[!is.na(dat$TF), , drop = FALSE]
    if (nrow(dat) == 0L) {
      next
    }

    rank_column <- first_existing_column(dat, c("Rank", "rank"))
    ranks <- if (!is.na(rank_column)) to_number(dat[[rank_column]]) else seq_len(nrow(dat))

    if (library_name == "Integrated--topRank") {
      score_column <- first_existing_column(dat, c("Score", "score"))
      library_column <- first_existing_column(dat, c("Library", "library"))
      integrated_top <- data.frame(
        TF = dat$TF,
        CheA3_Integrated_TopRank = ranks,
        CheA3_Integrated_Score = if (!is.na(score_column)) dat[[score_column]] else NA,
        CheA3_Integrated_Library = if (!is.na(library_column)) dat[[library_column]] else "",
        stringsAsFactors = FALSE
      )
      next
    }

    if (!library_name %in% CHEA3_EVIDENCE_LIBRARIES) {
      next
    }

    evidence_records[[length(evidence_records) + 1L]] <- data.frame(
      TF = dat$TF,
      Library_Name = library_name,
      Library_Rank = ranks,
      stringsAsFactors = FALSE
    )
  }

  evidence_table <- if (length(evidence_records) > 0L) {
    do.call(rbind, evidence_records)
  } else {
    data.frame(TF = character(0), Library_Name = character(0), Library_Rank = numeric(0))
  }

  evidence_summary <- if (nrow(evidence_table) > 0L) {
    split_records <- split(evidence_table, evidence_table$TF)
    do.call(rbind, lapply(names(split_records), function(tf) {
      x <- split_records[[tf]]
      x <- x[order(x$Library_Rank, x$Library_Name), , drop = FALSE]
      data.frame(
        TF = tf,
        CheA3_Library_Count = length(unique(x$Library_Name)),
        CheA3_Libraries = paste(
          paste0(x$Library_Name, "(rank=", x$Library_Rank, ")"),
          collapse = "; "
        ),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      TF = character(0),
      CheA3_Library_Count = integer(0),
      CheA3_Libraries = character(0)
    )
  }

  if (nrow(integrated_top) > 0L) {
    evidence_summary <- merge(
      evidence_summary,
      integrated_top,
      by = "TF",
      all = TRUE,
      sort = FALSE
    )
  }

  evidence_summary <- add_missing_columns(
    evidence_summary,
    c(
      "CheA3_Library_Count", "CheA3_Libraries",
      "CheA3_Integrated_TopRank", "CheA3_Integrated_Score",
      "CheA3_Integrated_Library"
    )
  )
  evidence_summary$CheA3_Library_Count[is.na(evidence_summary$CheA3_Library_Count)] <- 0
  evidence_summary$CheA3_Libraries[is.na(evidence_summary$CheA3_Libraries)] <- ""
  evidence_summary
}

add_chea3_evidence_columns <- function(dat, evidence) {
  dat$TF <- standardize_tf(dat$TF)
  evidence$TF <- standardize_tf(evidence$TF)

  evidence_columns <- c(
    "CheA3_Library_Count",
    "CheA3_Libraries",
    "CheA3_Integrated_TopRank",
    "CheA3_Integrated_Score",
    "CheA3_Integrated_Library"
  )
  evidence <- add_missing_columns(evidence, evidence_columns)

  match_index <- match(dat$TF, evidence$TF)
  for (column in evidence_columns) {
    dat[[column]] <- evidence[[column]][match_index]
  }

  dat$CheA3_Library_Count[is.na(dat$CheA3_Library_Count)] <- 0
  dat$CheA3_Libraries[is.na(dat$CheA3_Libraries)] <- ""
  dat
}


# 4. 单方法最终结果整理 --------------------------------------------------------

make_chea3_final <- function(input_dir, prefix, evidence) {
  target_file <- resolve_tf_result_file(file.path(
    input_dir,
    "chea3",
    paste0("chea3_", prefix, "_Integrated--topRank.csv")
  ))
  if (!file.exists(target_file)) {
    return(data.frame())
  }

  dat <- read_tf_csv(target_file)
  dat$TF <- standardize_tf(dat$TF)
  dat$Method <- "chea3"
  dat$Method_Rank <- to_number(dat$Rank)
  dat$Method_Score <- to_number(dat$Score)
  dat$Method_P_Value <- NA_real_
  dat$Method_Adjusted_P_Value <- NA_real_
  dat$Method_Direction <- ""

  dat <- add_chea3_evidence_columns(dat, evidence)
  deduplicate_tf_rows(dat)
}

make_enrichr_final <- function(input_dir, prefix, evidence) {
  enrichr_files <- get_method_csv_files(input_dir, "enrichr")
  enrichr_files <- enrichr_files[
    !grepl("Integrated--topRank[.]csv$", enrichr_files)
  ]
  if (length(enrichr_files) == 0L) {
    return(data.frame())
  }

  records <- list()
  for (file_name in enrichr_files) {
    library_name <- get_library_name_from_file(file_name, "enrichr", prefix)
    dat <- read_tf_csv(file_name)

    if (!"TF" %in% colnames(dat) && "Term" %in% colnames(dat)) {
      dat$TF <- sub("\\s.*$", "", trimws(as.character(dat$Term)))
    }
    if (!"TF" %in% colnames(dat)) {
      next
    }

    dat$TF <- standardize_tf(dat$TF)
    dat <- dat[!is.na(dat$TF), , drop = FALSE]
    if (nrow(dat) == 0L) {
      next
    }

    dat$P.value <- if ("P.value" %in% colnames(dat)) to_number(dat$P.value) else NA_real_
    dat$Adjusted.P.value <- if ("Adjusted.P.value" %in% colnames(dat)) {
      to_number(dat$Adjusted.P.value)
    } else {
      NA_real_
    }
    dat$Combined.Score <- if ("Combined.Score" %in% colnames(dat)) {
      to_number(dat$Combined.Score)
    } else {
      NA_real_
    }

    dat <- dat[order(
      dat$Adjusted.P.value,
      dat$P.value,
      -dat$Combined.Score,
      dat$TF,
      na.last = TRUE
    ), , drop = FALSE]
    dat$Library_Rank <- seq_len(nrow(dat))
    dat <- dat[!duplicated(dat$TF), , drop = FALSE]

    records[[length(records) + 1L]] <- data.frame(
      TF = dat$TF,
      Library_Name = library_name,
      Library_Rank = dat$Library_Rank,
      Library_Term = if ("Term" %in% colnames(dat)) dat$Term else "",
      Library_P_Value = dat$P.value,
      Library_Adjusted_P_Value = dat$Adjusted.P.value,
      Library_Combined_Score = dat$Combined.Score,
      stringsAsFactors = FALSE
    )
  }

  if (length(records) == 0L) {
    return(data.frame())
  }

  combined <- do.call(rbind, records)
  split_records <- split(combined, combined$TF)
  integrated <- do.call(rbind, lapply(names(split_records), function(tf) {
    x <- split_records[[tf]]
    x <- x[order(
      x$Library_Rank,
      x$Library_Adjusted_P_Value,
      x$Library_P_Value,
      -x$Library_Combined_Score,
      na.last = TRUE
    ), , drop = FALSE]

    data.frame(
      TF = tf,
      Best_Rank = min(x$Library_Rank, na.rm = TRUE),
      Mean_Rank = mean(x$Library_Rank, na.rm = TRUE),
      Library_Count = length(unique(x$Library_Name)),
      Best_P_Value = min(x$Library_P_Value, na.rm = TRUE),
      Best_Adjusted_P_Value = min(x$Library_Adjusted_P_Value, na.rm = TRUE),
      Best_Combined_Score = max(x$Library_Combined_Score, na.rm = TRUE),
      Library = paste(
        paste0(
          x$Library_Name,
          "(rank=", x$Library_Rank,
          ", adjP=", signif(x$Library_Adjusted_P_Value, 3),
          ")"
        ),
        collapse = "; "
      ),
      Top_Terms = paste(unique(x$Library_Term[seq_len(min(5, nrow(x)))]), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))

  integrated <- integrated[order(
    integrated$Best_Rank,
    integrated$Mean_Rank,
    -integrated$Library_Count,
    integrated$Best_Adjusted_P_Value,
    integrated$Best_P_Value,
    -integrated$Best_Combined_Score,
    integrated$TF,
    na.last = TRUE
  ), , drop = FALSE]

  integrated$Rank <- seq_len(nrow(integrated))
  integrated$Score <- integrated$Best_Combined_Score
  integrated$Method <- "enrichr"
  integrated$Method_Rank <- integrated$Rank
  integrated$Method_Score <- integrated$Score
  integrated$Method_P_Value <- integrated$Best_P_Value
  integrated$Method_Adjusted_P_Value <- integrated$Best_Adjusted_P_Value
  integrated$Method_Direction <- ""

  integrated <- add_chea3_evidence_columns(integrated, evidence)
  deduplicate_tf_rows(integrated)
}

make_ora_final <- function(input_dir, method_name, prefix, evidence) {
  target_file <- resolve_tf_result_file(file.path(
    input_dir,
    method_name,
    paste0(method_name, "_", prefix, "_tf_enrichment.csv")
  ))
  if (!file.exists(target_file)) {
    return(data.frame())
  }

  dat <- read_tf_csv(target_file)
  if (!"source" %in% colnames(dat)) {
    return(data.frame())
  }

  dat$TF <- standardize_tf(dat$source)
  dat$score <- if ("score" %in% colnames(dat)) to_number(dat$score) else NA_real_
  dat$p_value <- if ("p_value" %in% colnames(dat)) to_number(dat$p_value) else NA_real_
  dat <- dat[order(dat$p_value, -dat$score, dat$TF, na.last = TRUE), , drop = FALSE]
  dat$Rank <- seq_len(nrow(dat))
  dat$Method <- method_name
  dat$Method_Rank <- dat$Rank
  dat$Method_Score <- dat$score
  dat$Method_P_Value <- dat$p_value
  dat$Method_Adjusted_P_Value <- NA_real_
  dat$Method_Direction <- ""

  dat <- add_chea3_evidence_columns(dat, evidence)
  deduplicate_tf_rows(dat)
}

make_viper_final <- function(input_dir, prefix, evidence) {
  target_file <- resolve_tf_result_file(file.path(
    input_dir,
    "viper",
    paste0("viper_", prefix, "_tf_activity.csv")
  ))
  if (!file.exists(target_file)) {
    return(data.frame())
  }

  dat <- read_tf_csv(target_file)
  if (!"TF" %in% colnames(dat)) {
    return(data.frame())
  }

  dat$TF <- standardize_tf(dat$TF)
  score_columns <- setdiff(colnames(dat), "TF")
  score_matrix <- as.data.frame(lapply(dat[, score_columns, drop = FALSE], to_number))

  dat$Activity_Score_Mean <- rowMeans(score_matrix, na.rm = TRUE)
  dat$Score <- rowMeans(abs(as.matrix(score_matrix)), na.rm = TRUE)
  dat$Method_Direction <- ifelse(dat$Activity_Score_Mean > 0, "Activated", "Repressed")
  dat$Method_Direction[is.na(dat$Activity_Score_Mean)] <- ""
  dat <- dat[order(-dat$Score, dat$TF, na.last = TRUE), , drop = FALSE]
  dat$Rank <- seq_len(nrow(dat))
  dat$Method <- "viper"
  dat$Method_Rank <- dat$Rank
  dat$Method_Score <- dat$Score
  dat$Method_P_Value <- NA_real_
  dat$Method_Adjusted_P_Value <- NA_real_

  dat <- add_chea3_evidence_columns(dat, evidence)
  deduplicate_tf_rows(dat)
}

make_collectri_final <- function(input_dir, prefix, evidence) {
  target_file <- resolve_tf_result_file(file.path(
    input_dir,
    "collectri",
    paste0("collectri_", prefix, "_tf_activity.csv")
  ))
  if (!file.exists(target_file)) {
    return(data.frame())
  }

  dat <- read_tf_csv(target_file)
  if (!"source" %in% colnames(dat)) {
    return(data.frame())
  }

  dat$TF <- standardize_tf(dat$source)
  dat$score <- if ("score" %in% colnames(dat)) to_number(dat$score) else NA_real_
  dat$p_value <- if ("p_value" %in% colnames(dat)) to_number(dat$p_value) else NA_real_
  dat <- dat[!is.na(dat$TF), , drop = FALSE]
  if (nrow(dat) == 0L) {
    return(data.frame())
  }

  split_records <- split(dat, dat$TF)
  aggregated <- do.call(rbind, lapply(names(split_records), function(tf) {
    x <- split_records[[tf]]
    data.frame(
      TF = tf,
      Score = mean(abs(x$score), na.rm = TRUE),
      Activity_Score_Mean = mean(x$score, na.rm = TRUE),
      Method_P_Value = min(x$p_value, na.rm = TRUE),
      Condition_Scores = paste(
        paste0(x$condition, "=", signif(x$score, 4)),
        collapse = "; "
      ),
      stringsAsFactors = FALSE
    )
  }))

  aggregated$Method_P_Value[is.infinite(aggregated$Method_P_Value)] <- NA_real_
  aggregated$Method_Direction <- ifelse(
    aggregated$Activity_Score_Mean > 0,
    "Activated",
    "Repressed"
  )
  aggregated$Method_Direction[is.na(aggregated$Activity_Score_Mean)] <- ""
  aggregated <- aggregated[order(
    aggregated$Method_P_Value,
    -aggregated$Score,
    aggregated$TF,
    na.last = TRUE
  ), , drop = FALSE]
  aggregated$Rank <- seq_len(nrow(aggregated))
  aggregated$Method <- "collectri"
  aggregated$Method_Rank <- aggregated$Rank
  aggregated$Method_Score <- aggregated$Score
  aggregated$Method_Adjusted_P_Value <- NA_real_

  aggregated <- add_chea3_evidence_columns(aggregated, evidence)
  deduplicate_tf_rows(aggregated)
}

make_method_final_table <- function(input_dir, method_name, prefix, evidence) {
  if (method_name == "chea3") {
    return(make_chea3_final(input_dir, prefix, evidence))
  }
  if (method_name == "enrichr") {
    return(make_enrichr_final(input_dir, prefix, evidence))
  }
  if (method_name %in% c("dorothea", "trrust")) {
    return(make_ora_final(input_dir, method_name, prefix, evidence))
  }
  if (method_name == "viper") {
    return(make_viper_final(input_dir, prefix, evidence))
  }
  if (method_name == "collectri") {
    return(make_collectri_final(input_dir, prefix, evidence))
  }

  data.frame()
}

get_method_final_stem <- function(method_name, prefix) {
  if (method_name == "chea3") {
    return(paste0("chea3_", prefix, "_Integrated--topRank"))
  }
  if (method_name == "enrichr") {
    return(paste0("enrichr_", prefix, "_Integrated--topRank"))
  }

  paste0(method_name, "_", prefix, "_final_result")
}


# 5. 交集结果整理 --------------------------------------------------------------

get_tf_set <- function(method_table) {
  if (is.null(method_table) || nrow(method_table) == 0L || !"TF" %in% colnames(method_table)) {
    return(character(0))
  }

  method_table <- sort_by_rank_columns(method_table)
  unique(method_table$TF)
}

extract_method_values <- function(method_table, tf) {
  empty <- list(
    Rank = NA,
    Score = NA,
    P_Value = NA,
    Adjusted_P_Value = NA,
    Direction = ""
  )

  if (is.null(method_table) || nrow(method_table) == 0L) {
    return(empty)
  }

  idx <- match(tf, method_table$TF)
  if (is.na(idx)) {
    return(empty)
  }

  row <- method_table[idx, , drop = FALSE]
  list(
    Rank = row$Method_Rank,
    Score = row$Method_Score,
    P_Value = row$Method_P_Value,
    Adjusted_P_Value = row$Method_Adjusted_P_Value,
    Direction = row$Method_Direction
  )
}

format_number <- function(x, digits = 4) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || is.na(x)) {
    return("")
  }

  signif(x, digits)
}

format_method_evidence <- function(values) {
  parts <- character(0)

  if (!is.na(values$Rank)) {
    parts <- c(parts, paste0("rank=", format_number(values$Rank, 5)))
  }
  if (!is.na(values$Score)) {
    parts <- c(parts, paste0("score=", format_number(values$Score, 5)))
  }
  if (!is.na(values$P_Value)) {
    parts <- c(parts, paste0("p=", format_number(values$P_Value, 5)))
  }
  if (!is.na(values$Adjusted_P_Value)) {
    parts <- c(parts, paste0("adjP=", format_number(values$Adjusted_P_Value, 5)))
  }
  if (!is.na(values$Direction) && values$Direction != "") {
    parts <- c(parts, paste0("direction=", values$Direction))
  }

  paste(parts, collapse = "; ")
}

make_empty_intersection_table <- function(selected_methods = TF_METHODS) {
  base_columns <- data.frame(
    Consensus_Rank = integer(0),
    TF = character(0),
    Required_Methods = character(0),
    Source_Method_Count = integer(0),
    Source_Methods = character(0),
    Mean_Selected_Rank = numeric(0),
    Best_Selected_Rank = numeric(0),
    CheA3_Library_Count = integer(0),
    CheA3_Libraries = character(0),
    CheA3_Integrated_TopRank = numeric(0),
    stringsAsFactors = FALSE
  )

  for (method_name in selected_methods) {
    label <- TF_METHOD_LABELS[[method_name]]
    base_columns[[label]] <- character(0)
  }

  base_columns
}

make_intersection_table <- function(
    input_info,
    method_tables,
    method_sets,
    evidence,
    intersection_name,
    selected_methods) {
  selected_sets <- method_sets[selected_methods]
  selected_sets <- selected_sets[!vapply(selected_sets, is.null, logical(1))]

  if (length(selected_sets) == 0L || any(vapply(selected_sets, length, integer(1)) == 0L)) {
    return(make_empty_intersection_table(selected_methods))
  }

  intersected_tfs <- Reduce(intersect, selected_sets)
  if (length(intersected_tfs) == 0L) {
    return(make_empty_intersection_table(selected_methods))
  }

  all_method_sets <- method_sets[TF_METHODS]
  records <- lapply(sort(intersected_tfs), function(tf) {
    all_source_methods <- TF_METHODS[
      vapply(all_method_sets, function(x) tf %in% x, logical(1))
    ]
    selected_ranks <- vapply(selected_methods, function(method_name) {
      values <- extract_method_values(method_tables[[method_name]], tf)
      as.numeric(values$Rank)
    }, numeric(1))

    evidence_idx <- match(tf, evidence$TF)
    row <- data.frame(
      Consensus_Rank = NA_integer_,
      TF = tf,
      Required_Methods = paste(TF_METHOD_LABELS[selected_methods], collapse = ";"),
      Source_Method_Count = length(all_source_methods),
      Source_Methods = paste(TF_METHOD_LABELS[all_source_methods], collapse = ";"),
      Mean_Selected_Rank = mean(selected_ranks, na.rm = TRUE),
      Best_Selected_Rank = min(selected_ranks, na.rm = TRUE),
      CheA3_Library_Count = if (!is.na(evidence_idx)) evidence$CheA3_Library_Count[evidence_idx] else 0,
      CheA3_Libraries = if (!is.na(evidence_idx)) evidence$CheA3_Libraries[evidence_idx] else "",
      CheA3_Integrated_TopRank = if (!is.na(evidence_idx)) evidence$CheA3_Integrated_TopRank[evidence_idx] else NA,
      stringsAsFactors = FALSE
    )

    for (method_name in selected_methods) {
      label <- TF_METHOD_LABELS[[method_name]]
      values <- extract_method_values(method_tables[[method_name]], tf)
      row[[label]] <- format_method_evidence(values)
    }

    row
  })

  output <- do.call(rbind, records)
  output <- output[order(
    output$Mean_Selected_Rank,
    output$Best_Selected_Rank,
    -output$Source_Method_Count,
    -output$CheA3_Library_Count,
    output$TF,
    na.last = TRUE
  ), , drop = FALSE]
  output$Consensus_Rank <- seq_len(nrow(output))
  output[seq_len(min(INTERSECTION_TOP_N_TO_REPORT, nrow(output))), , drop = FALSE]
}

make_intersection_summary <- function(
    input_info,
    intersection_name,
    selected_methods,
    method_sets,
    result_table) {
  selected_sets <- method_sets[selected_methods]
  selected_counts <- vapply(selected_sets, length, integer(1))
  intersected_count <- if (length(selected_sets) == 0L || any(selected_counts == 0L)) {
    0L
  } else {
    length(Reduce(intersect, selected_sets))
  }

  data.frame(
    Input_Type = input_info$Input_Type,
    TF_Analysis_Name = input_info$Analysis_Name,
    Intersection_Name = intersection_name,
    Required_Methods = paste(TF_METHOD_LABELS[selected_methods], collapse = ";"),
    Number_Of_Methods = length(selected_methods),
    Intersected_TF_Count = intersected_count,
    Reported_Top_N = nrow(result_table),
    stringsAsFactors = FALSE
  )
}


# 6. 单个TF输入方案整合 --------------------------------------------------------

process_one_tf_input <- function(input_info) {
  prefix <- input_info$Method_Prefix
  output_base_dir <- file.path(
    TF_SUMMARY_ROOT,
    tolower(input_info$Input_Type),
    sanitize_file_name(input_info$Analysis_Name)
  )
  dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

  evidence <- build_chea3_evidence(input_info$Input_Dir, prefix)
  method_tables <- list()
  method_summary_records <- list()

  for (method_name in TF_METHODS) {
    method_table <- make_method_final_table(
      input_dir = input_info$Input_Dir,
      method_name = method_name,
      prefix = prefix,
      evidence = evidence
    )

    method_tables[[method_name]] <- method_table
    method_output_dir <- file.path(output_base_dir, "method_final", method_name)
    file_stem <- get_method_final_stem(method_name, prefix)
    method_csv_file <- write_csv_with_report_previews(
      dat = method_table,
      csv_file = file.path(method_output_dir, paste0(file_stem, ".csv")),
      n_rows = TABLE_PREVIEW_ROWS
    )

    method_summary_records[[method_name]] <- data.frame(
      Input_Type = input_info$Input_Type,
      TF_Analysis_Name = input_info$Analysis_Name,
      Method = method_name,
      Method_Label = TF_METHOD_LABELS[[method_name]],
      Final_TF_Count = nrow(method_table),
      Final_Result_File = method_csv_file,
      stringsAsFactors = FALSE
    )
  }

  method_summary <- do.call(rbind, method_summary_records)
  write_csv_with_report_previews(
    dat = method_summary,
    csv_file = file.path(output_base_dir, "method_final", "summary", "method_final_summary.csv"),
    n_rows = TABLE_PREVIEW_ROWS
  )

  method_sets <- lapply(method_tables, get_tf_set)
  method_sets <- method_sets[TF_METHODS]

  intersection_records <- list()
  for (intersection_name in INTERSECTION_SCHEMES_TO_RUN) {
    selected_methods <- TF_INTERSECTION_SCHEMES[[intersection_name]]
    selected_methods <- intersect(selected_methods, TF_METHODS)
    if (length(selected_methods) < 2L) {
      next
    }

    intersection_output_dir <- file.path(
      output_base_dir,
      "intersections",
      sanitize_file_name(intersection_name)
    )
    dir.create(intersection_output_dir, recursive = TRUE, showWarnings = FALSE)

    intersection_table <- make_intersection_table(
      input_info = input_info,
      method_tables = method_tables,
      method_sets = method_sets,
      evidence = evidence,
      intersection_name = intersection_name,
      selected_methods = selected_methods
    )
    intersection_csv_file <- write_csv_with_report_previews(
      dat = intersection_table,
      csv_file = file.path(intersection_output_dir, "candidates", "top10_tf_candidates.csv"),
      n_rows = TABLE_PREVIEW_ROWS
    )

    intersection_summary <- make_intersection_summary(
      input_info = input_info,
      intersection_name = intersection_name,
      selected_methods = selected_methods,
      method_sets = method_sets,
      result_table = intersection_table
    )
    intersection_summary$Result_File <- intersection_csv_file
    intersection_summary_file <- write_csv_with_report_previews(
      dat = intersection_summary,
      csv_file = file.path(intersection_output_dir, "summary", "summary.csv"),
      n_rows = TABLE_PREVIEW_ROWS
    )
    intersection_summary$Summary_File <- intersection_summary_file

    intersection_records[[intersection_name]] <- intersection_summary
  }

  intersection_summary_table <- if (length(intersection_records) > 0L) {
    do.call(rbind, intersection_records)
  } else {
    data.frame()
  }

  if (nrow(intersection_summary_table) > 0L) {
    write_csv_with_report_previews(
      dat = intersection_summary_table,
      csv_file = file.path(output_base_dir, "intersections", "summary", "intersection_summary.csv"),
      n_rows = TABLE_PREVIEW_ROWS
    )
  }

  data.frame(
    Input_Type = input_info$Input_Type,
    TF_Analysis_Name = input_info$Analysis_Name,
    Methods_Integrated = nrow(method_summary),
    Intersection_Schemes = nrow(intersection_summary_table),
    Output_Directory = output_base_dir,
    stringsAsFactors = FALSE
  )
}


# 7. 主流程 -------------------------------------------------------------------

clean_output_dirs()

tf_inputs <- get_tf_input_dirs(TF_RESULT_ROOT, INPUT_TYPES_TO_RUN)
if (length(tf_inputs) == 0L) {
  stop("No TF result directories were found under: ", TF_RESULT_ROOT)
}

cat("\nRunning TF enrichment result integration...\n")
cat("TF input schemes: ", length(tf_inputs), "\n", sep = "")
cat("Intersection schemes per input: ", length(INTERSECTION_SCHEMES_TO_RUN), "\n", sep = "")
cat("Intersection candidate sets: all TFs from each method_final result\n")

parallel_strategy <- setup_parallel_strategy(
  total_tasks = length(tf_inputs),
  inner_label = "TF integration inner workers",
  nested_label = "Nested workers"
)

summary_records <- run_indexed_tasks_with_progress(
  total_tasks = length(tf_inputs),
  workers = parallel_strategy$task_workers,
  progress_label = "TF integration",
  task_function = function(i) {
    process_one_tf_input(tf_inputs[[i]])
  }
)
stop_on_parallel_errors(summary_records, task_ids = seq_along(tf_inputs), label = "TF integration tasks")

summary_table <- do.call(rbind, summary_records)
rownames(summary_table) <- NULL
summary_csv_file <- write_csv_with_report_previews(
  dat = summary_table,
  csv_file = file.path(TF_SUMMARY_ROOT, "run_summary", "summary.csv"),
  n_rows = TABLE_PREVIEW_ROWS
)

cat("\nTF integration finished.\n")
cat("Summary table: ", summary_csv_file, "\n", sep = "")
cat("Table root:    ", TF_SUMMARY_ROOT, "\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
