# GSE114012 ATF3 TF rank and intersection search
#
# 目的：
# 1. 汇总每种DEG方案下，单个TF方法中ATF3在不同合理排序口径下的排名；
# 2. 枚举TF方法组合，寻找能让ATF3在最多DEG方案中成为交集候选top1的通用组合。
#
# 注意：这里仅评价TF富集/活性推断结果。DEG只作为输入方案，不对DEG基因表做排名。


# 0. Configuration ------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"
TARGET_TF <- "ATF3"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TF_SUMMARY_ROOT <- file.path(RESULT_ROOT, "TF_summary")
OUTPUT_ROOT <- file.path(RESULT_ROOT, "ATF3_evidence_summary", "tf_rank_intersection_search")

TF_METHODS <- c("dorothea", "chea3", "viper", "enrichr", "trrust", "collectri")
TF_METHOD_LABELS <- c(
  dorothea = "DoRothEA",
  chea3 = "ChEA3",
  viper = "VIPER",
  enrichr = "ENRICHR",
  trrust = "TRRUST",
  collectri = "CollecTRI"
)

# 单方法排序口径。
# standard_method_rank：09整合脚本已保存的Method_Rank/Rank，是默认主口径。
# method_p_value_asc / method_adj_p_value_asc：仅对有P值/FDR的方法有效。
# method_score_primary：按各方法分数方向排序；ChEA3分数越小越好，其他方法分数越大越好。
# chea3_external_rank_asc：所有method_final表中附带的ChEA3外部证据rank，仅作为外部证据辅助口径。
SINGLE_METHOD_SORT_RULES <- c(
  "standard_method_rank",
  "method_p_value_asc",
  "method_adj_p_value_asc",
  "method_score_primary",
  "chea3_external_rank_asc",
  "chea3_library_count_desc_then_rank"
)

# 交集共识排序口径。
# official_mean_rank与09号脚本一致。
INTERSECTION_SORT_RULES <- c(
  "official_mean_rank",
  "best_rank_first",
  "worst_rank_first",
  "chea3_integrated_rank_first",
  "chea3_library_then_integrated_rank",
  "source_count_then_mean_rank",
  "best_method_p_value_first",
  "fisher_p_value_first",
  "mean_neglog10_p_first"
)

TOP_N_TO_SAVE <- 10

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
options(width = 200)


# 1. Helpers ------------------------------------------------------------------

clean_path <- function(x) {
  gsub("\\\\", "/", x)
}

to_number <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

standardize_tf <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "" | is.na(x) | x == "NA"] <- NA_character_
  x
}

first_existing <- function(dat, candidates) {
  hit <- candidates[candidates %in% colnames(dat)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

cell_num <- function(row, candidates) {
  col <- first_existing(row, candidates)
  if (is.na(col)) {
    return(NA_real_)
  }
  to_number(row[[col]][1])
}

cell_chr <- function(row, candidates) {
  col <- first_existing(row, candidates)
  if (is.na(col)) {
    return(NA_character_)
  }
  as.character(row[[col]][1])
}

safe_read_csv <- function(file) {
  tryCatch(
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

collapse_unique <- function(x, sep = ";") {
  x <- unique(as.character(x[!is.na(x) & trimws(as.character(x)) != ""]))
  if (length(x) == 0) {
    return("")
  }
  paste(x, collapse = sep)
}

extract_method_file_meta <- function(file) {
  parts <- strsplit(clean_path(file), "/", fixed = TRUE)[[1]]
  idx <- match("TF_summary", parts)
  if (is.na(idx) || length(parts) < idx + 6) {
    return(NULL)
  }

  if (!identical(parts[idx + 1], "deg")) {
    return(NULL)
  }
  if (!identical(parts[idx + 3], "method_final")) {
    return(NULL)
  }

  list(
    Analysis_Name = parts[idx + 2],
    Method = parts[idx + 4]
  )
}

normalize_method_table <- function(file) {
  meta <- extract_method_file_meta(file)
  dat <- safe_read_csv(file)

  if (is.null(meta) || is.null(dat) || nrow(dat) == 0) {
    return(NULL)
  }

  tf_col <- first_existing(dat, c("TF", "source", "Source", "Transcription_Factor", "tf", "Term", "term"))
  if (is.na(tf_col)) {
    return(NULL)
  }

  out <- data.frame(
    Analysis_Name = meta$Analysis_Name,
    Method = meta$Method,
    TF = standardize_tf(dat[[tf_col]]),
    Method_Rank = to_number(dat[[first_existing(dat, c("Method_Rank", "Rank", "Best_Rank"))]]),
    Method_Score = to_number(dat[[first_existing(dat, c("Method_Score", "Score", "score", "Best_Combined_Score", "Activity_Score_Mean", "statistic"))]]),
    Method_P_Value = to_number(dat[[first_existing(dat, c("Method_P_Value", "p_value", "P.Value", "p.value", "Best_P_Value"))]]),
    Method_Adjusted_P_Value = to_number(dat[[first_existing(dat, c("Method_Adjusted_P_Value", "Best_Adjusted_P_Value", "adj.P.Val", "p.adjust", "FDR"))]]),
    Method_Direction = as.character(dat[[first_existing(dat, c("Method_Direction", "Direction", "direction"))]]),
    CheA3_Library_Count = to_number(dat[[first_existing(dat, c("CheA3_Library_Count"))]]),
    CheA3_Libraries = as.character(dat[[first_existing(dat, c("CheA3_Libraries"))]]),
    CheA3_Integrated_TopRank = to_number(dat[[first_existing(dat, c("CheA3_Integrated_TopRank"))]]),
    Result_File = file,
    stringsAsFactors = FALSE
  )

  out$Method_Rank[is.na(out$Method_Rank)] <- seq_len(nrow(out))[is.na(out$Method_Rank)]
  out$Method_Direction[is.na(out$Method_Direction)] <- ""
  out$CheA3_Library_Count[is.na(out$CheA3_Library_Count)] <- 0
  out <- out[!is.na(out$TF) & out$TF != "", , drop = FALSE]

  out <- out[order(
    out$Method_Rank,
    out$Method_P_Value,
    -abs(out$Method_Score),
    out$TF,
    na.last = TRUE
  ), , drop = FALSE]
  out[!duplicated(out$TF), , drop = FALSE]
}

rank_one_table <- function(dat, method, sort_rule) {
  if (nrow(dat) == 0) {
    return(NULL)
  }

  dat <- dat
  dat$Sort_Rule <- sort_rule

  if (sort_rule == "standard_method_rank") {
    ord <- order(dat$Method_Rank, dat$Method_P_Value, -abs(dat$Method_Score), dat$TF, na.last = TRUE)
  } else if (sort_rule == "method_p_value_asc") {
    if (all(is.na(dat$Method_P_Value))) {
      return(NULL)
    }
    ord <- order(dat$Method_P_Value, dat$Method_Rank, -abs(dat$Method_Score), dat$TF, na.last = TRUE)
  } else if (sort_rule == "method_adj_p_value_asc") {
    if (all(is.na(dat$Method_Adjusted_P_Value))) {
      return(NULL)
    }
    ord <- order(dat$Method_Adjusted_P_Value, dat$Method_P_Value, dat$Method_Rank, dat$TF, na.last = TRUE)
  } else if (sort_rule == "method_score_primary") {
    if (all(is.na(dat$Method_Score))) {
      return(NULL)
    }
    if (identical(method, "chea3")) {
      ord <- order(dat$Method_Score, dat$Method_Rank, dat$TF, na.last = TRUE)
    } else {
      ord <- order(-dat$Method_Score, dat$Method_Rank, dat$TF, na.last = TRUE)
    }
  } else if (sort_rule == "chea3_external_rank_asc") {
    if (all(is.na(dat$CheA3_Integrated_TopRank))) {
      return(NULL)
    }
    ord <- order(dat$CheA3_Integrated_TopRank, dat$Method_Rank, dat$TF, na.last = TRUE)
  } else if (sort_rule == "chea3_library_count_desc_then_rank") {
    if (all(dat$CheA3_Library_Count == 0) && all(is.na(dat$CheA3_Integrated_TopRank))) {
      return(NULL)
    }
    ord <- order(-dat$CheA3_Library_Count, dat$CheA3_Integrated_TopRank, dat$Method_Rank, dat$TF, na.last = TRUE)
  } else {
    stop("Unknown single method sort rule: ", sort_rule)
  }

  ranked <- dat[ord, , drop = FALSE]
  ranked$Rank_By_Sort_Rule <- seq_len(nrow(ranked))
  ranked
}

get_tf_row <- function(dat, tf) {
  row <- dat[dat$TF == standardize_tf(tf), , drop = FALSE]
  if (nrow(row) == 0) {
    return(NULL)
  }
  row[1, , drop = FALSE]
}

extract_method_values <- function(method_tables, method, tf) {
  dat <- method_tables[[method]]
  if (is.null(dat)) {
    return(NULL)
  }
  get_tf_row(dat, tf)
}

fisher_p <- function(p_values) {
  p_values <- p_values[is.finite(p_values) & p_values > 0 & p_values <= 1]
  if (length(p_values) == 0) {
    return(NA_real_)
  }
  stats::pchisq(-2 * sum(log(p_values)), df = 2 * length(p_values), lower.tail = FALSE)
}


# 2. Load DEG method_final tables --------------------------------------------

method_files <- list.files(
  TF_SUMMARY_ROOT,
  pattern = "[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
method_files <- method_files[
  grepl("/deg/", clean_path(method_files)) &
    grepl("/method_final/", clean_path(method_files)) &
    !grepl("/method_final/summary/", clean_path(method_files))
]

method_tables_long <- do.call(
  rbind,
  Filter(Negate(is.null), lapply(method_files, normalize_method_table))
)

stopifnot(nrow(method_tables_long) > 0)

analysis_names <- sort(unique(method_tables_long$Analysis_Name))

method_tables_by_analysis <- lapply(analysis_names, function(analysis_name) {
  by_method <- split(
    method_tables_long[method_tables_long$Analysis_Name == analysis_name, , drop = FALSE],
    method_tables_long$Method[method_tables_long$Analysis_Name == analysis_name]
  )
  by_method[TF_METHODS]
})
names(method_tables_by_analysis) <- analysis_names


# 3. Single method ATF3 rank under different sort rules ------------------------

single_rank_records <- list()

for (analysis_name in analysis_names) {
  for (method in TF_METHODS) {
    dat <- method_tables_by_analysis[[analysis_name]][[method]]
    if (is.null(dat) || nrow(dat) == 0) {
      next
    }

    for (sort_rule in SINGLE_METHOD_SORT_RULES) {
      ranked <- rank_one_table(dat, method, sort_rule)
      if (is.null(ranked)) {
        next
      }

      atf3 <- get_tf_row(ranked, TARGET_TF)
      if (is.null(atf3)) {
        next
      }

      single_rank_records[[length(single_rank_records) + 1L]] <- data.frame(
        Analysis_Name = analysis_name,
        Method = method,
        Method_Label = TF_METHOD_LABELS[[method]],
        Sort_Rule = sort_rule,
        ATF3_Rank = atf3$Rank_By_Sort_Rule,
        ATF3_Standard_Method_Rank = atf3$Method_Rank,
        ATF3_Method_Score = atf3$Method_Score,
        ATF3_Method_P_Value = atf3$Method_P_Value,
        ATF3_Method_Adjusted_P_Value = atf3$Method_Adjusted_P_Value,
        ATF3_Method_Direction = atf3$Method_Direction,
        ATF3_CheA3_Integrated_TopRank = atf3$CheA3_Integrated_TopRank,
        ATF3_CheA3_Library_Count = atf3$CheA3_Library_Count,
        ATF3_Top1 = atf3$Rank_By_Sort_Rule == 1,
        stringsAsFactors = FALSE
      )
    }
  }
}

single_method_rank <- do.call(rbind, single_rank_records)
single_method_rank <- single_method_rank[order(
  single_method_rank$Analysis_Name,
  single_method_rank$Method,
  single_method_rank$Sort_Rule
), , drop = FALSE]

single_method_rank_top1 <- single_method_rank[
  single_method_rank$ATF3_Top1,
  ,
  drop = FALSE
]

single_method_top1_summary <- aggregate(
  ATF3_Top1 ~ Method + Method_Label + Sort_Rule,
  data = single_method_rank,
  FUN = sum
)
colnames(single_method_top1_summary)[colnames(single_method_top1_summary) == "ATF3_Top1"] <- "ATF3_Top1_Analysis_Count"
single_method_top1_summary$Top1_Analyses <- vapply(seq_len(nrow(single_method_top1_summary)), function(i) {
  row <- single_method_top1_summary[i, , drop = FALSE]
  collapse_unique(single_method_rank$Analysis_Name[
    single_method_rank$Method == row$Method &
      single_method_rank$Sort_Rule == row$Sort_Rule &
      single_method_rank$ATF3_Top1
  ])
}, character(1))
single_method_top1_summary <- single_method_top1_summary[order(
  -single_method_top1_summary$ATF3_Top1_Analysis_Count,
  single_method_top1_summary$Method,
  single_method_top1_summary$Sort_Rule
), , drop = FALSE]


# 4. Intersection combination search -----------------------------------------

make_intersection_candidates <- function(method_tables, selected_methods) {
  selected_tables <- method_tables[selected_methods]
  if (any(vapply(selected_tables, is.null, logical(1)))) {
    return(NULL)
  }

  selected_sets <- lapply(selected_tables, function(dat) dat$TF)
  if (length(selected_sets) == 0 || any(vapply(selected_sets, length, integer(1)) == 0L)) {
    return(NULL)
  }

  intersected_tfs <- Reduce(intersect, selected_sets)
  if (length(intersected_tfs) == 0) {
    return(NULL)
  }

  all_sets <- lapply(method_tables[TF_METHODS], function(dat) {
    if (is.null(dat)) character(0) else dat$TF
  })

  records <- lapply(sort(intersected_tfs), function(tf) {
    selected_rows <- Filter(
      Negate(is.null),
      lapply(selected_methods, function(method) extract_method_values(method_tables, method, tf))
    )
    selected_ranks <- vapply(selected_rows, function(row) row$Method_Rank, numeric(1))
    selected_scores <- vapply(selected_rows, function(row) row$Method_Score, numeric(1))
    selected_pvalues <- vapply(selected_rows, function(row) row$Method_P_Value, numeric(1))
    selected_adjp <- vapply(selected_rows, function(row) row$Method_Adjusted_P_Value, numeric(1))

    all_source_methods <- TF_METHODS[
      vapply(all_sets, function(x) tf %in% x, logical(1))
    ]

    evidence_rows <- Filter(
      Negate(is.null),
      lapply(TF_METHODS, function(method) extract_method_values(method_tables, method, tf))
    )
    evidence_row <- if (length(evidence_rows) > 0) evidence_rows[[1]] else NULL

    best_p <- suppressWarnings(min(selected_pvalues, na.rm = TRUE))
    if (!is.finite(best_p)) best_p <- NA_real_
    best_adjp <- suppressWarnings(min(selected_adjp, na.rm = TRUE))
    if (!is.finite(best_adjp)) best_adjp <- NA_real_

    data.frame(
      TF = tf,
      Required_Methods = paste(TF_METHOD_LABELS[selected_methods], collapse = ";"),
      Required_Method_Keys = paste(selected_methods, collapse = "+"),
      Source_Method_Count = length(all_source_methods),
      Source_Methods = paste(TF_METHOD_LABELS[all_source_methods], collapse = ";"),
      Mean_Selected_Rank = mean(selected_ranks, na.rm = TRUE),
      Median_Selected_Rank = stats::median(selected_ranks, na.rm = TRUE),
      Best_Selected_Rank = min(selected_ranks, na.rm = TRUE),
      Worst_Selected_Rank = max(selected_ranks, na.rm = TRUE),
      Mean_Selected_Score = mean(selected_scores, na.rm = TRUE),
      Best_Method_P_Value = best_p,
      Best_Method_Adjusted_P_Value = best_adjp,
      Fisher_P_Value = fisher_p(selected_pvalues),
      Mean_NegLog10_P = ifelse(
        any(is.finite(selected_pvalues) & selected_pvalues > 0),
        mean(-log10(selected_pvalues[is.finite(selected_pvalues) & selected_pvalues > 0])),
        NA_real_
      ),
      CheA3_Library_Count = if (is.null(evidence_row)) 0 else evidence_row$CheA3_Library_Count,
      CheA3_Integrated_TopRank = if (is.null(evidence_row)) NA_real_ else evidence_row$CheA3_Integrated_TopRank,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, records)
}

rank_intersection_candidates <- function(candidates, sort_rule) {
  if (is.null(candidates) || nrow(candidates) == 0) {
    return(NULL)
  }

  dat <- candidates

  if (sort_rule == "official_mean_rank") {
    ord <- order(
      dat$Mean_Selected_Rank,
      dat$Best_Selected_Rank,
      -dat$Source_Method_Count,
      -dat$CheA3_Library_Count,
      dat$TF,
      na.last = TRUE
    )
  } else if (sort_rule == "best_rank_first") {
    ord <- order(
      dat$Best_Selected_Rank,
      dat$Mean_Selected_Rank,
      -dat$Source_Method_Count,
      -dat$CheA3_Library_Count,
      dat$TF,
      na.last = TRUE
    )
  } else if (sort_rule == "worst_rank_first") {
    ord <- order(
      dat$Worst_Selected_Rank,
      dat$Mean_Selected_Rank,
      dat$Best_Selected_Rank,
      -dat$Source_Method_Count,
      -dat$CheA3_Library_Count,
      dat$TF,
      na.last = TRUE
    )
  } else if (sort_rule == "chea3_integrated_rank_first") {
    ord <- order(
      dat$CheA3_Integrated_TopRank,
      dat$Mean_Selected_Rank,
      dat$Best_Selected_Rank,
      -dat$Source_Method_Count,
      -dat$CheA3_Library_Count,
      dat$TF,
      na.last = TRUE
    )
  } else if (sort_rule == "chea3_library_then_integrated_rank") {
    ord <- order(
      -dat$CheA3_Library_Count,
      dat$CheA3_Integrated_TopRank,
      dat$Mean_Selected_Rank,
      dat$Best_Selected_Rank,
      dat$TF,
      na.last = TRUE
    )
  } else if (sort_rule == "source_count_then_mean_rank") {
    ord <- order(
      -dat$Source_Method_Count,
      dat$Mean_Selected_Rank,
      dat$Best_Selected_Rank,
      -dat$CheA3_Library_Count,
      dat$TF,
      na.last = TRUE
    )
  } else if (sort_rule == "best_method_p_value_first") {
    ord <- order(
      dat$Best_Method_P_Value,
      dat$Mean_Selected_Rank,
      dat$Best_Selected_Rank,
      -dat$Source_Method_Count,
      dat$TF,
      na.last = TRUE
    )
  } else if (sort_rule == "fisher_p_value_first") {
    ord <- order(
      dat$Fisher_P_Value,
      dat$Mean_Selected_Rank,
      dat$Best_Selected_Rank,
      -dat$Source_Method_Count,
      dat$TF,
      na.last = TRUE
    )
  } else if (sort_rule == "mean_neglog10_p_first") {
    ord <- order(
      -dat$Mean_NegLog10_P,
      dat$Mean_Selected_Rank,
      dat$Best_Selected_Rank,
      -dat$Source_Method_Count,
      dat$TF,
      na.last = TRUE
    )
  } else {
    stop("Unknown intersection sort rule: ", sort_rule)
  }

  out <- dat[ord, , drop = FALSE]
  out$Consensus_Rank <- seq_len(nrow(out))
  out
}

combo_records <- list()
detail_records <- list()
top_candidate_records <- list()

all_combinations <- unlist(
  lapply(2:length(TF_METHODS), function(k) {
    utils::combn(TF_METHODS, k, simplify = FALSE)
  }),
  recursive = FALSE
)

for (selected_methods in all_combinations) {
  combination_key <- paste(selected_methods, collapse = "+")
  combination_label <- paste(TF_METHOD_LABELS[selected_methods], collapse = "+")

  for (analysis_name in analysis_names) {
    method_tables <- method_tables_by_analysis[[analysis_name]]
    candidates <- make_intersection_candidates(method_tables, selected_methods)
    if (is.null(candidates)) {
      next
    }

    for (sort_rule in INTERSECTION_SORT_RULES) {
      ranked <- rank_intersection_candidates(candidates, sort_rule)
      if (is.null(ranked) || nrow(ranked) == 0) {
        next
      }

      atf3 <- ranked[ranked$TF == standardize_tf(TARGET_TF), , drop = FALSE]
      if (nrow(atf3) == 0) {
        atf3_rank <- NA_real_
        atf3_mean_rank <- NA_real_
        atf3_best_rank <- NA_real_
        atf3_worst_rank <- NA_real_
        atf3_best_p <- NA_real_
        atf3_fisher_p <- NA_real_
        atf3_chea3_rank <- NA_real_
      } else {
        atf3 <- atf3[1, , drop = FALSE]
        atf3_rank <- atf3$Consensus_Rank
        atf3_mean_rank <- atf3$Mean_Selected_Rank
        atf3_best_rank <- atf3$Best_Selected_Rank
        atf3_worst_rank <- atf3$Worst_Selected_Rank
        atf3_best_p <- atf3$Best_Method_P_Value
        atf3_fisher_p <- atf3$Fisher_P_Value
        atf3_chea3_rank <- atf3$CheA3_Integrated_TopRank
      }

      detail_records[[length(detail_records) + 1L]] <- data.frame(
        Analysis_Name = analysis_name,
        Combination_Key = combination_key,
        Combination_Label = combination_label,
        Number_Of_Methods = length(selected_methods),
        Sort_Rule = sort_rule,
        Intersected_TF_Count = nrow(candidates),
        ATF3_In_Intersection = nrow(atf3) > 0,
        ATF3_Consensus_Rank = atf3_rank,
        ATF3_Mean_Selected_Rank = atf3_mean_rank,
        ATF3_Best_Selected_Rank = atf3_best_rank,
        ATF3_Worst_Selected_Rank = atf3_worst_rank,
        ATF3_Best_Method_P_Value = atf3_best_p,
        ATF3_Fisher_P_Value = atf3_fisher_p,
        ATF3_CheA3_Integrated_TopRank = atf3_chea3_rank,
        stringsAsFactors = FALSE
      )

      top_save <- ranked[seq_len(min(TOP_N_TO_SAVE, nrow(ranked))), , drop = FALSE]
      top_save$Analysis_Name <- analysis_name
      top_save$Combination_Key <- combination_key
      top_save$Combination_Label <- combination_label
      top_save$Sort_Rule <- sort_rule
      top_candidate_records[[length(top_candidate_records) + 1L]] <- top_save
    }
  }
}

intersection_detail <- do.call(rbind, detail_records)
intersection_top_candidates <- do.call(rbind, top_candidate_records)

summary_split <- split(
  intersection_detail,
  paste(intersection_detail$Combination_Key, intersection_detail$Sort_Rule, sep = "___")
)

intersection_summary <- do.call(
  rbind,
  lapply(summary_split, function(dat) {
    ranks <- dat$ATF3_Consensus_Rank
    data.frame(
      Combination_Key = dat$Combination_Key[1],
      Combination_Label = dat$Combination_Label[1],
      Number_Of_Methods = dat$Number_Of_Methods[1],
      Sort_Rule = dat$Sort_Rule[1],
      Analyses_Checked = length(unique(dat$Analysis_Name)),
      ATF3_In_Intersection_Count = sum(dat$ATF3_In_Intersection, na.rm = TRUE),
      ATF3_Top1_Count = sum(ranks == 1, na.rm = TRUE),
      ATF3_Top3_Count = sum(ranks <= 3, na.rm = TRUE),
      ATF3_Top10_Count = sum(ranks <= 10, na.rm = TRUE),
      ATF3_Top1_Analyses = collapse_unique(dat$Analysis_Name[ranks == 1]),
      ATF3_Top3_Analyses = collapse_unique(dat$Analysis_Name[ranks <= 3]),
      ATF3_Top10_Analyses = collapse_unique(dat$Analysis_Name[ranks <= 10]),
      Mean_ATF3_Consensus_Rank = ifelse(any(!is.na(ranks)), mean(ranks, na.rm = TRUE), NA_real_),
      Median_ATF3_Consensus_Rank = ifelse(any(!is.na(ranks)), stats::median(ranks, na.rm = TRUE), NA_real_),
      Max_ATF3_Consensus_Rank = ifelse(any(!is.na(ranks)), max(ranks, na.rm = TRUE), NA_real_),
      stringsAsFactors = FALSE
    )
  })
)

intersection_summary <- intersection_summary[order(
  -intersection_summary$ATF3_Top1_Count,
  -intersection_summary$ATF3_Top3_Count,
  -intersection_summary$ATF3_Top10_Count,
  intersection_summary$Number_Of_Methods,
  intersection_summary$Mean_ATF3_Consensus_Rank,
  intersection_summary$Combination_Key,
  intersection_summary$Sort_Rule
), , drop = FALSE]

best_any <- intersection_summary[1, , drop = FALSE]
best_official <- intersection_summary[
  intersection_summary$Sort_Rule == "official_mean_rank",
  ,
  drop = FALSE
]
best_official <- best_official[order(
  -best_official$ATF3_Top1_Count,
  -best_official$ATF3_Top3_Count,
  -best_official$ATF3_Top10_Count,
  best_official$Number_Of_Methods,
  best_official$Mean_ATF3_Consensus_Rank
), , drop = FALSE][1, , drop = FALSE]

best_any_top_candidates <- intersection_top_candidates[
  intersection_top_candidates$Combination_Key == best_any$Combination_Key &
    intersection_top_candidates$Sort_Rule == best_any$Sort_Rule,
  ,
  drop = FALSE
]

best_official_top_candidates <- intersection_top_candidates[
  intersection_top_candidates$Combination_Key == best_official$Combination_Key &
    intersection_top_candidates$Sort_Rule == best_official$Sort_Rule,
  ,
  drop = FALSE
]


# 5. Save outputs -------------------------------------------------------------

write.csv(
  single_method_rank,
  file.path(OUTPUT_ROOT, "atf3_single_method_rank_by_sort_rule.csv"),
  row.names = FALSE
)
write.csv(
  single_method_top1_summary,
  file.path(OUTPUT_ROOT, "atf3_single_method_top1_summary.csv"),
  row.names = FALSE
)
write.csv(
  intersection_detail,
  file.path(OUTPUT_ROOT, "intersection_combination_atf3_rank_detail.csv"),
  row.names = FALSE
)
write.csv(
  intersection_summary,
  file.path(OUTPUT_ROOT, "intersection_combination_search_summary.csv"),
  row.names = FALSE
)
write.csv(
  best_any_top_candidates,
  file.path(OUTPUT_ROOT, "best_any_sort_rule_top10_candidates.csv"),
  row.names = FALSE
)
write.csv(
  best_official_top_candidates,
  file.path(OUTPUT_ROOT, "best_official_mean_rank_top10_candidates.csv"),
  row.names = FALSE
)

single_standard <- single_method_rank[
  single_method_rank$Sort_Rule == "standard_method_rank",
  c(
    "Analysis_Name", "Method_Label", "ATF3_Rank", "ATF3_Method_Score",
    "ATF3_Method_P_Value", "ATF3_Method_Adjusted_P_Value",
    "ATF3_Method_Direction", "ATF3_Top1"
  ),
  drop = FALSE
]

standard_lines <- c(
  "# ATF3 TF rank and intersection search",
  "",
  "## Single TF method: standard ranks",
  "",
  paste(
    "| Analysis | Method | ATF3 rank | Score | P value | adj.P | Direction | Top1 |",
    "|---|---|---:|---:|---:|---:|---|---|",
    sep = "\n"
  )
)

for (i in seq_len(nrow(single_standard))) {
  row <- single_standard[i, , drop = FALSE]
  standard_lines <- c(
    standard_lines,
    paste0(
      "| ", row$Analysis_Name,
      " | ", row$Method_Label,
      " | ", row$ATF3_Rank,
      " | ", signif(row$ATF3_Method_Score, 4),
      " | ", signif(row$ATF3_Method_P_Value, 4),
      " | ", signif(row$ATF3_Method_Adjusted_P_Value, 4),
      " | ", row$ATF3_Method_Direction,
      " | ", row$ATF3_Top1,
      " |"
    )
  )
}

markdown_lines <- c(
  standard_lines,
  "",
  "## Best intersection schemes",
  "",
  paste0(
    "- Best using existing official consensus sorting: ",
    best_official$Combination_Label,
    " / sort=", best_official$Sort_Rule,
    "; ATF3 top1 in ", best_official$ATF3_Top1_Count,
    " analyses: ", best_official$ATF3_Top1_Analyses,
    "."
  ),
  paste0(
    "- Best across all tested sorting rules: ",
    best_any$Combination_Label,
    " / sort=", best_any$Sort_Rule,
    "; ATF3 top1 in ", best_any$ATF3_Top1_Count,
    " analyses: ", best_any$ATF3_Top1_Analyses,
    "."
  ),
  "",
  "## Reproducible ranking definitions",
  "",
  "- official_mean_rank: intersect TFs present in every selected method, then sort by Mean_Selected_Rank ascending, Best_Selected_Rank ascending, Source_Method_Count descending, CheA3_Library_Count descending, TF symbol.",
  "- best_rank_first: same intersection, but prioritize the best rank achieved in any selected method before the mean rank.",
  "- chea3_integrated_rank_first: same intersection, but prioritize ChEA3_Integrated_TopRank before the selected-method mean rank.",
  "- chea3_library_then_integrated_rank: prioritize broader ChEA3 library support count, then ChEA3 integrated rank, then selected-method mean rank.",
  "- P-value rules use available method P values only; ChEA3 and VIPER often do not provide comparable P values in the integrated table, so these are supplementary rather than preferred.",
  "",
  "## Output files",
  "",
  "- atf3_single_method_rank_by_sort_rule.csv",
  "- atf3_single_method_top1_summary.csv",
  "- intersection_combination_atf3_rank_detail.csv",
  "- intersection_combination_search_summary.csv",
  "- best_any_sort_rule_top10_candidates.csv",
  "- best_official_mean_rank_top10_candidates.csv",
  ""
)

writeLines(markdown_lines, file.path(OUTPUT_ROOT, "ATF3_TF_rank_intersection_search.md"))

cat("ATF3 TF rank/intersection search saved to: ", OUTPUT_ROOT, "\n", sep = "")
cat("Best official mean-rank scheme: ", best_official$Combination_Label, " | top1 analyses=", best_official$ATF3_Top1_Count, " | ", best_official$ATF3_Top1_Analyses, "\n", sep = "")
cat("Best any-rule scheme: ", best_any$Combination_Label, " | sort=", best_any$Sort_Rule, " | top1 analyses=", best_any$ATF3_Top1_Count, " | ", best_any$ATF3_Top1_Analyses, "\n", sep = "")
