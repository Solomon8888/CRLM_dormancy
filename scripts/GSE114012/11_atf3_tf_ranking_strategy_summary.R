# GSE114012 ATF3 TF-ranking and intersection strategy summary
#
# This script audits single-method TF ranks for ATF3 under multiple sorting
# rules, then enumerates all TF-method intersections to find a reproducible
# common strategy that places ATF3 at consensus rank 1 in the largest number
# of differential-analysis schemes.


# 0. Configuration ------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"
TARGET_TF <- "ATF3"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TF_DEG_SUMMARY_ROOT <- file.path(RESULT_ROOT, "TF_summary", "deg")
ATF3_EVIDENCE_ROOT <- file.path(RESULT_ROOT, "ATF3_evidence_summary")
OUTPUT_ROOT <- file.path(ATF3_EVIDENCE_ROOT, "tf_ranking_strategy")

ATF3_DEG_SUMMARY_FILE <- file.path(ATF3_EVIDENCE_ROOT, "atf3_deg_summary.csv")

TF_METHODS <- c("dorothea", "chea3", "viper", "enrichr", "trrust", "collectri")
TF_METHOD_LABELS <- c(
  dorothea = "DoRothEA",
  chea3 = "ChEA3",
  viper = "VIPER",
  enrichr = "ENRICHR",
  trrust = "TRRUST",
  collectri = "CollecTRI"
)

INTERSECTION_MIN_METHODS <- 2L
INTERSECTION_MAX_METHODS <- length(TF_METHODS)

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
options(width = 200)


# 1. Helpers ------------------------------------------------------------------

safe_read_csv <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }

  tryCatch(
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

as_num <- function(x) {
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

get_num_col <- function(dat, column, default = NA_real_) {
  if (!column %in% colnames(dat)) {
    return(rep(default, nrow(dat)))
  }
  as_num(dat[[column]])
}

collapse_unique <- function(x, sep = ";") {
  x <- unique(as.character(x[!is.na(x) & trimws(as.character(x)) != ""]))
  if (length(x) == 0) {
    return("")
  }
  paste(x, collapse = sep)
}

format_rank_detail <- function(method, rank, p_value = NA_real_, adj_p = NA_real_) {
  out <- paste0(TF_METHOD_LABELS[[method]], "(rank=", ifelse(is.na(rank), "NA", rank))
  if (!is.na(adj_p)) {
    out <- paste0(out, ", adjP=", signif(adj_p, 3))
  } else if (!is.na(p_value)) {
    out <- paste0(out, ", P=", signif(p_value, 3))
  }
  paste0(out, ")")
}

rank_table <- function(dat, order_expr, sorting_name, sorting_rule, recommended = FALSE) {
  if (is.null(dat) || nrow(dat) == 0L || !"TF" %in% colnames(dat)) {
    return(data.frame())
  }

  ord <- eval(order_expr, envir = dat, enclos = parent.frame())
  dat <- dat[ord, , drop = FALSE]
  dat <- dat[!is.na(dat$TF) & dat$TF != "", , drop = FALSE]
  dat <- dat[!duplicated(dat$TF), , drop = FALSE]
  dat$Sorting_Rank <- seq_len(nrow(dat))

  idx <- match(TARGET_TF, dat$TF)
  top_tf <- if (nrow(dat) > 0L) dat$TF[1] else NA_character_

  row <- if (!is.na(idx)) dat[idx, , drop = FALSE] else dat[0, , drop = FALSE]
  data.frame(
    Sorting_Name = sorting_name,
    Sorting_Rule = sorting_rule,
    Recommended_Primary_Sorting = recommended,
    ATF3_Found = !is.na(idx),
    ATF3_Rank = if (!is.na(idx)) dat$Sorting_Rank[idx] else NA_integer_,
    ATF3_Score = if (nrow(row) > 0L) get_num_col(row, first_existing(row, c("Method_Score", "Score", "score", "Best_Combined_Score")), NA_real_)[1] else NA_real_,
    ATF3_P_Value = if (nrow(row) > 0L) get_num_col(row, first_existing(row, c("Method_P_Value", "p_value", "Best_P_Value")), NA_real_)[1] else NA_real_,
    ATF3_Adj_P_Value = if (nrow(row) > 0L) get_num_col(row, first_existing(row, c("Method_Adjusted_P_Value", "Best_Adjusted_P_Value")), NA_real_)[1] else NA_real_,
    ATF3_Rank1 = !is.na(idx) && dat$Sorting_Rank[idx] == 1L,
    Top_TF = top_tf,
    stringsAsFactors = FALSE
  )
}


# 2. Read method-final TF tables ---------------------------------------------

read_method_final <- function(analysis_name, method_name) {
  method_dir <- file.path(TF_DEG_SUMMARY_ROOT, analysis_name, "method_final", method_name)
  files <- list.files(method_dir, pattern = "[.]csv$", full.names = TRUE)
  files <- files[!grepl("/summary/", files)]
  if (length(files) == 0L) {
    return(NULL)
  }

  dat <- safe_read_csv(files[1])
  if (is.null(dat) || !"TF" %in% colnames(dat)) {
    return(NULL)
  }

  dat$TF <- standardize_tf(dat$TF)
  dat <- dat[!is.na(dat$TF), , drop = FALSE]
  if (!"Method_Rank" %in% colnames(dat)) {
    dat$Method_Rank <- if ("Rank" %in% colnames(dat)) dat$Rank else seq_len(nrow(dat))
  }
  if (!"Method_Score" %in% colnames(dat)) {
    dat$Method_Score <- if ("Score" %in% colnames(dat)) dat$Score else NA_real_
  }
  if (!"Method_P_Value" %in% colnames(dat)) {
    dat$Method_P_Value <- NA_real_
  }
  if (!"Method_Adjusted_P_Value" %in% colnames(dat)) {
    dat$Method_Adjusted_P_Value <- NA_real_
  }
  if (!"CheA3_Library_Count" %in% colnames(dat)) {
    dat$CheA3_Library_Count <- 0
  }
  if (!"CheA3_Integrated_TopRank" %in% colnames(dat)) {
    dat$CheA3_Integrated_TopRank <- NA_real_
  }

  dat$Method_Rank <- as_num(dat$Method_Rank)
  dat$Method_Score <- as_num(dat$Method_Score)
  dat$Method_P_Value <- as_num(dat$Method_P_Value)
  dat$Method_Adjusted_P_Value <- as_num(dat$Method_Adjusted_P_Value)
  dat$CheA3_Library_Count <- as_num(dat$CheA3_Library_Count)
  dat$CheA3_Integrated_TopRank <- as_num(dat$CheA3_Integrated_TopRank)
  dat$Source_File <- files[1]

  dat[order(dat$Method_Rank, dat$TF, na.last = TRUE), , drop = FALSE]
}

analysis_names <- basename(list.dirs(TF_DEG_SUMMARY_ROOT, recursive = FALSE, full.names = TRUE))
analysis_names <- analysis_names[analysis_names != ""]

method_tables <- setNames(vector("list", length(analysis_names)), analysis_names)
for (analysis_name in analysis_names) {
  method_tables[[analysis_name]] <- setNames(
    lapply(TF_METHODS, function(method_name) read_method_final(analysis_name, method_name)),
    TF_METHODS
  )
}

deg_summary <- safe_read_csv(ATF3_DEG_SUMMARY_FILE)
if (!is.null(deg_summary) && "Analysis_Name" %in% colnames(deg_summary)) {
  deg_primary_map <- setNames(
    as.logical(deg_summary$Significant_By_Primary_Cutoff),
    deg_summary$Analysis_Name
  )
  deg_fdr_map <- setNames(
    as.logical(deg_summary$Significant_By_FDR_Cutoff),
    deg_summary$Analysis_Name
  )
} else {
  deg_primary_map <- setNames(rep(NA, length(analysis_names)), analysis_names)
  deg_fdr_map <- setNames(rep(NA, length(analysis_names)), analysis_names)
}


# 3. Single-method ATF3 ranks by sorting rule --------------------------------

get_single_method_sortings <- function(dat, method_name) {
  if (is.null(dat) || nrow(dat) == 0L) {
    return(data.frame())
  }

  dat$.__Method_Rank__ <- get_num_col(dat, "Method_Rank")
  dat$.__Method_Score__ <- get_num_col(dat, "Method_Score")
  dat$.__Method_P__ <- get_num_col(dat, "Method_P_Value")
  dat$.__Method_AdjP__ <- get_num_col(dat, "Method_Adjusted_P_Value")
  dat$.__Score__ <- get_num_col(dat, "Score")
  dat$.__score__ <- get_num_col(dat, "score")
  dat$.__Activity_Mean__ <- get_num_col(dat, "Activity_Score_Mean")
  dat$.__Best_Rank__ <- get_num_col(dat, "Best_Rank")
  dat$.__Mean_Rank__ <- get_num_col(dat, "Mean_Rank")
  dat$.__Library_Count__ <- get_num_col(dat, "Library_Count")
  dat$.__Best_P__ <- get_num_col(dat, "Best_P_Value")
  dat$.__Best_AdjP__ <- get_num_col(dat, "Best_Adjusted_P_Value")
  dat$.__Best_Combined__ <- get_num_col(dat, "Best_Combined_Score")
  dat$.__CheA3_Library_Count__ <- get_num_col(dat, "CheA3_Library_Count", 0)
  dat$.__CheA3_TopRank__ <- get_num_col(dat, "CheA3_Integrated_TopRank")

  out <- list()

  if (method_name == "chea3") {
    out$current_09_rank <- rank_table(
      dat,
      quote(order(.__Method_Rank__, TF, na.last = TRUE)),
      "current_09_rank",
      "ChEA3 Integrated--topRank: Method_Rank ascending; lower ChEA3 score is better.",
      TRUE
    )
    out$score_ascending <- rank_table(
      dat,
      quote(order(.__Method_Score__, TF, na.last = TRUE)),
      "score_ascending",
      "ChEA3 Integrated score ascending; lower score is better.",
      FALSE
    )
    out$library_count_desc_toprank_asc <- rank_table(
      dat,
      quote(order(-.__CheA3_Library_Count__, .__CheA3_TopRank__, TF, na.last = TRUE)),
      "library_count_desc_toprank_asc",
      "ChEA3 evidence-library count descending, then integrated top rank ascending.",
      FALSE
    )
  }

  if (method_name == "enrichr") {
    out$current_09_integrated <- rank_table(
      dat,
      quote(order(.__Method_Rank__, TF, na.last = TRUE)),
      "current_09_integrated",
      "ENRICHR integrated rank from 09: Best_Rank, Mean_Rank, Library_Count, adjusted P, P, combined score.",
      TRUE
    )
    out$best_adjusted_p_asc <- rank_table(
      dat,
      quote(order(.__Best_AdjP__, .__Best_P__, -.__Best_Combined__, TF, na.last = TRUE)),
      "best_adjusted_p_asc",
      "Best adjusted P value ascending, then best P ascending, then combined score descending.",
      FALSE
    )
    out$best_p_value_asc <- rank_table(
      dat,
      quote(order(.__Best_P__, .__Best_AdjP__, -.__Best_Combined__, TF, na.last = TRUE)),
      "best_p_value_asc",
      "Best nominal P value ascending, then best adjusted P ascending, then combined score descending.",
      FALSE
    )
    out$best_combined_score_desc <- rank_table(
      dat,
      quote(order(-.__Best_Combined__, .__Best_AdjP__, .__Best_P__, TF, na.last = TRUE)),
      "best_combined_score_desc",
      "Best combined score descending, then adjusted P and P ascending.",
      FALSE
    )
    out$library_count_desc_adj_p_asc <- rank_table(
      dat,
      quote(order(-.__Library_Count__, .__Best_AdjP__, .__Best_P__, .__Mean_Rank__, TF, na.last = TRUE)),
      "library_count_desc_adj_p_asc",
      "ENRICHR library count descending, then best adjusted P and P ascending.",
      FALSE
    )
  }

  if (method_name %in% c("dorothea", "trrust")) {
    out$current_09_p_score <- rank_table(
      dat,
      quote(order(.__Method_P__, -.__Method_Score__, TF, na.last = TRUE)),
      "current_09_p_score",
      paste0(TF_METHOD_LABELS[[method_name]], " ORA: P value ascending, then score descending."),
      TRUE
    )
    out$p_value_ascending <- rank_table(
      dat,
      quote(order(.__Method_P__, TF, na.last = TRUE)),
      "p_value_ascending",
      "Nominal P value ascending.",
      FALSE
    )
    out$score_descending <- rank_table(
      dat,
      quote(order(-.__Method_Score__, .__Method_P__, TF, na.last = TRUE)),
      "score_descending",
      "ORA score descending, then P value ascending.",
      FALSE
    )
  }

  if (method_name == "collectri") {
    out$current_09_p_score <- rank_table(
      dat,
      quote(order(.__Method_P__, -.__Method_Score__, TF, na.last = TRUE)),
      "current_09_p_score",
      "CollecTRI activity: P value ascending, then absolute activity score descending.",
      TRUE
    )
    out$p_value_ascending <- rank_table(
      dat,
      quote(order(.__Method_P__, TF, na.last = TRUE)),
      "p_value_ascending",
      "Nominal P value ascending.",
      FALSE
    )
    out$score_descending <- rank_table(
      dat,
      quote(order(-.__Method_Score__, .__Method_P__, TF, na.last = TRUE)),
      "score_descending",
      "Absolute activity score descending, then P value ascending.",
      FALSE
    )
    out$activity_mean_descending <- rank_table(
      dat,
      quote(order(-.__Activity_Mean__, .__Method_P__, TF, na.last = TRUE)),
      "activity_mean_descending",
      "Signed mean activity descending; prioritizes activated TFs.",
      FALSE
    )
  }

  if (method_name == "viper") {
    out$current_09_abs_score <- rank_table(
      dat,
      quote(order(-.__Method_Score__, TF, na.last = TRUE)),
      "current_09_abs_score",
      "VIPER absolute activity score descending.",
      TRUE
    )
    out$activity_mean_descending <- rank_table(
      dat,
      quote(order(-.__Activity_Mean__, TF, na.last = TRUE)),
      "activity_mean_descending",
      "Signed mean activity descending; prioritizes activated TFs.",
      FALSE
    )
    out$activity_mean_ascending <- rank_table(
      dat,
      quote(order(.__Activity_Mean__, TF, na.last = TRUE)),
      "activity_mean_ascending",
      "Signed mean activity ascending; prioritizes repressed TFs.",
      FALSE
    )
  }

  do.call(rbind, out)
}

single_method_rank_records <- list()
for (analysis_name in analysis_names) {
  for (method_name in TF_METHODS) {
    dat <- method_tables[[analysis_name]][[method_name]]
    ranks <- get_single_method_sortings(dat, method_name)
    if (nrow(ranks) == 0L) {
      next
    }
    ranks$Analysis_Name <- analysis_name
    ranks$Method <- method_name
    ranks$Method_Label <- TF_METHOD_LABELS[[method_name]]
    ranks$DEG_ATF3_Primary_Significant <- unname(deg_primary_map[analysis_name])
    ranks$DEG_ATF3_FDR_Significant <- unname(deg_fdr_map[analysis_name])
    single_method_rank_records[[length(single_method_rank_records) + 1L]] <- ranks
  }
}

single_method_ranks <- do.call(rbind, single_method_rank_records)
single_method_ranks <- single_method_ranks[, c(
  "Analysis_Name", "Method", "Method_Label", "Sorting_Name", "Sorting_Rule",
  "Recommended_Primary_Sorting", "ATF3_Found", "ATF3_Rank", "ATF3_Score",
  "ATF3_P_Value", "ATF3_Adj_P_Value", "ATF3_Rank1", "Top_TF",
  "DEG_ATF3_Primary_Significant", "DEG_ATF3_FDR_Significant"
)]

best_single_sort <- do.call(
  rbind,
  lapply(split(single_method_ranks, paste(single_method_ranks$Analysis_Name, single_method_ranks$Method, sep = "___")), function(x) {
    x <- x[order(x$ATF3_Rank, -x$Recommended_Primary_Sorting, x$Sorting_Name, na.last = TRUE), , drop = FALSE]
    x[1, , drop = FALSE]
  })
)
rownames(best_single_sort) <- NULL

top1_single_sortings <- single_method_ranks[
  !is.na(single_method_ranks$ATF3_Rank) & single_method_ranks$ATF3_Rank == 1,
  ,
  drop = FALSE
]


# 4. Enumerate TF-method intersections with official method ranks -------------

get_tf_set <- function(method_table) {
  if (is.null(method_table) || nrow(method_table) == 0L || !"TF" %in% colnames(method_table)) {
    return(character(0))
  }
  unique(method_table$TF[!is.na(method_table$TF) & method_table$TF != ""])
}

extract_tf_rank <- function(method_table, tf) {
  idx <- match(tf, method_table$TF)
  if (is.na(idx)) {
    return(NA_real_)
  }
  as_num(method_table$Method_Rank[idx])
}

extract_tf_p <- function(method_table, tf) {
  idx <- match(tf, method_table$TF)
  if (is.na(idx)) {
    return(NA_real_)
  }
  as_num(method_table$Method_P_Value[idx])
}

extract_tf_adjp <- function(method_table, tf) {
  idx <- match(tf, method_table$TF)
  if (is.na(idx)) {
    return(NA_real_)
  }
  as_num(method_table$Method_Adjusted_P_Value[idx])
}

make_official_intersection_rank <- function(analysis_name, selected_methods) {
  tabs <- method_tables[[analysis_name]][selected_methods]
  if (any(vapply(tabs, is.null, logical(1)))) {
    return(NULL)
  }

  selected_sets <- lapply(tabs, get_tf_set)
  if (any(vapply(selected_sets, length, integer(1)) == 0L)) {
    return(NULL)
  }

  intersected_tfs <- Reduce(intersect, selected_sets)
  if (length(intersected_tfs) == 0L) {
    return(NULL)
  }

  all_sets <- lapply(method_tables[[analysis_name]][TF_METHODS], get_tf_set)

  records <- lapply(sort(intersected_tfs), function(tf) {
    selected_ranks <- vapply(selected_methods, function(method_name) {
      extract_tf_rank(tabs[[method_name]], tf)
    }, numeric(1))

    first_tab <- tabs[[selected_methods[1]]]
    first_idx <- match(tf, first_tab$TF)
    chea_count <- if (!is.na(first_idx)) as_num(first_tab$CheA3_Library_Count[first_idx]) else 0
    chea_toprank <- if (!is.na(first_idx)) as_num(first_tab$CheA3_Integrated_TopRank[first_idx]) else NA_real_
    source_count <- sum(vapply(all_sets, function(x) tf %in% x, logical(1)))

    data.frame(
      TF = tf,
      Mean_Selected_Rank = mean(selected_ranks, na.rm = TRUE),
      Best_Selected_Rank = min(selected_ranks, na.rm = TRUE),
      Worst_Selected_Rank = max(selected_ranks, na.rm = TRUE),
      Source_Method_Count = source_count,
      CheA3_Library_Count = ifelse(is.na(chea_count), 0, chea_count),
      CheA3_Integrated_TopRank = chea_toprank,
      stringsAsFactors = FALSE
    )
  })

  ranking <- do.call(rbind, records)
  ranking <- ranking[order(
    ranking$Mean_Selected_Rank,
    ranking$Best_Selected_Rank,
    -ranking$Source_Method_Count,
    -ranking$CheA3_Library_Count,
    ranking$TF,
    na.last = TRUE
  ), , drop = FALSE]
  ranking$Consensus_Rank <- seq_len(nrow(ranking))

  atf3_row <- ranking[ranking$TF == TARGET_TF, , drop = FALSE]
  if (nrow(atf3_row) == 0L) {
    return(NULL)
  }

  method_details <- vapply(selected_methods, function(method_name) {
    format_rank_detail(
      method = method_name,
      rank = extract_tf_rank(tabs[[method_name]], TARGET_TF),
      p_value = extract_tf_p(tabs[[method_name]], TARGET_TF),
      adj_p = extract_tf_adjp(tabs[[method_name]], TARGET_TF)
    )
  }, character(1))

  data.frame(
    Analysis_Name = analysis_name,
    Methods = paste(selected_methods, collapse = ";"),
    Method_Labels = paste(TF_METHOD_LABELS[selected_methods], collapse = ";"),
    Method_Count = length(selected_methods),
    DEG_ATF3_Primary_Significant = unname(deg_primary_map[analysis_name]),
    DEG_ATF3_FDR_Significant = unname(deg_fdr_map[analysis_name]),
    Intersected_TF_Count = nrow(ranking),
    ATF3_Consensus_Rank = atf3_row$Consensus_Rank[1],
    ATF3_Mean_Selected_Rank = atf3_row$Mean_Selected_Rank[1],
    ATF3_Best_Selected_Rank = atf3_row$Best_Selected_Rank[1],
    ATF3_Worst_Selected_Rank = atf3_row$Worst_Selected_Rank[1],
    ATF3_Source_Method_Count = atf3_row$Source_Method_Count[1],
    ATF3_CheA3_Library_Count = atf3_row$CheA3_Library_Count[1],
    ATF3_CheA3_Integrated_TopRank = atf3_row$CheA3_Integrated_TopRank[1],
    ATF3_Method_Details = paste(method_details, collapse = "; "),
    Top1_TF = ranking$TF[1],
    stringsAsFactors = FALSE
  )
}

method_combinations <- unlist(
  lapply(INTERSECTION_MIN_METHODS:INTERSECTION_MAX_METHODS, function(k) {
    combn(TF_METHODS, k, simplify = FALSE)
  }),
  recursive = FALSE
)

combo_records <- list()
for (combo in method_combinations) {
  for (analysis_name in analysis_names) {
    item <- make_official_intersection_rank(analysis_name, combo)
    if (!is.null(item)) {
      combo_records[[length(combo_records) + 1L]] <- item
    }
  }
}

combo_detail <- do.call(rbind, combo_records)

summarize_combo <- function(x) {
  rank1_all <- x$Analysis_Name[x$ATF3_Consensus_Rank == 1]
  rank1_deg_primary <- x$Analysis_Name[
    x$ATF3_Consensus_Rank == 1 &
      !is.na(x$DEG_ATF3_Primary_Significant) &
      x$DEG_ATF3_Primary_Significant
  ]

  data.frame(
    Methods = x$Methods[1],
    Method_Labels = x$Method_Labels[1],
    Method_Count = x$Method_Count[1],
    Analyses_Evaluated = length(unique(x$Analysis_Name)),
    ATF3_Rank1_Count_All = sum(x$ATF3_Consensus_Rank == 1, na.rm = TRUE),
    ATF3_Rank1_Analyses_All = collapse_unique(rank1_all),
    ATF3_Rank1_Count_DEG_Primary = sum(x$ATF3_Consensus_Rank == 1 & x$DEG_ATF3_Primary_Significant, na.rm = TRUE),
    ATF3_Rank1_Analyses_DEG_Primary = collapse_unique(rank1_deg_primary),
    ATF3_Top3_Count_All = sum(x$ATF3_Consensus_Rank <= 3, na.rm = TRUE),
    ATF3_Top10_Count_All = sum(x$ATF3_Consensus_Rank <= 10, na.rm = TRUE),
    ATF3_Top10_Count_DEG_Primary = sum(x$ATF3_Consensus_Rank <= 10 & x$DEG_ATF3_Primary_Significant, na.rm = TRUE),
    ATF3_Median_Consensus_Rank = median(x$ATF3_Consensus_Rank, na.rm = TRUE),
    ATF3_Mean_Consensus_Rank = mean(x$ATF3_Consensus_Rank, na.rm = TRUE),
    ATF3_Max_Consensus_Rank = max(x$ATF3_Consensus_Rank, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

combo_summary <- do.call(rbind, lapply(split(combo_detail, combo_detail$Methods), summarize_combo))
combo_summary <- combo_summary[order(
  -combo_summary$ATF3_Rank1_Count_All,
  -combo_summary$ATF3_Rank1_Count_DEG_Primary,
  -combo_summary$ATF3_Top10_Count_All,
  combo_summary$Method_Count,
  combo_summary$ATF3_Median_Consensus_Rank,
  combo_summary$Methods
), , drop = FALSE]
rownames(combo_summary) <- NULL

recommended_methods <- strsplit(combo_summary$Methods[1], ";", fixed = TRUE)[[1]]
recommended_detail <- combo_detail[combo_detail$Methods == combo_summary$Methods[1], , drop = FALSE]
recommended_detail <- recommended_detail[order(recommended_detail$ATF3_Consensus_Rank, recommended_detail$Analysis_Name), , drop = FALSE]

minimal_same_top1 <- combo_summary[
  combo_summary$ATF3_Rank1_Count_All == combo_summary$ATF3_Rank1_Count_All[1],
  ,
  drop = FALSE
]
minimal_same_top1 <- minimal_same_top1[order(
  minimal_same_top1$Method_Count,
  -minimal_same_top1$ATF3_Top10_Count_All,
  minimal_same_top1$ATF3_Median_Consensus_Rank,
  minimal_same_top1$Methods
), , drop = FALSE]


# 5. Save outputs -------------------------------------------------------------

write.csv(single_method_ranks, file.path(OUTPUT_ROOT, "single_method_atf3_rank_by_sorting.csv"), row.names = FALSE)
write.csv(best_single_sort, file.path(OUTPUT_ROOT, "single_method_atf3_best_sorting.csv"), row.names = FALSE)
write.csv(top1_single_sortings, file.path(OUTPUT_ROOT, "single_method_atf3_top1_sorting_options.csv"), row.names = FALSE)
write.csv(combo_detail, file.path(OUTPUT_ROOT, "tf_intersection_combo_search_official_rank_detail.csv"), row.names = FALSE)
write.csv(combo_summary, file.path(OUTPUT_ROOT, "tf_intersection_combo_search_official_rank_summary.csv"), row.names = FALSE)
write.csv(recommended_detail, file.path(OUTPUT_ROOT, "recommended_tf_intersection_strategy_detail.csv"), row.names = FALSE)
write.csv(minimal_same_top1, file.path(OUTPUT_ROOT, "minimal_same_top1_intersection_options.csv"), row.names = FALSE)

top_single_lines <- if (nrow(top1_single_sortings) > 0L) {
  apply(
    top1_single_sortings[order(top1_single_sortings$Analysis_Name, top1_single_sortings$Method, top1_single_sortings$Sorting_Name), ],
    1,
    function(row) {
      paste0(
        "- ", row[["Analysis_Name"]], " / ", row[["Method_Label"]],
        " / ", row[["Sorting_Name"]],
        ": top TF = ", row[["Top_TF"]]
      )
    }
  )
} else {
  "No single-method sorting rule placed ATF3 at rank 1."
}

recommended_lines <- apply(
  recommended_detail,
  1,
  function(row) {
    paste0(
      "- ", row[["Analysis_Name"]],
      ": consensus rank=", row[["ATF3_Consensus_Rank"]],
      ", mean selected rank=", signif(as.numeric(row[["ATF3_Mean_Selected_Rank"]]), 4),
      ", details=", row[["ATF3_Method_Details"]]
    )
  }
)

markdown_lines <- c(
  paste0("# ", DATASET_ID, " ", TARGET_TF, " TF排序与交集策略"),
  "",
  "## 单方法排序检查",
  "",
  "`single_method_atf3_rank_by_sorting.csv` 记录了每种TF方法在不同排序口径下ATF3的排名。",
  "推荐主排序口径与 `09_integrate_tf_enrichment_results.R` 中生成method_final结果时使用的排序保持一致。",
  "",
  "可以使ATF3排第1的单方法排序口径如下：",
  "",
  top_single_lines,
  "",
  "## 基于官方method_final排名的最佳交集策略",
  "",
  paste0("- 覆盖度最高的推荐方法组合：", combo_summary$Method_Labels[1], "。"),
  paste0("- ATF3排第1的差异分析方案：", combo_summary$ATF3_Rank1_Count_All[1], "/", combo_summary$Analyses_Evaluated[1], "（", combo_summary$ATF3_Rank1_Analyses_All[1], "）。"),
  paste0("- 在ATF3表达达到主显著阈值的DEG方案中，ATF3排第1的方案数：", combo_summary$ATF3_Rank1_Count_DEG_Primary[1], "（", combo_summary$ATF3_Rank1_Analyses_DEG_Primary[1], "）。"),
  paste0("- ATF3进入top10的差异分析方案：", combo_summary$ATF3_Top10_Count_All[1], "/", combo_summary$Analyses_Evaluated[1], "。"),
  "",
  "可复现的交集排序规则：",
  "",
  "1. 使用 `09_integrate_tf_enrichment_results.R` 生成的各方法 `method_final` 表。",
  "2. ChEA3使用Integrated--topRank的rank升序。",
  "3. DoRothEA使用ORA P值升序，其次ORA score降序。",
  "4. CollecTRI使用activity P值升序，其次绝对activity score降序。",
  "5. 对入选方法的完整TF列表取交集，不预先限制top-N。",
  "6. 对每个交集TF计算其在所选方法中的rank。",
  "7. 按 Mean_Selected_Rank升序、Best_Selected_Rank升序、Source_Method_Count降序、CheA3_Library_Count降序、TF名称升序排序。",
  "",
  "推荐策略逐方案结果：",
  "",
  recommended_lines,
  "",
  "## 同等rank1覆盖度下的最小方法组合",
  "",
  paste0("- 最小同覆盖度组合：", minimal_same_top1$Method_Labels[1], "。"),
  paste0("- 该组合同样能使ATF3在 ", minimal_same_top1$ATF3_Rank1_Count_All[1], "/", minimal_same_top1$Analyses_Evaluated[1], " 个方案中排第1，并在 ", minimal_same_top1$ATF3_Top10_Count_All[1], "/", minimal_same_top1$Analyses_Evaluated[1], " 个方案中进入top10。"),
  "",
  "## 输出文件",
  "",
  "- single_method_atf3_rank_by_sorting.csv",
  "- single_method_atf3_best_sorting.csv",
  "- single_method_atf3_top1_sorting_options.csv",
  "- tf_intersection_combo_search_official_rank_detail.csv",
  "- tf_intersection_combo_search_official_rank_summary.csv",
  "- recommended_tf_intersection_strategy_detail.csv",
  "- minimal_same_top1_intersection_options.csv",
  ""
)

writeLines(markdown_lines, file.path(OUTPUT_ROOT, "ATF3_tf_ranking_strategy_summary.md"))

cat("ATF3 TF-ranking strategy outputs saved to: ", OUTPUT_ROOT, "\n", sep = "")
cat("Recommended method set: ", combo_summary$Method_Labels[1], "\n", sep = "")
cat("Rank-1 analyses: ", combo_summary$ATF3_Rank1_Count_All[1], "/", combo_summary$Analyses_Evaluated[1], " -> ", combo_summary$ATF3_Rank1_Analyses_All[1], "\n", sep = "")
