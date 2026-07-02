# GSE114012 ATF3 evidence summary
#
# 整理每套analysis设计中ATF3的差异表达证据、TF富集证据和TF交集证据。


# 0. Configuration ------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"
TARGET_GENE <- "ATF3"

CLINICAL_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_clinical_edit.csv"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
DEG_ROOT <- file.path(RESULT_ROOT, "tables")
TF_SUMMARY_ROOT <- file.path(RESULT_ROOT, "TF_summary")
OUTPUT_ROOT <- file.path(RESULT_ROOT, "ATF3_evidence_summary")

P_VALUE_CUTOFF <- 0.05
P_VALUE_COLUMN <- "P.Value"
ADJ_P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5
TF_TOP_RANK_CUTOFF <- 10

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

options(width = 200)


# 1. Helper functions ---------------------------------------------------------

source(FUNCTION_FILE)

safe_read_csv <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }

  tryCatch(
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

as_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% colnames(x)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

collapse_unique <- function(x, sep = ";") {
  x <- unique(as.character(x[!is.na(x) & trimws(as.character(x)) != ""]))
  if (length(x) == 0) {
    return("")
  }
  paste(x, collapse = sep)
}

cell_chr <- function(row, candidates) {
  col <- first_existing(row, candidates)
  if (is.na(col)) {
    return(NA_character_)
  }
  as.character(row[[col]][1])
}

cell_num <- function(row, candidates) {
  as_num(cell_chr(row, candidates))
}

clean_path <- function(x) {
  gsub("\\\\", "/", x)
}

extract_tfsummary_meta <- function(file, expected_section = NULL) {
  parts <- strsplit(clean_path(file), "/", fixed = TRUE)[[1]]
  idx <- match("TF_summary", parts)
  if (is.na(idx) || length(parts) < idx + 3) {
    return(NULL)
  }

  section <- parts[idx + 3]
  if (!is.null(expected_section) && !identical(section, expected_section)) {
    return(NULL)
  }

  meta <- list(
    Input_Type = toupper(parts[idx + 1]),
    TF_Analysis_Name = parts[idx + 2]
  )

  if (identical(section, "method_final")) {
    meta$Method <- parts[idx + 4]
  }

  if (identical(section, "intersections")) {
    meta$Intersection_Name <- parts[idx + 4]
  }

  meta
}

extract_rank_from_method_row <- function(row) {
  rank_value <- cell_num(row, c("Rank", "Method_Rank", "Best_Rank", "Consensus_Rank"))
  if (!is.na(rank_value)) {
    return(rank_value)
  }
  NA_real_
}

format_method_detail <- function(df) {
  if (nrow(df) == 0) {
    return("")
  }

  details <- mapply(
    function(method, rank, p, adjp) {
      item <- paste0(method, "(rank=", ifelse(is.na(rank), "NA", rank))
      if (!is.na(adjp)) {
        item <- paste0(item, ",adjP=", signif(adjp, 3))
      } else if (!is.na(p)) {
        item <- paste0(item, ",P=", signif(p, 3))
      }
      paste0(item, ")")
    },
    df$Method,
    df$Rank,
    df$P_Value,
    df$Adj_P_Value,
    USE.NAMES = FALSE
  )

  paste(details, collapse = ";")
}

format_intersection_detail <- function(df) {
  if (nrow(df) == 0) {
    return("")
  }

  details <- mapply(
    function(name, rank, methods) {
      paste0(name, "(rank=", rank, ",methods=", methods, ")")
    },
    df$Intersection_Name,
    df$Consensus_Rank,
    df$Required_Methods,
    USE.NAMES = FALSE
  )

  paste(details, collapse = ";")
}


# 2. Analysis design summary --------------------------------------------------

clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

analysis_designs <- get_analysis_designs(clinical_data)

design_summary <- do.call(
  rbind,
  lapply(seq_len(nrow(analysis_designs)), function(i) {
    design_samples <- prepare_design_samples(
      sample_info = clinical_data,
      group_column_index = analysis_designs$Column_Index[i],
      experiment_group = analysis_designs$Experiment_Group[i]
    )

    group_list <- design_samples$group_list
    control_group <- design_samples$control_group
    experiment_group <- analysis_designs$Experiment_Group[i]

    data.frame(
      Dataset = DATASET_ID,
      Analysis_Order = analysis_designs$Analysis_Order[i],
      Analysis_Name = analysis_designs$Analysis_Name[i],
      Design_Column = analysis_designs$Column_Name[i],
      Experiment_Group = experiment_group,
      Control_Group = control_group,
      Contrast = paste0(experiment_group, "_vs_", control_group),
      Samples_Used = length(group_list),
      Experiment_N = sum(group_list == experiment_group),
      Control_N = sum(group_list == control_group),
      stringsAsFactors = FALSE
    )
  })
)


# 3. ATF3 differential expression summary ------------------------------------

extract_atf3_deg <- function(analysis_name) {
  all_file <- file.path(DEG_ROOT, analysis_name, "DEG", "all_genes.csv")
  sig_file <- file.path(DEG_ROOT, analysis_name, "DEG", "significant_genes.csv")
  summary_file <- file.path(DEG_ROOT, analysis_name, "DEG", "summary.csv")

  all_genes <- safe_read_csv(all_file)
  sig_genes <- safe_read_csv(sig_file)
  deg_summary <- safe_read_csv(summary_file)

  if (is.null(all_genes) || !"Symbol" %in% colnames(all_genes)) {
    return(data.frame(
      Analysis_Name = analysis_name,
      ATF3_Found_In_DEG_Table = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  idx <- which(toupper(trimws(as.character(all_genes$Symbol))) == toupper(TARGET_GENE))

  if (length(idx) == 0) {
    return(data.frame(
      Analysis_Name = analysis_name,
      ATF3_Found_In_DEG_Table = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  row <- all_genes[idx[which.min(all_genes$P.Value[idx])], , drop = FALSE]
  atf3_in_sig <- FALSE
  if (!is.null(sig_genes) && "Symbol" %in% colnames(sig_genes)) {
    atf3_in_sig <- any(toupper(trimws(as.character(sig_genes$Symbol))) == toupper(TARGET_GENE))
  }

  logfc <- cell_num(row, "logFC")
  pvalue <- cell_num(row, "P.Value")
  adjp <- cell_num(row, "adj.P.Val")

  data.frame(
    Analysis_Name = analysis_name,
    ATF3_Found_In_DEG_Table = TRUE,
    GeneID = cell_chr(row, "GeneID"),
    Symbol = cell_chr(row, "Symbol"),
    logFC = logfc,
    AveExpr = cell_num(row, "AveExpr"),
    t = cell_num(row, "t"),
    P.Value = pvalue,
    adj.P.Val = adjp,
    B = cell_num(row, "B"),
    ATF3_Direction = ifelse(is.na(logfc), NA_character_, ifelse(logfc > 0, "Up", "Down")),
    Significant_By_Primary_Cutoff = !is.na(pvalue) && !is.na(logfc) &&
      pvalue < P_VALUE_CUTOFF && abs(logfc) >= LOGFC_CUTOFF,
    Significant_By_FDR_Cutoff = !is.na(adjp) && !is.na(logfc) &&
      adjp < ADJ_P_VALUE_CUTOFF && abs(logfc) >= LOGFC_CUTOFF,
    ATF3_In_Significant_Genes_Table = atf3_in_sig,
    DEG_Total_Significant_Genes = if (!is.null(deg_summary) && "Total_Significant_Genes" %in% colnames(deg_summary)) deg_summary$Total_Significant_Genes[1] else NA,
    DEG_Up = if (!is.null(deg_summary) && "Up" %in% colnames(deg_summary)) deg_summary$Up[1] else NA,
    DEG_Down = if (!is.null(deg_summary) && "Down" %in% colnames(deg_summary)) deg_summary$Down[1] else NA,
    DEG_Result_File = all_file,
    stringsAsFactors = FALSE
  )
}

atf3_deg_summary <- do.call(
  rbind,
  lapply(design_summary$Analysis_Name, extract_atf3_deg)
)

atf3_deg_summary <- merge(
  design_summary,
  atf3_deg_summary,
  by = "Analysis_Name",
  all.x = TRUE,
  sort = FALSE
)


# 4. ATF3 TF method evidence --------------------------------------------------

method_files <- list.files(
  TF_SUMMARY_ROOT,
  pattern = "[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
method_files <- method_files[
  grepl("/method_final/", clean_path(method_files)) &
    !grepl("/method_final/summary/", clean_path(method_files))
]

extract_atf3_method <- function(file) {
  meta <- extract_tfsummary_meta(file, expected_section = "method_final")
  if (is.null(meta)) {
    return(NULL)
  }

  dat <- safe_read_csv(file)
  if (is.null(dat) || nrow(dat) == 0) {
    return(NULL)
  }

  tf_col <- first_existing(dat, c("TF", "source", "Source", "Transcription_Factor", "tf", "term", "Term"))
  if (is.na(tf_col)) {
    return(NULL)
  }

  tf_values <- toupper(trimws(as.character(dat[[tf_col]])))
  idx <- which(tf_values == toupper(TARGET_GENE))

  if (length(idx) == 0) {
    return(data.frame(
      Input_Type = meta$Input_Type,
      TF_Analysis_Name = meta$TF_Analysis_Name,
      Method = meta$Method,
      ATF3_Found = FALSE,
      Rank = NA_real_,
      Score = NA_real_,
      P_Value = NA_real_,
      Adj_P_Value = NA_real_,
      Direction = NA_character_,
      ATF3_Rank1 = FALSE,
      ATF3_TopN = FALSE,
      ATF3_P_Significant = FALSE,
      ATF3_AdjP_Significant = FALSE,
      ATF3_Method_Supported = FALSE,
      Result_File = file,
      stringsAsFactors = FALSE
    ))
  }

  dat$.__rank__ <- as_num(dat[[first_existing(dat, c("Rank", "Method_Rank", "Best_Rank"))]])
  dat$.__rank__[is.na(dat$.__rank__)] <- seq_len(nrow(dat))[is.na(dat$.__rank__)]
  best_idx <- idx[which.min(dat$.__rank__[idx])]
  row <- dat[best_idx, , drop = FALSE]

  rank <- extract_rank_from_method_row(row)
  if (is.na(rank)) {
    rank <- dat$.__rank__[best_idx]
  }
  score <- cell_num(row, c("Method_Score", "Score", "score", "Best_Combined_Score", "Activity_Score_Mean", "statistic"))
  pvalue <- cell_num(row, c("Method_P_Value", "p_value", "P.Value", "p.value", "Best_P_Value"))
  adjp <- cell_num(row, c("Method_Adjusted_P_Value", "Best_Adjusted_P_Value", "adj.P.Val", "p.adjust", "FDR"))
  direction <- cell_chr(row, c("Method_Direction", "Direction", "direction"))

  rank1 <- !is.na(rank) && rank == 1
  topn <- !is.na(rank) && rank <= TF_TOP_RANK_CUTOFF
  psig <- !is.na(pvalue) && pvalue < P_VALUE_CUTOFF
  adjpsig <- !is.na(adjp) && adjp < ADJ_P_VALUE_CUTOFF

  data.frame(
    Input_Type = meta$Input_Type,
    TF_Analysis_Name = meta$TF_Analysis_Name,
    Method = meta$Method,
    ATF3_Found = TRUE,
    Rank = rank,
    Score = score,
    P_Value = pvalue,
    Adj_P_Value = adjp,
    Direction = direction,
    ATF3_Rank1 = rank1,
    ATF3_TopN = topn,
    ATF3_P_Significant = psig,
    ATF3_AdjP_Significant = adjpsig,
    ATF3_Method_Supported = topn || psig || adjpsig,
    Result_File = file,
    stringsAsFactors = FALSE
  )
}

atf3_method_evidence <- do.call(
  rbind,
  Filter(Negate(is.null), lapply(method_files, extract_atf3_method))
)

summarize_method_evidence <- function(df) {
  if (nrow(df) == 0) {
    return(NULL)
  }

  supported_df <- df[df$ATF3_Method_Supported, , drop = FALSE]
  topn_df <- df[df$ATF3_TopN, , drop = FALSE]
  rank1_df <- df[df$ATF3_Rank1, , drop = FALSE]
  psig_df <- df[df$ATF3_P_Significant | df$ATF3_AdjP_Significant, , drop = FALSE]

  data.frame(
    Input_Type = df$Input_Type[1],
    TF_Analysis_Name = df$TF_Analysis_Name[1],
    Methods_Checked = length(unique(df$Method)),
    ATF3_Method_Found_Count = length(unique(df$Method[df$ATF3_Found])),
    ATF3_Method_Supported_Count = length(unique(supported_df$Method)),
    ATF3_Method_TopN_Count = length(unique(topn_df$Method)),
    ATF3_Method_Rank1_Count = length(unique(rank1_df$Method)),
    ATF3_Method_P_Significant_Count = length(unique(psig_df$Method)),
    ATF3_Supported_Methods = format_method_detail(supported_df[order(supported_df$Rank), , drop = FALSE]),
    ATF3_Rank1_Methods = collapse_unique(rank1_df$Method),
    stringsAsFactors = FALSE
  )
}

method_split <- split(
  atf3_method_evidence,
  paste(atf3_method_evidence$Input_Type, atf3_method_evidence$TF_Analysis_Name, sep = "___")
)

atf3_method_summary <- do.call(
  rbind,
  Filter(Negate(is.null), lapply(method_split, summarize_method_evidence))
)


# 5. ATF3 TF intersection evidence -------------------------------------------

candidate_files <- list.files(
  TF_SUMMARY_ROOT,
  pattern = "top10_tf_candidates[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)

extract_atf3_intersection <- function(file) {
  meta <- extract_tfsummary_meta(file, expected_section = "intersections")
  if (is.null(meta)) {
    return(NULL)
  }

  dat <- safe_read_csv(file)
  if (is.null(dat) || nrow(dat) == 0 || !"TF" %in% colnames(dat)) {
    return(NULL)
  }

  idx <- which(toupper(trimws(as.character(dat$TF))) == toupper(TARGET_GENE))

  if (length(idx) == 0) {
    return(data.frame(
      Input_Type = meta$Input_Type,
      TF_Analysis_Name = meta$TF_Analysis_Name,
      Intersection_Name = meta$Intersection_Name,
      ATF3_Found_In_Top_Candidates = FALSE,
      Consensus_Rank = NA_real_,
      Required_Methods = NA_character_,
      Number_Of_Methods = NA_integer_,
      Source_Method_Count = NA_integer_,
      Mean_Selected_Rank = NA_real_,
      Best_Selected_Rank = NA_real_,
      CheA3_Integrated_TopRank = NA_real_,
      Result_File = file,
      stringsAsFactors = FALSE
    ))
  }

  row <- dat[idx[which.min(as_num(dat$Consensus_Rank[idx]))], , drop = FALSE]
  required_methods <- cell_chr(row, "Required_Methods")
  number_of_methods <- ifelse(
    is.na(required_methods),
    NA_integer_,
    length(strsplit(required_methods, ";", fixed = TRUE)[[1]])
  )
  consensus_rank <- cell_num(row, "Consensus_Rank")

  data.frame(
    Input_Type = meta$Input_Type,
    TF_Analysis_Name = meta$TF_Analysis_Name,
    Intersection_Name = meta$Intersection_Name,
    ATF3_Found_In_Top_Candidates = TRUE,
    Consensus_Rank = consensus_rank,
    Required_Methods = required_methods,
    Number_Of_Methods = number_of_methods,
    Source_Method_Count = cell_num(row, "Source_Method_Count"),
    Mean_Selected_Rank = cell_num(row, "Mean_Selected_Rank"),
    Best_Selected_Rank = cell_num(row, "Best_Selected_Rank"),
    CheA3_Integrated_TopRank = cell_num(row, "CheA3_Integrated_TopRank"),
    Result_File = file,
    stringsAsFactors = FALSE
  )
}

atf3_intersection_evidence <- do.call(
  rbind,
  Filter(Negate(is.null), lapply(candidate_files, extract_atf3_intersection))
)

summarize_intersection_evidence <- function(df) {
  if (nrow(df) == 0) {
    return(NULL)
  }

  found_df <- df[df$ATF3_Found_In_Top_Candidates, , drop = FALSE]
  topn_df <- found_df[
    !is.na(found_df$Consensus_Rank) & found_df$Consensus_Rank <= TF_TOP_RANK_CUTOFF,
    ,
    drop = FALSE
  ]
  rank1_df <- found_df[
    !is.na(found_df$Consensus_Rank) & found_df$Consensus_Rank == 1,
    ,
    drop = FALSE
  ]

  data.frame(
    Input_Type = df$Input_Type[1],
    TF_Analysis_Name = df$TF_Analysis_Name[1],
    Intersections_Checked = length(unique(df$Intersection_Name)),
    ATF3_Intersection_TopN_Count = length(unique(topn_df$Intersection_Name)),
    ATF3_Intersection_Rank1_Count = length(unique(rank1_df$Intersection_Name)),
    ATF3_Best_Consensus_Rank = if (nrow(topn_df) > 0) min(topn_df$Consensus_Rank, na.rm = TRUE) else NA_real_,
    ATF3_Max_Required_Methods = if (nrow(topn_df) > 0) max(topn_df$Number_Of_Methods, na.rm = TRUE) else NA_real_,
    ATF3_TopN_Intersections = format_intersection_detail(topn_df[order(topn_df$Consensus_Rank), , drop = FALSE]),
    ATF3_Rank1_Intersections = collapse_unique(rank1_df$Intersection_Name),
    stringsAsFactors = FALSE
  )
}

intersection_split <- split(
  atf3_intersection_evidence,
  paste(atf3_intersection_evidence$Input_Type, atf3_intersection_evidence$TF_Analysis_Name, sep = "___")
)

atf3_intersection_summary <- do.call(
  rbind,
  Filter(Negate(is.null), lapply(intersection_split, summarize_intersection_evidence))
)


# 6. Integrated decision table ------------------------------------------------

deg_method_summary <- atf3_method_summary[atf3_method_summary$Input_Type == "DEG", , drop = FALSE]
deg_intersection_summary <- atf3_intersection_summary[
  atf3_intersection_summary$Input_Type == "DEG",
  ,
  drop = FALSE
]

integrated <- merge(
  atf3_deg_summary,
  deg_method_summary,
  by.x = "Analysis_Name",
  by.y = "TF_Analysis_Name",
  all.x = TRUE,
  sort = FALSE
)

integrated <- merge(
  integrated,
  deg_intersection_summary,
  by.x = "Analysis_Name",
  by.y = "TF_Analysis_Name",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("_Method", "_Intersection")
)

zero_if_na <- function(x) {
  x[is.na(x)] <- 0
  x
}

integrated$ATF3_Method_Supported_Count <- zero_if_na(integrated$ATF3_Method_Supported_Count)
integrated$ATF3_Method_Rank1_Count <- zero_if_na(integrated$ATF3_Method_Rank1_Count)
integrated$ATF3_Intersection_TopN_Count <- zero_if_na(integrated$ATF3_Intersection_TopN_Count)
integrated$ATF3_Intersection_Rank1_Count <- zero_if_na(integrated$ATF3_Intersection_Rank1_Count)

integrated$Meets_Primary_DEG_Criterion <- integrated$Significant_By_Primary_Cutoff
integrated$Meets_Strict_FDR_DEG_Criterion <- integrated$Significant_By_FDR_Cutoff
integrated$Meets_TF_Multi_Method_Criterion <- integrated$ATF3_Method_Supported_Count >= 2
integrated$Meets_TF_Intersection_Criterion <- integrated$ATF3_Intersection_TopN_Count >= 2
integrated$Meets_User_Primary_Criterion <- integrated$Meets_Primary_DEG_Criterion &
  (integrated$Meets_TF_Multi_Method_Criterion | integrated$Meets_TF_Intersection_Criterion)
integrated$Meets_User_Strict_Criterion <- integrated$Meets_Strict_FDR_DEG_Criterion &
  integrated$Meets_TF_Multi_Method_Criterion

integrated <- integrated[order(
  -integrated$Meets_User_Primary_Criterion,
  -integrated$Meets_User_Strict_Criterion,
  -integrated$Significant_By_FDR_Cutoff,
  -integrated$ATF3_Method_Supported_Count,
  integrated$P.Value
), ]


# 7. Current project result inventory ----------------------------------------

result_datasets <- list.dirs(file.path("results", DATA_TYPE), recursive = FALSE, full.names = TRUE)
result_datasets <- result_datasets[file.info(result_datasets)$isdir]

project_result_inventory <- do.call(
  rbind,
  lapply(result_datasets, function(dataset_dir) {
    files <- list.files(dataset_dir, recursive = TRUE, full.names = TRUE)
    data.frame(
      Dataset = basename(dataset_dir),
      Total_Files = length(files),
      CSV_Tables = sum(grepl("[.]csv$", files, ignore.case = TRUE)),
      PDF_Plots = sum(grepl("[.]pdf$", files, ignore.case = TRUE)),
      PNG_Plots = sum(grepl("[.]png$", files, ignore.case = TRUE)),
      TF_Result_Files = sum(grepl("/TF|/TF_summary", clean_path(files))),
      DEG_Analysis_Count = length(list.files(file.path(dataset_dir, "tables"), pattern = "all_genes[.]csv$", recursive = TRUE)),
      stringsAsFactors = FALSE
    )
  })
)


# 8. Save outputs -------------------------------------------------------------

write.csv(design_summary, file.path(OUTPUT_ROOT, "analysis_design_summary.csv"), row.names = FALSE)
write.csv(atf3_deg_summary, file.path(OUTPUT_ROOT, "atf3_deg_summary.csv"), row.names = FALSE)
write.csv(atf3_method_evidence, file.path(OUTPUT_ROOT, "atf3_tf_method_evidence.csv"), row.names = FALSE)
write.csv(atf3_method_summary, file.path(OUTPUT_ROOT, "atf3_tf_method_summary.csv"), row.names = FALSE)
write.csv(atf3_intersection_evidence, file.path(OUTPUT_ROOT, "atf3_tf_intersection_evidence.csv"), row.names = FALSE)
write.csv(atf3_intersection_summary, file.path(OUTPUT_ROOT, "atf3_tf_intersection_summary.csv"), row.names = FALSE)
write.csv(integrated, file.path(OUTPUT_ROOT, "atf3_integrated_evidence_summary.csv"), row.names = FALSE)
write.csv(project_result_inventory, file.path(OUTPUT_ROOT, "project_result_inventory.csv"), row.names = FALSE)

primary_hits <- integrated[integrated$Meets_User_Primary_Criterion, , drop = FALSE]
strict_hits <- integrated[integrated$Meets_User_Strict_Criterion, , drop = FALSE]

markdown_lines <- c(
  paste0("# ", DATASET_ID, " ", TARGET_GENE, " evidence summary"),
  "",
  "## Decision rule",
  "",
  paste0("- DEG primary cutoff: ", P_VALUE_COLUMN, " < ", P_VALUE_CUTOFF, " and |logFC| >= ", LOGFC_CUTOFF, "."),
  paste0("- DEG strict cutoff: adj.P.Val < ", ADJ_P_VALUE_CUTOFF, " and |logFC| >= ", LOGFC_CUTOFF, "."),
  paste0("- TF method support: ", TARGET_GENE, " appears within rank <= ", TF_TOP_RANK_CUTOFF, " or passes method P/FDR cutoff where available."),
  "- User primary criterion: DEG primary cutoff plus at least two TF methods supporting ATF3, or at least two consensus TF-intersection candidate lists containing ATF3.",
  "- User strict criterion: DEG strict FDR cutoff plus at least two TF methods supporting ATF3.",
  "",
  "## Primary matched analysis designs",
  ""
)

if (nrow(primary_hits) == 0) {
  markdown_lines <- c(markdown_lines, "No analysis design met the primary criterion.", "")
} else {
  for (i in seq_len(nrow(primary_hits))) {
    row <- primary_hits[i, , drop = FALSE]
    markdown_lines <- c(
      markdown_lines,
      paste0(
        "- ", row$Analysis_Name,
        ": logFC=", signif(row$logFC, 4),
        ", P=", signif(row$P.Value, 4),
        ", adj.P=", signif(row$adj.P.Val, 4),
        ", TF_methods=", row$ATF3_Method_Supported_Count,
        ", TF_intersections=", row$ATF3_Intersection_TopN_Count,
        ", supported_methods=", row$ATF3_Supported_Methods
      )
    )
  }
  markdown_lines <- c(markdown_lines, "")
}

markdown_lines <- c(markdown_lines, "## Strict FDR matched analysis designs", "")

if (nrow(strict_hits) == 0) {
  markdown_lines <- c(markdown_lines, "No analysis design met the strict FDR criterion.", "")
} else {
  for (i in seq_len(nrow(strict_hits))) {
    row <- strict_hits[i, , drop = FALSE]
    markdown_lines <- c(
      markdown_lines,
      paste0(
        "- ", row$Analysis_Name,
        ": logFC=", signif(row$logFC, 4),
        ", P=", signif(row$P.Value, 4),
        ", adj.P=", signif(row$adj.P.Val, 4),
        ", TF_methods=", row$ATF3_Method_Supported_Count,
        ", rank1_methods=", row$ATF3_Rank1_Methods
      )
    )
  }
  markdown_lines <- c(markdown_lines, "")
}

markdown_lines <- c(
  markdown_lines,
  "## Output files",
  "",
  "- analysis_design_summary.csv",
  "- atf3_deg_summary.csv",
  "- atf3_tf_method_evidence.csv",
  "- atf3_tf_method_summary.csv",
  "- atf3_tf_intersection_evidence.csv",
  "- atf3_tf_intersection_summary.csv",
  "- atf3_integrated_evidence_summary.csv",
  "- project_result_inventory.csv",
  ""
)

writeLines(markdown_lines, file.path(OUTPUT_ROOT, "ATF3_evidence_summary.md"))

cat("ATF3 evidence summary saved to: ", OUTPUT_ROOT, "\n", sep = "")
cat("Primary matched designs: ", paste(primary_hits$Analysis_Name, collapse = ", "), "\n", sep = "")
cat("Strict FDR matched designs: ", paste(strict_hits$Analysis_Name, collapse = ", "), "\n", sep = "")
