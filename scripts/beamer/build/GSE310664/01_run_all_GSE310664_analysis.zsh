#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

DATASET_ID="GSE310664"
DATA_TYPE="ngs"
SCRIPT_ROOT="scripts/${DATASET_ID}"
RESULT_ROOT="results/${DATA_TYPE}/${DATASET_ID}"
TEMP_ROOT="temporary/${DATA_TYPE}/${DATASET_ID}"
CLINICAL_FILE="data/${DATA_TYPE}/${DATASET_ID}/data_prepare/${DATASET_ID}_clinical_edit.csv"
GSEA_QS2_CACHE_ROOT="${TEMP_ROOT}/GSEA_qs2_cache"
OMNIPATHR_LOG_ROOT="${TEMP_ROOT}/omnipathr-log"

R_BIN="${R_BIN:-Rscript}"
VALIDATE_ONLY="${GSE310664_VALIDATE_ONLY:-0}"

# 保留R脚本自己的实时进度条；这里仅负责串联、预检和清理。
export PARALLEL_RUNTIME_BACKEND="${PARALLEL_RUNTIME_BACKEND:-auto}"
export GSEA_REFRESH_QS2_CACHE="${GSEA_REFRESH_QS2_CACHE:-1}"

R_SCRIPTS=(
  "${SCRIPT_ROOT}/00_sample_clustering_heatmap.R"
  "${SCRIPT_ROOT}/01_limma_differential_expression.R"
  "${SCRIPT_ROOT}/02_intersect_significant_genes.R"
  "${SCRIPT_ROOT}/03_volcano_plot.R"
  "${SCRIPT_ROOT}/04_multiple_volcano_plot.R"
  "${SCRIPT_ROOT}/05_top_deg_gene_heatmap.R"
  "${SCRIPT_ROOT}/06_gsea_analysis.R"
  "${SCRIPT_ROOT}/07_gsea_plot.R"
  "${SCRIPT_ROOT}/08_tf_enrichment_analysis.R"
  "${SCRIPT_ROOT}/09_integrate_tf_enrichment_results.R"
)

validate_inputs() {
  "$R_BIN" --vanilla - <<'RSCRIPT'
source("scripts/functions/limma_de_functions.R")

clinical_file <- "data/ngs/GSE310664/data_prepare/GSE310664_clinical_edit.csv"
control_group <- "CTRL"

if (!file.exists(clinical_file)) {
  stop("Missing clinical file: ", clinical_file)
}

clinical <- read.csv(
  clinical_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot("Sample_ID" %in% colnames(clinical))
stopifnot(!any(duplicated(clinical$Sample_ID)))

analysis_designs <- get_analysis_designs(clinical)
errors <- character(0)
warnings <- character(0)
summary_list <- vector("list", nrow(analysis_designs))

for (i in seq_len(nrow(analysis_designs))) {
  design <- analysis_designs[i, , drop = FALSE]
  values <- trimws(as.character(clinical[[design$Column_Name]]))
  values[is.na(values)] <- ""
  used_values <- values[nzchar(values)]
  groups <- unique(used_values)
  experiment_group <- design$Experiment_Group

  if (!experiment_group %in% groups) {
    errors <- c(
      errors,
      paste0(
        design$Column_Name,
        " must contain experiment group '",
        experiment_group,
        "'. Detected groups: ",
        paste(groups, collapse = ", ")
      )
    )
  }

  if (!control_group %in% groups) {
    errors <- c(
      errors,
      paste0(design$Column_Name, " must contain control group '", control_group, "'.")
    )
  }

  non_control_groups <- setdiff(groups, control_group)
  if (length(non_control_groups) != 1) {
    errors <- c(
      errors,
      paste0(
        design$Column_Name,
        " must contain exactly one non-CTRL experiment group. Detected non-CTRL groups: ",
        paste(non_control_groups, collapse = ", ")
      )
    )
  }

  experiment_n <- sum(values == experiment_group)
  control_n <- sum(values == control_group)
  blank_n <- sum(!nzchar(values))

  if (experiment_n < 2 || control_n < 2) {
    warnings <- c(
      warnings,
      paste0(
        design$Analysis_Name,
        " has small group size: ",
        experiment_group,
        " n=",
        experiment_n,
        ", ",
        control_group,
        " n=",
        control_n
      )
    )
  }

  summary_list[[i]] <- data.frame(
    Analysis_Name = design$Analysis_Name,
    Column_Name = design$Column_Name,
    Experiment_Group = experiment_group,
    Experiment_N = experiment_n,
    Control_Group = control_group,
    Control_N = control_n,
    Blank_N = blank_n,
    stringsAsFactors = FALSE
  )
}

cat("\nGSE310664 clinical design check:\n")
print(do.call(rbind, summary_list), row.names = FALSE)

if (length(warnings) > 0) {
  cat("\nWarnings:\n")
  cat(paste0("- ", warnings, collapse = "\n"), "\n", sep = "")
}

if (length(errors) > 0) {
  cat("\nErrors:\n")
  cat(paste0("- ", errors, collapse = "\n"), "\n", sep = "")
  stop("Clinical design check failed.")
}
RSCRIPT
}

clean_previous_outputs() {
  rm -rf \
    "${RESULT_ROOT}/tables" \
    "${RESULT_ROOT}/plots" \
    "${RESULT_ROOT}/intersect" \
    "${RESULT_ROOT}/TF" \
    "${RESULT_ROOT}/TF_summary" \
    "${GSEA_QS2_CACHE_ROOT}" \
    "${OMNIPATHR_LOG_ROOT}" \
    omnipathr-log \
    omipathr-log \
    manuscripts

  mkdir -p "$RESULT_ROOT" "$TEMP_ROOT"
}

run_r_script() {
  local script="$1"
  local index="$2"
  local total="$3"

  if [[ ! -f "$script" ]]; then
    printf "Missing script: %s\n" "$script" >&2
    return 1
  fi

  printf "\n============================================================\n"
  printf "[%02d/%02d] Running: %s\n" "$index" "$total" "$script"
  printf "============================================================\n"
  "$R_BIN" "$script"
}

printf "\nRunning %s full analysis pipeline...\n" "$DATASET_ID"
printf "Project root: %s\n" "$PROJECT_ROOT"
printf "Clinical file: %s\n" "$CLINICAL_FILE"
printf "Parallel backend: %s\n" "$PARALLEL_RUNTIME_BACKEND"
printf "Refresh GSEA qs2 cache: %s\n" "$GSEA_REFRESH_QS2_CACHE"

validate_inputs

if [[ "$VALIDATE_ONLY" == "1" ]]; then
  printf "\nGSE310664_VALIDATE_ONLY=1, stop after clinical design check.\n"
  exit 0
fi

printf "\nCleaning previous %s outputs...\n" "$DATASET_ID"
clean_previous_outputs

start_time="$(date +%s)"
total_scripts="${#R_SCRIPTS[@]}"

for (( i = 1; i <= total_scripts; i++ )); do
  run_r_script "${R_SCRIPTS[$i]}" "$i" "$total_scripts"
done

end_time="$(date +%s)"
elapsed="$(( end_time - start_time ))"
printf "\n%s full analysis pipeline finished in %02d:%02d.\n" \
  "$DATASET_ID" \
  "$(( elapsed / 60 ))" \
  "$(( elapsed % 60 ))"
