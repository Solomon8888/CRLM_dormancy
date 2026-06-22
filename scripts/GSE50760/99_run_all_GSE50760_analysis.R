# GSE50760 analysis pipeline runner.
#
# Default order:
# 00 sample design -> 01 ATF3 expression -> 02 paired DE ->
# 03 within-tissue ATF3 high/low DE/correlation/GSEA -> 04 DE-GSEA ->
# 05 volcano -> 06 multiple volcano -> 07 GSEA dotplots.


# 0. Config -------------------------------------------------------------------

SCRIPTS_TO_RUN <- c(
  "scripts/GSE50760/00_prepare_sample_design.R",
  "scripts/GSE50760/01_atf3_expression_overview.R",
  "scripts/GSE50760/02_limma_differential_expression.R",
  "scripts/GSE50760/03_atf3_correlation_gsea.R",
  "scripts/GSE50760/04_gsea_analysis.R",
  "scripts/GSE50760/05_volcano_plot.R",
  "scripts/GSE50760/06_multiple_volcano_plot.R",
  "scripts/GSE50760/07_gsea_dotplot.R"
)

options(width = 200)


# 1. Run ----------------------------------------------------------------------

cat("\nRunning GSE50760 analysis pipeline...\n")
cat("Scripts: ", length(SCRIPTS_TO_RUN), "\n\n", sep = "")

rscript_bin <- file.path(R.home("bin"), "Rscript")

for (script_file in SCRIPTS_TO_RUN) {
  if (!file.exists(script_file)) {
    stop("Missing script: ", script_file)
  }

  cat("\n============================================================\n")
  cat("Running: ", script_file, "\n", sep = "")
  cat("============================================================\n")

  exit_status <- system2(rscript_bin, args = script_file)

  if (!identical(exit_status, 0L)) {
    stop("Script failed: ", script_file, " (exit status ", exit_status, ")")
  }
}

cat("\nGSE50760 analysis pipeline finished.\n")
