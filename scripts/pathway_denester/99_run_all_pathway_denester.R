# PathwayDenester all-project runner.
#
# 00 export MSigDB SYMBOL GMT -> 01 run official Python PathwayDenester ->
# 02 plot summaries.


# 0. Config -------------------------------------------------------------------

SCRIPTS_TO_RUN <- c(
  "scripts/pathway_denester/00_export_msigdb_symbol_gmt.R",
  "scripts/pathway_denester/01_run_pathway_denester.py",
  "scripts/pathway_denester/02_plot_pathway_denester_results.R"
)

PYTHON_ARGS <- c(
  "--clone-if-missing",
  "--refresh"
)

options(width = 200)


# 1. Run ----------------------------------------------------------------------

cat("\nRunning all-project PathwayDenester workflow...\n")
cat("Scripts: ", length(SCRIPTS_TO_RUN), "\n\n", sep = "")

rscript_bin <- file.path(R.home("bin"), "Rscript")
python_bin <- Sys.getenv("PATHWAY_DENESTER_PYTHON", unset = "")
venv_python <- file.path("temporary", "pathway_denester_venv", "bin", "python")
if (python_bin == "" && file.exists(venv_python)) {
  python_bin <- venv_python
}
if (python_bin == "") {
  python_bin <- Sys.which("python3")
}
if (python_bin == "") {
  python_bin <- Sys.which("python")
}
if (python_bin == "") {
  stop("No python executable was found in PATH.")
}

for (script_file in SCRIPTS_TO_RUN) {
  if (!file.exists(script_file)) {
    stop("Missing script: ", script_file)
  }

  cat("\n============================================================\n")
  cat("Running: ", script_file, "\n", sep = "")
  cat("============================================================\n")

  if (grepl("[.]py$", script_file)) {
    exit_status <- system2(python_bin, args = c(script_file, PYTHON_ARGS))
  } else {
    exit_status <- system2(rscript_bin, args = script_file)
  }

  if (!identical(exit_status, 0L)) {
    stop("Script failed: ", script_file, " (exit status ", exit_status, ")")
  }
}

cat("\nAll-project PathwayDenester workflow finished.\n")
