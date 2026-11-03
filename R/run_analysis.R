# Rscript R/run_analysis.R --data /path/analysis.csv --stage association
if (.Platform$OS.type=="windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
args <- commandArgs(trailingOnly=TRUE)
arg <- function(key, default="") {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", key)
  args[i+1L]
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."), winslash="/", mustWork=TRUE)
stage <- arg("--stage", "association")
if (!stage %in% c("prepare", "association", "incremental", "sem", "observed", "all")) stop("Unknown stage: ", stage)
input <- arg("--data", Sys.getenv("CHP_DATA"))
work <- arg("--work", file.path(tempdir(), "chp_analysis"))
dir.create(work, recursive=TRUE, showWarnings=FALSE)
work <- normalizePath(work, winslash="/", mustWork=TRUE)
if (stage == "prepare") input <- file.path(work, "ChpSTREM2AD_analysis_dataset.csv")
if (!nzchar(input) || (stage != "prepare" && !file.exists(input))) stop("Supply --data with the harmonized analysis CSV.")
output <- arg("--output", file.path(root, "statistical", "result"))
dir.create(output, recursive=TRUE, showWarnings=FALSE)
Sys.setenv(CHP_ROOT=root, CHP_DATA=normalizePath(input,winslash="/",mustWork=FALSE),
           CHP_WORK=work, CHP_OUTPUT=normalizePath(output,winslash="/"),
           CHP_RAW_DATA=arg("--raw-data"), CHP_RAW_ALL=arg("--raw-all"),
           CHP_PREPARE=if(stage=="prepare") "1" else "0", CHP_BOOTSTRAP=arg("--bootstrap", "no"))
setwd(file.path(root, "statistical", "code"))
run <- function(file) { message("Running ", file); source(file, encoding="UTF-8") }
if (stage == "prepare") {
  if (!file.exists(Sys.getenv("CHP_RAW_DATA")) || !file.exists(Sys.getenv("CHP_RAW_ALL"))) stop("Supply --raw-data and --raw-all.")
  run("01_clean_project_data.R")
  run("02_prepare_analysis_dataset.R")
}
if (stage %in% c("association", "all")) {
  if (!file.exists(Sys.getenv("CHP_RAW_ALL"))) stop("Supply --raw-all for visit, assay, site and protocol metadata.")
  # Descriptive comparisons precede matching; all sensitivity models use the same transformations.
  for (f in c("02b_descriptive_statistics.R", "03_group_comparisons.R", "04_balance_and_matching.R",
              "05_chp_strem2_analysis.R", "06_tau_analysis.R", "06_ptau_analysis.R",
              "06b_advanced_biomarker_models.R", "06c_abeta_truncation_sensitivity.R",
              "06d_robustness_extensions.R", "06e_overlap_weighting_sensitivity.R",
              "06f_diagnosis_biologic_sensitivity.R", "06g_phase_site_protocol_sensitivity.R",
              "06h_timewindow_lv_sensitivity.R", "02c_sample_flow_missingness_alignment.R",
              "08_plot_estimates.R", "09_matched_models.R")) run(f)
}
if (stage %in% c("incremental", "all")) run("10_incremental_models.R")
if (stage %in% c("sem", "all")) run("20_sem_models.R")
if (stage %in% c("observed", "all")) run("21_observed_scores.R")
capture.output(sessionInfo(), file=file.path(Sys.getenv("CHP_OUTPUT"), "sessionInfo.txt"))
message("Finished: ", stage)
