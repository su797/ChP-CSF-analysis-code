# Shared paths and model settings. The runner sets CHP_ROOT once for all stages.
repo_root <- Sys.getenv("CHP_ROOT")
if (!nzchar(repo_root)) stop("Run R/run_analysis.R from the repository root.")
script_root <- file.path(repo_root, "statistical", "code")
statistics_root <- script_root
project_root <- repo_root
source(file.path(repo_root, "tools", "utils_io.R"))
source(file.path(repo_root, "tools", "config_utils.R"))
helpers <- list.files(file.path(repo_root, "R"), pattern="\\.R$", full.names=TRUE)
helpers <- helpers[!basename(helpers) %in% c("run_analysis.R", "bootstrap_sem.R", "matched_models.R")]
invisible(lapply(helpers, source))
source(file.path(repo_root, "config", "analysis_config.R"))
runtime_settings <- read_json_config(file.path(repo_root, "config", "analysis_settings.json"))
project_config <- merge_lists_recursive(project_config, runtime_settings$project_config_override %||% list())
project_config$results$write_latest_clean_copy <- FALSE
report_digits <- 4
report_language <- "en"
locale_bundle <- list()
project_report_title <- project_config$project_name
set.seed(project_config$seed)

result_base_dir <- Sys.getenv("CHP_OUTPUT", file.path(repo_root, "statistical", "result"))
result_run_dir <- result_base_dir
result_summary_dir <- ensure_dir(file.path(result_base_dir, "association"))
result_tables_dir <- ensure_dir(file.path(result_base_dir, "tables"))
result_report_dir <- ensure_dir(file.path(result_base_dir, "diagnostics"))
# Matched samples and intermediate datasets stay in a local working directory.
result_data_clean_dir <- ensure_dir(Sys.getenv("CHP_WORK"))
result_figures_dir <- ensure_dir(file.path(result_data_clean_dir, "diagnostic_plots"))
latest_clean_dir <- result_data_clean_dir
result_version <- "analysis"
raw_data_path <- Sys.getenv("CHP_RAW_DATA")
raw_data_all_path <- Sys.getenv("CHP_RAW_ALL")
analysis_data_path <- Sys.getenv("CHP_DATA")
selected_data_path <- file.path(result_data_clean_dir, "ChpSTREM2AD_selected_dataset.csv")
selected_data_latest_path <- selected_data_path
analysis_data_latest_path <- analysis_data_path
normality_table_path <- file.path(result_summary_dir, "normality_screening.csv")
conversion_log_path <- file.path(result_summary_dir, "numeric_conversion_log.csv")
transformation_table_path <- file.path(repo_root, "config", "transformation_plan.csv")
if (Sys.getenv("CHP_PREPARE") == "1") transformation_table_path <- file.path(result_summary_dir, "transformation_plan.csv")
cognition_map_path <- file.path(result_summary_dir, "cognition_component_map.csv")
variable_dictionary_path <- file.path(result_summary_dir, "selected_variable_dictionary.csv")
# Keep run logs local as they may contain input paths.
append_analysis_log <- function(...) invisible(NULL)
