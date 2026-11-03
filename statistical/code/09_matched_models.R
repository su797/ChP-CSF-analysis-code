source(file.path(getwd(), "00_setup.R"))
source(file.path(repo_root, "R", "matched_models.R"), encoding="UTF-8")
project_config <- dget(file.path(repo_root, "config", "matched_models.R"))
analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
rows <- run_matched_key_models()
# Apply the same standardization to the coefficient and both CI endpoints.
ratio <- rows$std_beta / rows$estimate
rows$std_low <- pmin(rows$conf.low*ratio, rows$conf.high*ratio)
rows$std_high <- pmax(rows$conf.low*ratio, rows$conf.high*ratio)
write.csv(rows, file.path(result_summary_dir, "matched_key_model_sensitivity.csv"), row.names=FALSE)
