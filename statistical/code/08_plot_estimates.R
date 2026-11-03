source(file.path(getwd(),"00_setup.R"))
d <- read_project_data(analysis_data_path)
trans <- read_project_data(transformation_table_path)
covars <- project_config$variables$covariates
factors <- project_config$variables$categorical_covariates
stars <- function(p) if(is.na(p)) "" else if(p<.001) "***" else if(p<.01) "**" else if(p<.05) "*" else ""
biomarkers <- c("MSD_STREM2CORRECTED","TAU","PTAU","S_ABETA")
labels <- c("sTREM2","Tau","P-Tau","Aβ")
rows <- list()
for(i in seq_along(biomarkers)) for(g in c("Overall","CN","MCI","AD")) {
  dd <- if(g=="Overall") d else d[d$S_DX_label==g,,drop=FALSE]
  fit <- fit_linear_model_base(dd,"ChPICV",biomarkers[i],covars,factors,trans)
  fit <- fit[fit$term=="exposure",]
  rows[[length(rows)+1]] <- data.frame(panel=LETTERS[i],panel_title=paste(labels[i],"vs ChP/ICV"),stage=g,
       x=biomarkers[i],y="ChPICV",x_label=labels[i],y_label="ChP/ICV",n=fit$n,beta=fit$std_beta,p_value=fit$p.value)
}
write.csv(do.call(rbind,rows),file.path(result_summary_dir,"combined_regression_stage_panel_stats.csv"),row.names=FALSE)
vars <- c("ChPICV",biomarkers,"MMSE","MOCA")
labels <- c("ChP/ICV","sTREM2","Tau","P-Tau","CSF Aβ42","MMSE","MoCA")
rows <- list()
for(i in seq_along(vars)) for(j in seq_along(vars)) if(i!=j) {
  fit <- run_partial_correlation(d,vars[j],vars[i],covars,factors,trans)
  rows[[length(rows)+1]] <- data.frame(row_var=vars[i],col_var=vars[j],row_label=labels[i],col_label=labels[j],
       estimate=fit$estimate[1],p_value=fit$p_value[1],stars=stars(fit$p_value[1]))
}
write.csv(do.call(rbind,rows),file.path(result_summary_dir,"association_heatmap_partial_correlations.csv"),row.names=FALSE)
rows <- list()
for(pair in project_config$matching$pairs) {
  pair_name <- paste(pair,collapse="_vs_")
  unmatched <- d[d$S_DX_label %in% pair,,drop=FALSE]
  matched <- read_project_data(file.path(result_data_clean_dir,paste0("matched_",pair_name,".csv")))
  for(sample in c("Unmatched","Matched")) {
    dd <- if(sample=="Matched") matched else unmatched
    for(variable in c("Age","Female sex")) {
      x <- if(variable=="Age") dd$S_AGE else as.numeric(dd$S_PTGENDER==2)
      levels <- unique(as.character(dd$S_DX_label))
      a <- x[dd$S_DX_label==levels[1]]; b <- x[dd$S_DX_label==levels[2]]
      pooled_p <- (mean(a)+mean(b))/2
      denominator <- if(variable=="Age") sqrt((var(a)+var(b))/2) else sqrt(pooled_p*(1-pooled_p))
      rows[[length(rows)+1]] <- data.frame(pair=gsub("AD","Dementia",paste(pair,collapse=" vs ")),sample=sample,variable=variable,smd=(mean(a)-mean(b))/denominator)
    }
  }
}
write.csv(do.call(rbind,rows),file.path(result_summary_dir,"matching_smd_summary.csv"),row.names=FALSE)
