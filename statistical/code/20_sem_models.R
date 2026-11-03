source(file.path(getwd(), "00_setup.R"))
suppressPackageStartupMessages(library(lavaan))
d <- read_project_data(analysis_data_path)
stopifnot(!anyDuplicated(d$PTID))
deep <- ensure_dir(file.path(result_base_dir, "sem"))
groups <- c("CN", "MCI", "AD")
lookup <- c(ChP="L_ChPICV", sTREM2="L_MSD_STREM2CORRECTED", TAU="L_TAU", PTAU="L_PTAU", ABETA="L_S_ABETA")
scale_all <- function(x) as.data.frame(lapply(x, function(v) as.numeric(scale(v))))
make_input <- function(kind, model, definition="raw") {
  v <- strsplit(model,"_to_",fixed=TRUE)[[1]]
  if (kind=="single") x <- data.frame(x=d[[lookup[[v[1]]]]], mediator=d[[lookup[[v[2]]]]])
  if (kind=="parallel") x <- data.frame(x=d$L_ChPICV, mediator_1=d$L_MSD_STREM2CORRECTED,
      mediator_2=d[[if(grepl("PTAU",model)) "L_PTAU" else "L_TAU"]], mediator_3=d$L_S_ABETA)
  if (kind=="serial") x <- data.frame(x=d[[lookup[[v[1]]]]], mediator_1=d$L_MSD_STREM2CORRECTED, mediator_2=d$L_ChPICV)
  # Reverse ADAS13 so all three indicators point toward better cognition.
  indicators <- switch(definition, raw=c("MMSE","MOCA","ADAS13"),
                       mixed=c("L_MMSE","L_MOCA","ADAS13"), log=c("L_MMSE","L_MOCA","L_ADAS13"))
  x$indicator_1 <- d[[indicators[1]]]
  x$indicator_2 <- d[[indicators[2]]]
  x$indicator_3 <- -d[[indicators[3]]]
  x$S_PTGENDER <- d$S_PTGENDER; x$S_AGE <- d$S_AGE
  x$PTEDUCAT <- d$L_PTEDUCAT; x$APOE401 <- d$APOE401
  # Standardize over the full cohort before adding diagnostic group labels.
  scale_all(x)
}
estimates <- diagnostics <- tests <- jobs <- list()
for (file in list.files(file.path(repo_root,"config/models"),pattern="\\.lavaan$",full.names=TRUE)) {
  id <- sub("\\.lavaan$","",basename(file)); parts <- strsplit(id,"__",fixed=TRUE)[[1]]
  kind <- parts[2]; model <- parts[3]; gr <- parts[4]
  input <- make_input(kind,model)
  syntax <- paste(readLines(file,warn=FALSE),collapse="\n")
  a <- list(model=syntax,data=input,estimator="MLR",missing="fiml",fixed.x=TRUE,std.lv=TRUE,meanstructure=TRUE)
  if(gr=="Multigroup") {
    input$.group <- d$S_DX_label; a$data <- input; a$group <- ".group"
    a$group.label <- groups; a$group.equal <- "loadings"
  }
  message("  ",id)
  fit <- do.call(sem,a)
  p <- parameterEstimates(fit,standardized=TRUE); pp <- p
  pp$id <- id; pp$variant <- parts[1]; pp$family <- model; pp$sample <- gr
  estimates[[id]] <- pp
  fm <- fitMeasures(fit,c("df","cfi.robust","rmsea.robust","srmr"))
  diagnostics[[id]] <- data.frame(id=id,n=sum(lavInspect(fit,"nobs")),converged=lavInspect(fit,"converged"),t(fm),check.names=FALSE)
  jobs[[id]] <- list(id=id,name=model,kind=kind,definition=parts[1],group=gr,
       group_equal=if(gr=="Multigroup") "loadings" else NULL,data=input,model=syntax,point=p,start=coef(fit))
  if(gr=="Multigroup") {
    effects <- if(kind=="single") "indirect" else c(paste0("indirect_",1:3),"total_indirect")
    for(effect in effects) {
      constraints <- paste0(effect,"_CN == ",effect,"_MCI\n",effect,"_CN == ",effect,"_AD")
      test <- lavTestWald(fit,constraints=constraints)
      tests[[length(tests)+1]] <- data.frame(kind=kind,model=model,effect=effect,stat=test$stat,df=test$df,p=test$p.value)
    }
  }
}
bind_rows <- function(rows) {
  cols <- unique(unlist(lapply(rows,names)))
  do.call(rbind,lapply(rows,function(p){for(n in setdiff(cols,names(p))) p[[n]]<-if(n=="group")1 else NA;p[,cols]}))
}
write.csv(bind_rows(estimates),file.path(deep,"point_estimates.csv"),row.names=FALSE)
write.csv(bind_rows(diagnostics),file.path(deep,"fit_indices.csv"),row.names=FALSE)
tests <- bind_rows(tests)
# One family: eight single-path and eight parallel-effect omnibus tests.
stopifnot(nrow(tests)==16)
tests$q_BH16 <- p.adjust(tests$p,"BH")
write.csv(tests,file.path(deep,"group_omnibus_tests.csv"),row.names=FALSE)

# Compare measurement scales in the reference path and test nested CFA constraints.
reference <- "ChP_to_sTREM2_to_Cognition"
reference_syntax <- jobs[[paste("raw","single",reference,"Overall",sep="__")]]$model
scale_estimates <- invariance <- list()
for(definition in c("raw","mixed","log")) {
  input <- make_input("single",reference,definition)
  fit <- sem(reference_syntax,data=input,estimator="MLR",missing="fiml",fixed.x=TRUE,std.lv=TRUE,meanstructure=TRUE)
  pe <- parameterEstimates(fit,standardized=TRUE); pe$definition <- definition
  scale_estimates[[definition]] <- pe
  cd <- input[paste0("indicator_",1:3)]; cd$.group <- d$S_DX_label
  fits <- list()
  for(level in c("configural","metric","scalar","strict")) {
    eq <- switch(level,configural=NULL,metric="loadings",scalar=c("loadings","intercepts"),strict=c("loadings","intercepts","residuals"))
    fits[[level]] <- cfa("Cog =~ indicator_1 + indicator_2 + indicator_3",data=cd,estimator="MLR",
       missing="fiml",std.lv=TRUE,meanstructure=TRUE,group=".group",group.label=groups,group.equal=eq)
  }
  for(pair in list(c("configural","metric"),c("metric","scalar"),c("scalar","strict"))) {
    lr <- lavTestLRT(fits[[pair[1]]],fits[[pair[2]]])
    invariance[[length(invariance)+1]] <- data.frame(definition=definition,step=paste(pair,collapse="->"),
       chisq_diff=lr[2,"Chisq diff"],df_diff=lr[2,"Df diff"],p=lr[2,"Pr(>Chisq)"])
  }
}
write.csv(bind_rows(scale_estimates),file.path(deep,"measurement_scale_estimates.csv"),row.names=FALSE)
write.csv(bind_rows(invariance),file.path(deep,"invariance_lrt.csv"),row.names=FALSE)
distribution <- correlations <- list()
for(group in c("Overall",groups)) {
  dd <- if(group=="Overall") d else d[d$S_DX_label==group,,drop=FALSE]
  for(variable in c("MMSE","MOCA","ADAS13")) {
    x <- dd[[variable]]; xx <- x[!is.na(x)]
    distribution[[length(distribution)+1]] <- data.frame(group=group,variable=variable,n=length(xx),missing=sum(is.na(x)),mean=mean(xx),sd=sd(xx),
       min=min(xx),q25=quantile(xx,.25),median=median(xx),q75=quantile(xx,.75),max=max(xx),ceiling_30=if(variable!="ADAS13")sum(xx==30) else NA,zero=sum(xx==0))
  }
  cm <- cor(data.frame(MMSE=dd$MMSE,MOCA=dd$MOCA,NegaADAS13=-dd$ADAS13),use="pairwise.complete.obs")
  for(i in 1:2) for(k in (i+1):3) correlations[[length(correlations)+1]] <- data.frame(group=group,var1=rownames(cm)[i],var2=colnames(cm)[k],cor=cm[i,k])
}
write.csv(bind_rows(distribution),file.path(deep,"score_distribution.csv"),row.names=FALSE)
write.csv(bind_rows(correlations),file.path(deep,"score_correlations.csv"),row.names=FALSE)
if(Sys.getenv("CHP_BOOTSTRAP")=="yes") {
  source(file.path(repo_root,"R/bootstrap_sem.R"))
  ensure_dir(file.path(deep,"bootstrap"))
  files <- sort(list.files(file.path(repo_root,"config/models"),full.names=TRUE),method="radix")
  signature <- unname(tools::md5sum(c(analysis_data_path,files)))
  signature_path <- file.path(deep,"bootstrap_signature.txt")
  if(length(list.files(file.path(deep,"bootstrap")))>0) {
    if(!file.exists(signature_path) || !identical(readLines(signature_path),signature))
      stop("Bootstrap inputs differ or have no signature. Use a new --output directory.")
  }
  writeLines(signature,signature_path)
  # Pairing the seeds across models is needed for the serial/single-path ratio.
  status <- lapply(jobs,fit_job)
  write.csv(bind_rows(status),file.path(deep,"bootstrap_status.csv"),row.names=FALSE)
}
