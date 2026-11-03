source(file.path(getwd(),"00_setup.R"))
suppressPackageStartupMessages(library(lavaan))
d <- read_project_data(analysis_data_path)
deep <- ensure_dir(file.path(result_base_dir,"sem"))
z <- function(x) as.numeric(scale(x))
Y <- data.frame(MMSE=z(d$MMSE),MoCA=z(d$MOCA),ADAS13_reversed=z(-d$ADAS13))
Y$two_score_mean <- (Y$MMSE+Y$MoCA)/2
Y$three_score_mean <- (Y$MMSE+Y$MoCA+Y$ADAS13_reversed)/3
# All score definitions use the same participants, so missingness does not change the comparison.
common <- complete.cases(Y); stopifnot(sum(common)>0)
lookup <- c(ChP="L_ChPICV",sTREM2="L_MSD_STREM2CORRECTED",TAU="L_TAU",PTAU="L_PTAU",ABETA="L_S_ABETA")
covariates <- c("S_PTGENDER","S_AGE","PTEDUCAT","APOE401")
groups <- c("CN","MCI","AD")
rows <- tests <- list()
files <- list.files(file.path(repo_root,"config/models"),pattern="^raw_fullmetric__single__",full.names=TRUE)
for(file in files) {
  model <- strsplit(basename(file),"__",fixed=TRUE)[[1]][3]
  v <- strsplit(model,"_to_",fixed=TRUE)[[1]]
  dat <- data.frame(x=z(d[[lookup[[v[1]]]]]),mediator=z(d[[lookup[[v[2]]]]]),
       S_PTGENDER=z(d$S_PTGENDER),S_AGE=z(d$S_AGE),PTEDUCAT=z(d$L_PTEDUCAT),APOE401=z(d$APOE401))
  for(outcome in names(Y)) {
    for(group in c("Overall",groups)) {
      keep <- common & if(group=="Overall") TRUE else d$S_DX_label==group
      dd <- dat[keep,,drop=FALSE]; dd$Cog <- Y[[outcome]][keep]
      xa <- cbind(1,as.matrix(dd[,c("x",covariates)]))
      xb <- cbind(1,as.matrix(dd[,c("mediator","x",covariates)]))
      indirect <- function(i) unname(lm.fit(xa[i,,drop=FALSE],dd$mediator[i])$coefficients[2]*lm.fit(xb[i,,drop=FALSE],dd$Cog[i])$coefficients[2])
      set.seed(20260926+match(group,c("Overall",groups))*1000)
      draws <- replicate(1000,indirect(sample.int(nrow(dd),nrow(dd),replace=TRUE)))
      ci <- quantile(draws,c(.025,.975),na.rm=TRUE)
      rows[[length(rows)+1]] <- data.frame(id=paste(outcome,model,group,sep="__"),outcome=outcome,model=model,group=group,
          n=nrow(dd),estimate=indirect(seq_len(nrow(dd))),ci_low=ci[1],ci_high=ci[2],valid=sum(is.finite(draws)),B=1000)
    }
    dd <- dat[common,,drop=FALSE]; dd$Cog <- Y[[outcome]][common]; dd$.group <- d$S_DX_label[common]
    # Removing the measurement equation leaves the same structural paths with an observed outcome.
    syntax <- paste(readLines(file,warn=FALSE)[-1],collapse="\n")
    fit <- sem(syntax,data=dd,estimator="MLR",missing="fiml",fixed.x=TRUE,std.lv=TRUE,meanstructure=TRUE,group=".group",group.label=groups)
    test <- lavTestWald(fit,constraints="a_CN*b_CN == a_MCI*b_MCI\na_CN*b_CN == a_AD*b_AD")
    tests[[length(tests)+1]] <- data.frame(outcome=outcome,model=model,n=nrow(dd),stat=test$stat,df=test$df,p=test$p.value)
  }
}
write.csv(do.call(rbind,rows),file.path(deep,"observed_bootstrap.csv"),row.names=FALSE)
tests <- do.call(rbind,tests); tests$q_BH_all40 <- p.adjust(tests$p,"BH")
write.csv(tests,file.path(deep,"observed_group_wald.csv"),row.names=FALSE)
