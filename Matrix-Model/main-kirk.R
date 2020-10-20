rm(list=ls())

# Requires:
# readr
# ggplot2
# gridExtra
library(foreach)
library(doMC)
cores <- detectCores()
registerDoMC(cores)

# Code dir
cdir <- "/data/home/koen/Fur_Seals/Code/"
sources <- c("matrix_model_kirk_version.R","logs.R")
# Out dir
odir <- "/data/home/koen/Fur_Seals/out-kirklike-mvm-mvf-scan/"
odir.raw <- paste(odir,"raw",sep="")

setwd(cdir)
sapply(sources,source)

if(!dir.exists(odir)){
  dir.create(odir)
  dir.create(odir.raw)
}else{
  stop("Refusing to overwrite existing directory")
}

# factor 2 codes for the fact that the matrix model looks at total number of males (typically 0 - 0.5), instead of their frequency.
# hence, a correction has to be made, since the old model looked at frequency of type 2 males (0-1) >> not true the model looks at frequency of males!
sfuns <- list(constant = function(t2,sval,sslope) rep(sval,length(t2)),
             linear = function(t2,sval,sslope) sval*t2,
              logistic = function(t2,sval,sslope) sval*plogis(sslope*(t2-0.5)))

# determine focal
focal <- list(a2=3,sval=0.05,sslope=7.5,stype=c('logistic'),wm=0.5,d=0.5,d2=0.5,mvm=0.05,mvf=0.05)
var.vars <- list(stype="constant",mvm=seq(0,1,0.05),mvf=seq(0,1,0.05)) # variable levels for the variables

pre_sims <- do.call(rbind,lapply(names(var.vars),FUN = function(curvar){
  if(!curvar %in% names(focal)){stop("unknown variable supplied to var.vars")}

  varlevels <- var.vars[[curvar]]
  if(focal[[curvar]] %in% varlevels){
    varlevels <- varlevels[! varlevels == focal[[curvar]]]
  }
  if(length(varlevels) > 0){
  varlevels <- varlevels[order(varlevels)]
  curdf <- focal
  curdf[[curvar]] <- varlevels
  as.data.frame(curdf,stringsAsFactors=FALSE)
  }else{
    return(NULL)
  }
}))

pre_sims <- rbind(pre_sims,as.data.frame(focal))

Nrep <- 20
pre_sims_rep <- expand.grid(rown=1:Nrep,rep=1:Nrep)
sims <- pre_sims[pre_sims_rep$rown,]
sims$rep <-pre_sims_rep$rep
sims$wf <- 1-sims$wm

Nt <- 1e5
# parameter values

sims$seed <- sample(.Machine$integer.max,nrow(sims),replace=FALSE)#1:nrow(sims)
sims$outfile <- paste("out",formatC(1:nrow(sims),width=3,flag="0"),sep="")
setwd(odir)

write.csv(sims,"simruns.csv")
write.csv(as.data.frame(focal),"focal.csv")
save(sfuns,file = "sfuns")

setwd(cdir)
construct_log(odir,c(sources,"main-kirk.R"),sims)


#tmp <- statistics(empty=TRUE)
#out.stat <- as.data.frame(matrix(NA,ncol=ncol(tmp)+1,nrow=nrow(sims)))
#colnames(out.stat) <- c(colnames(tmp),'outfile')
#out.stat$outfile <- sims$outfile

setwd(odir.raw)
out.stat <- foreach(i = 1:nrow(sims)) %dopar% {
 set.seed(sims$seed[i])
 with(sims[i,],run_sim(outfile,surv=0,surv_off=function(n) sfuns[[stype]](n,sval,sslope),A.adv=a2,wm=wm,wf=wf, min_val_m = mvm,min_val_f = mvf,Nt=Nt,d=d,d2=d2,maxfreq = 1))
  # stats <- statistics(sims$outfile[i],surv = sims$surv[i],A.adv=sims$A.adv[i],dens_reg=sims$dens_reg[i],Tp=10,maxfreq = sims$maxfreq[i])
return(0)
}
# setwd(odir)
# out.stat <- do.call(rbind,out.stat)
# write.csv(out.stat,"summary.stat")
# log files (incl. the code names)
# construct_log()
source("plotting_kirk_fun.R")
plot_one_ax(odir)
# 
