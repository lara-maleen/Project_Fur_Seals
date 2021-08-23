rm(list=ls())

library(foreach)
library(doMC)
cores <- detectCores()
registerDoMC(cores)

# Code dir
cdir <- "/data/home/koen/Fur_Seals/Code/"
sources <- c("matrix_model_kirk_version.R","logs.R","stat2.R")
# Out dir
odir <- "/data/home/koen/Fur_Seals/out-absdens-2/"
odir.raw <- paste(odir,"raw",sep="")

setwd(cdir)
sapply(sources,source)

if(!dir.exists(odir)){
  dir.create(odir)
  dir.create(odir.raw)
}else{
  stop("Refusing to overwrite existing directory")
}

sfuns <- list(constant = function(t2,sval,sslope) rep(sval,length(t2)),
             linear = function(t2,sval,sslope) sval*t2,
              logistic = function(t2,sval,sslope) sval*plogis(sslope*(t2-0.25)),
             allee = function(t2,sval,sslope) sval*(exp(-sslope*t2) + plogis(sslope*(t2-0.5)))
             )


### Generating the table with all the scenarios that will be run
# run_sim(N0=runif(nrow(dummy)),min_val_m=0.05,min_val_f=0.05,Nt=100,surv_off=sfun2,test=TRUE,tol=-10,A.adv = 1.8,Apenalty = 0.3,wm = 1,wf=0,d=0.5,d2 = 0.5)

# determine focal
focal <- list(a2=1.8,sval=0.05,sslope=2*7.5,stype=c('logistic'),wm=0.5,d=0.5,d2=0.5,mvm=0.05,mvf=0.05,Apenalty=0.3,random_father=FALSE,rel_dens=FALSE)
# potential graphs:
# 1. a2 x Apenalty
# 2. sval x Apenalty
# 3. mvm x mvf
# first (1,2)
var.vars <- list(mvm=seq(0,0.5,0.05),mvf=seq(0,0.5,0.05))#,sval=seq(0,0.09,0.01)) # variable levels for the variables

Nrep <- 15

# option A: full factorial
if(TRUE){
allvar <- lapply(names(focal),FUN = function(x){
   tmp <- c(focal[[x]],var.vars[[x]])
   if(is.numeric(tmp)){
     tmp <- tmp[order(tmp)]
     tmp <- tmp[c(TRUE,abs(diff(tmp))>1e-14)]
     tmp
   }else{
     tmp[!duplicated(tmp)]
   }
  })
names(allvar) <- names(focal)
allvar[['replicate']] <- 1:Nrep
pre_sims <- do.call(expand.grid,allvar)
# remove multi slopes for constant survival functions
sims <- pre_sims[!(pre_sims$stype%in%c("constant","linear") & pre_sims$sslope!=min(pre_sims$sslope)),]
}else{
# option B: focal and varieties on that focal
 
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
pre_sims_rep <- expand.grid(rown=1:nrow(pre_sims),rep=1:Nrep)
sims <- pre_sims[pre_sims_rep$rown,]
sims$replicate <-pre_sims_rep$rep
}
sims$wf <- 1-sims$wm

Nt <- 1e4
# parameter values

sims$seed <- sample(.Machine$integer.max,nrow(sims),replace=FALSE)#1:nrow(sims)
sims$outfile <- paste("out",formatC(1:nrow(sims),width=4,flag="0"),sep="")
setwd(odir)

write.csv(sims,"simruns.csv")
write.csv(as.data.frame(focal),"focal.csv")
save(sfuns,file = "sfuns")

setwd(cdir)
construct_log(odir,c(sources,"main-kirk-2.R"),sims)


#tmp <- statistics(empty=TRUE)
#out.stat <- as.data.frame(matrix(NA,ncol=ncol(tmp)+1,nrow=nrow(sims)))
#colnames(out.stat) <- c(colnames(tmp),'outfile')
#out.stat$outfile <- sims$outfile

setwd(odir.raw)
out.stat <- foreach(i = 1:nrow(sims)) %dopar% {
 set.seed(sims$seed[i])
 with(sims[i,],run_sim(outfile,surv=0,surv_off=function(n) sfuns[[as.character(stype)]](n,sval,sslope),A.adv=a2,wm=wm,wf=wf, min_val_m = mvm,min_val_f = mvf,Nt=Nt,d=d,d2=d2,maxfreq = 1,Apenalty = Apenalty,random_father=random_father,rel_dens=rel_dens))
  # stats <- statistics(sims$outfile[i],surv = sims$surv[i],A.adv=sims$A.adv[i],dens_reg=sims$dens_reg[i],Tp=10,maxfreq = sims$maxfreq[i])
return(0)
}

get_stats(odir)
# setwd(odir)
# out.stat <- do.call(rbind,out.stat)
# write.csv(out.stat,"summary.stat")
# log files (incl. the code names)
# construct_log()
# source(file.path(cdir,"plotting_kirk_fun.R"))
# plot_one_ax(odir)
# 
