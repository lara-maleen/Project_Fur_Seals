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
odir <- "/data/home/koen/Fur_Seals/out-kirklike-additive-geom-het/"
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

a2s <- c(1.2,3)
svals <- c(0.05,0.5)
sslopes <- c(7.5)
wms <- c(0,0.5,1)
ds <- c(0,0.5,1)
d2s <- c(0,0.5,1)
mvm <- c(0.05)
mvf <- c(0.05)
# start <- expand.grid(x1 = seq(0.05,0.85,0.2), x2 = seq(0.05,0.85,0.2),x3 = seq(0.05,0.85,0.2))
# start$x4 <- 1 - rowSums(start)
# start <- start[!apply(start,1,FUN = function(x) any(x <0 | x > 1)),]
Nrep <- 20
tmp1 <- expand.grid(a2=a2s,wm=wms,stype= c('constant'),sval = svals,sslope=NA,replicate=1:Nrep,d=ds,d2=d2s,mvm=mvm,mvf=mvf,stringsAsFactors = FALSE)
tmp2 <- expand.grid(a2=a2s,wm=wms,stype=c('logistic'),sval=svals,sslope=sslopes,replicate=1:Nrep,d=ds,d2=d2s,mvm=mvm,mvf=mvf, stringsAsFactors = FALSE)
sims <- rbind(tmp1,tmp2)
sims$wf <- 1-sims$wm
Nt <- 1e5
# parameter values

sims$seed <- sample(.Machine$integer.max,nrow(sims),replace=FALSE)#1:nrow(sims)
sims$outfile <- paste("out",formatC(1:nrow(sims),width=3,flag="0"),sep="")
setwd(odir)

write.csv(sims,"simruns.csv")
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


# 
