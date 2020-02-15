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
sources <- c("matrix_model.R","stat.R","logs.R")
# Out dir
odir <- "/data/home/koen/Fur_Seals/out2/"
odir.raw <- paste(odir,"raw",sep="")

setwd(cdir)
sapply(sources,source)

if(!dir.exists(odir)){
  dir.create(odir)
  dir.create(odir.raw)
}else{
  stop("Refusing to overwrite existing directory")
}


# parameter values
sims <- expand.grid(surv=c(0,0.5,0.8,0.9),A.adv=c(1,1.1,1.3,1.5),dens_reg=c(0,0.05,0.1,0.5),rep=1:10,min_val_m = c(0,0.1,0.2,0.3),min_val_f = c(0,0.1,0.2,0.3))
sims$seed <- 1:nrow(sims)
sims$outfile <- paste("out",formatC(1:nrow(sims),width=3,flag="0"),sep="")
setwd(odir)

write.csv(sims,"simruns.csv")

setwd(cdir)
construct_log(odir,c(sources,"main.R"),sims)


#tmp <- statistics(empty=TRUE)
#out.stat <- as.data.frame(matrix(NA,ncol=ncol(tmp)+1,nrow=nrow(sims)))
#colnames(out.stat) <- c(colnames(tmp),'outfile')
#out.stat$outfile <- sims$outfile

setwd(odir.raw)
out.stat <- foreach(i = 1:nrow(sims)) %dopar% {
 set.seed(sims$seed[i])
 run_sim(sims$outfile[i],surv=sims$surv[i],A.adv=sims$A.adv[i],dens_reg = sims$dens_reg[i],min_val_m = sims$min_val_m[i],min_val_f = sims$min_val_f[i],Nt=1e4)
  stats <- statistics(sims$outfile[i],surv = sims$surv[i],A.adv=sims$A.adv[i],dens_reg=sims$dens_reg[i],Tp=10)
return(stats)
}
setwd(odir)
out.stat <- do.call(rbind,out.stat)
write.csv(out.stat,"summary.stat")
# log files (incl. the code names)
# construct_log()


# 
