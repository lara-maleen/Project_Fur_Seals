rm(list=ls())
# Requires:
# readr
# ggplot2
# gridExtra

# Code dir
cdir <- "/home/koen/Documents/projects/Project_Fur_Seals/Matrix-Model/"
sources <- c("matrix_model.R","stat.R","logs.R")
# Out dir
odir <- "/home/koen/Documents/projects/Project_Fur_Seals/Matrix-Model/out1/"
odir.raw <- paste(odir,"raw",sep="")
setwd(cdir)
sapply(sources,source)

if(!dir.exists(odir)){
  dir.create(odir)
  dir.create(odir.raw)
}else{
  # stop("Refusing to overwrite existing directory")
}


# parameter values
sims <- expand.grid(a=1:3,b=2)
sims$outfile <- paste("out",formatC(1:nrow(sims),width=3,flag="0"),sep="")
setwd(odir)
write.csv(sims,"simruns.csv")

setwd(cdir)
construct_log(odir,c(sources,"main.R"),sims)


tmp <- statistics(empty=TRUE)
out.stat <- as.data.frame(matrix(NA,ncol=ncol(tmp)+1,nrow=nrow(sims)))
colnames(out.stat) <- c(colnames(tmp),'outfile')
out.stat$outfile <- sims$outfile

setwd(odir.raw)
for(i in 1:nrow(sims)){
  run_sim(sims$outfile[i],Nt=10)
  out.stat[i,-ncol(out.stat)] <- statistics(sims$outfile[i],surv = 0,A.adv=1.5)
}
setwd(odir)
write.csv(out.stat,"summary.stat")
# log files (incl. the code names)
# construct_log()


# 