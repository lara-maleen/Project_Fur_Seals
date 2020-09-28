rm(list=ls())
wdir <- "data/home/koen/Fur_Seals/out-kirklike-additive-geom-het" #"/data/home/koen/Kirk_Fur_Seals/out1/"
idirs <- c("data/home/koen/Fur_Seals/out-kirklike-additive-geom-het","data/home/koen/Fur_Seals/out-kirklike-additive-geom-het-dom")
setwd(wdir)

simruns <- lapply(idirs, FUN = function(x) read.csv(paste(x,"simruns.csv",sep="/")))
lapply(1:length(simruns), FUN = function(x){x$idir <- i; x})
sims_all <- do.call(rbind,simruns)#read.csv("simruns.csv")

dumdat <- lapply(idirs,FUN=function(x) read.csv(paste(x,"raw/out001.dum",sep="/")))
if(!identical(dumdat[[1]],dumdat[[2]])){warning("dumdats different!\n")}

ress <- matrix(NA,nrow=nrow(sims_all),ncol=nrow(dumdat[[1]]))
for(i in 1:nrow(sims_all)){
  dm <- read.csv(paste(idir[sims_all$idir[i]],"/raw/",sims_all$outfile[i],".csv",sep=""))[,-1]
  ress[i,]  <- as.numeric(dm[nrow(dm),-1])
}
sumdat <- cbind(sims_all,data.frame(
  Nm1 = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$sex=='m')*x*dumdat$p1)), # Nm1
  Nm2 = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$sex=='m')*x*(1-dumdat$p1))), # Nm2
  Nf1 = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$sex=='f')*x*dumdat$p1)), # Nf1
  Nf2 = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$sex=='f')*x*(1-dumdat$p1))), # Nf2
  NmAA = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Ascore==2 & dumdat$sex=='m')*x)), # NmAA
  NmAa = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Ascore==1 & dumdat$sex=='m')*x)),
  Nmaa = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Ascore==0 & dumdat$sex=='m')*x)),
  NfBB = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==2 & dumdat$sex=='f')*x)), # NfBB
  NfBb = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==1 & dumdat$sex=='f')*x)),
  Nfbb = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==0 & dumdat$sex=='f')*x))
))
library(ggplot2)
library(tidyr)
sumdat$dnf <- sumdat$Nf1 - sumdat$Nf2
sumdat$dnm <- sumdat$Nm1 - sumdat$Nm2
sumdat$A <- sumdat$NmAA + 0.5*sumdat$NmAa
sumdat$B <- sumdat$NfBB + 0.5*sumdat$NfBb
sumdat2 <- gather(sumdat,type,val,dnf:dnm)
sumdat3 <- gather(sumdat,type,val,A:B)

cats <- expand.grid(d=unique(sims_all$d),mvm=unique(sims_all$mvm),mvf=unique(sims_all$mvf))
pdf("out-graph-2.pdf")
for(i in nrow(cats)){
d <- cats$d[i]
mvm <- cats$mvm[i]
mvf <- cat$mvf[i]
bla0 <- ggplot(sumdat2[sumdat2$d==d & sumdat2$mvm==mvm & sumdat2$mvf==mvf,],aes(y=val,x=factor(wf),col=type)) + geom_boxplot() + facet_grid(a2~paste(stype,sval,sslope)) +theme_bw()  + ggtitle(paste("d =",i,", mvm =",mvm,", mvf =",mvf))
bla <- ggplot(sumdat2[sumdat2$d==d & sumdat2$mvm==mvm & sumdat2$mvf==mvf,],aes(y=val,x=paste(stype,sval,sslope),col=type)) + geom_boxplot() + facet_grid(a2~wf) +theme_bw()  + ggtitle(paste("d =",i,", mvm =",mvm,", mvf =",mvf))
bla1 <- ggplot(sumdat3[sumdat3$d==d & sumdat3$mvm==mvm & sumdat3$mvf==mvf,],aes(y=val,x=factor(wf),col=type)) + geom_boxplot() + facet_grid(a2~paste(stype,sval,sslope)) +theme_bw()  + ggtitle(paste("d =",i,", mvm =",mvm,", mvf =",mvf))
bla2 <- ggplot(sumdat3[sumdat3$d==d & sumdat3$mvm==mvm & sumdat3$mvf==mvf,],aes(y=val,x=paste(stype,sval,sslope),col=type)) + geom_boxplot() + facet_grid(a2~wf) +theme_bw() + ggtitle(paste("d =",i,", mvm =",mvm,", mvf =",mvf))

# -------------------------- #

print(bla0)
print(bla)
print(bla1)
print(bla2)
}
dev.off()
