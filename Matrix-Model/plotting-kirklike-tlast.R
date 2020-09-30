rm(list=ls())
wdir <- "/data/home/koen/Fur_Seals/out-kirklike-additive-geom-het" #"/data/home/koen/Kirk_Fur_Seals/out1/"
idirs <- c("/data/home/koen/Fur_Seals/out-kirklike-additive-geom-het","/data/home/koen/Fur_Seals/out-kirklike-additive-geom-het-dom")
setwd(wdir)

simruns <- lapply(idirs, FUN = function(x) read.csv(paste(x,"simruns.csv",sep="/")))
simruns <- lapply(1:length(simruns), FUN = function(x){tmp <- simruns[[x]]; tmp$idir <- x; tmp})
sims_all <- do.call(rbind,simruns)#read.csv("simruns.csv")

dumdat <- lapply(idirs,FUN=function(x) read.csv(paste(x,"raw/out001.dum",sep="/")))
if(!identical(dumdat[[1]],dumdat[[2]])){warning("dumdats different!\n")}
dumdat <- dumdat[[1]]

ress <- matrix(NA,nrow=nrow(sims_all),ncol=nrow(dumdat))
for(i in 1:nrow(sims_all)){
  dm <- read.csv(paste(idirs[sims_all$idir[i]],"/raw/",sims_all$outfile[i],".csv",sep=""))[,-1]
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
sumdat$ff1 <- sumdat$Nf1/(sumdat$Nf1 + sumdat$Nf2)
sumdat$fm1 <- sumdat$Nm1/(sumdat$Nm1 + sumdat$Nm2)
sumdat$A <- sumdat$NmAA + 0.5*sumdat$NmAa
sumdat$B <- sumdat$NfBB + 0.5*sumdat$NfBb
sumdat$AA <- sumdat$NmAA/(sumdat$NmAA + sumdat$NmAa +sumdat$Nmaa)
sumdat$Aa <- sumdat$NmAa/(sumdat$NmAA + sumdat$NmAa +sumdat$Nmaa)
sumdat$aa <- sumdat$Nmaa/(sumdat$NmAA + sumdat$NmAa +sumdat$Nmaa)
sumdat$BB <- sumdat$NmBB/(sumdat$NmBB + sumdat$NmBb +sumdat$Nmbb)
sumdat$Bb <- sumdat$NmBa/(sumdat$NmBB + sumdat$NmBb +sumdat$Nmbb)
sumdat$bb <- sumdat$Nmbb/(sumdat$NmBB + sumdat$NmBb +sumdat$Nmbb)

write.csv(sumdat,"sumdat.csv")
# 
# sumdat2 <- gather(sumdat,type,val,ff1:fm1)
# sumdat3 <- gather(sumdat,type,val,AA:aa)
# sumdat4 <- gather(sumdat,type,val,BB:bb)
# 
# # central case, to which to compare things
# cent <- list(mvm=0,mvf=0,d=0.5,stype=linear,sval=0.3,wf=0.5)
# 
# #cats <- expand.grid(d=unique(sims_all$d),mvm=unique(sims_all$mvm),mvf=unique(sims_all$mvf))
# 
# pdf("out-graph-2.pdf")
# for(i in nrow(cats)){
# d <- cats$d[i]
# mvm <- cats$mvm[i]
# mvf <- cats$mvf[i]
# 
# bla0 <- ggplot(sumdat2[sumdat2$d==d & sumdat2$mvm==mvm & sumdat2$mvf==mvf,],aes(y=val,x=factor(wf),col=type)) +geom_hline(yintercept = 0.5) + ylab("Fraction on island 1") + xlab("wf") + geom_boxplot() + facet_grid(a2~paste(stype,sval,sslope)) +theme_bw()  + ggtitle(paste("d =",i,", mvm =",mvm,", mvf =",mvf))
# bla <- ggplot(sumdat2[sumdat2$d==d & sumdat2$mvm==mvm & sumdat2$mvf==mvf,],aes(y=val,x=paste(stype,sval,sslope),col=type)) + geom_boxplot() + facet_grid(a2~wf) +theme_bw()  + ggtitle(paste("d =",i,", mvm =",mvm,", mvf =",mvf))
# bla1 <- ggplot(sumdat3[sumdat3$d==d & sumdat3$mvm==mvm & sumdat3$mvf==mvf,],aes(y=val,x=factor(wf),col=type)) + geom_boxplot() + facet_grid(a2~paste(stype,sval,sslope)) +theme_bw()  + ggtitle(paste("d =",i,", mvm =",mvm,", mvf =",mvf))
# bla2 <- ggplot(sumdat3[sumdat3$d==d & sumdat3$mvm==mvm & sumdat3$mvf==mvf,],aes(y=val,x=paste(stype,sval,sslope),col=type)) + geom_boxplot() + facet_grid(a2~wf) +theme_bw() + ggtitle(paste("d =",i,", mvm =",mvm,", mvf =",mvf))
# 
# # -------------------------- #
# 
# print(bla0)
# print(bla)
# print(bla1)
# print(bla2)
# }
# 
# cats2 <- expand.grid(wf=unique(sims_all$wf),mvm=unique(sims_all$mvm),mvf=unique(sims_all$mvf))
# ggplot(sumdat2p[],aes(y=val,x=d))
# dev.off()
