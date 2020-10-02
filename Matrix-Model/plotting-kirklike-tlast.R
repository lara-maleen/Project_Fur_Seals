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
summ <- data.frame(Nm1=rep(NA,nrow(ress)),Nm2=NA,Nf1=NA,Nf2=NA,NmAA=NA,NmAa=NA,Nmaa=NA,NfBB=NA,NfBb=NA,Nfbb=NA)
for(i in 1:nrow(sims_all)){
  dm <- read.csv(paste(idirs[sims_all$idir[i]],"/raw/",sims_all$outfile[i],".csv",sep=""))[,-1]
  dumdat <- read.csv(paste(idirs[sims_all$idir[i]],"/raw/",sims_all$outfile[i],".csv",sep=""))
  Nlast <- as.numeric(dumdat[nrow(dm),-1])
  summ$Nm1[i] <- sum(as.numeric(dumdat$sex=='m')*Nlast*dumdat$p1)
  summ$Nm2[i] <- sum(as.numeric(dumdat$sex=='m')*Nlast*(1-dumdat$p1))
  summ$Nf1[i] <- sum(as.numeric(dumdat$sex=='f')*Nlast*dumdat$p1)
  summ$Nf2[i] <- sum(as.numeric(dumdat$sex=='f')*Nlast*(1-dumdat$p1))
  summ$NmAA[i] <-  sum(as.numeric(dumdat$sex=='m' & dumdat$Ascore==2)*Nlast)
  summ$NmAa[i] <-  sum(as.numeric(dumdat$sex=='m' & dumdat$Ascore==1)*Nlast)
  summ$Nmaa[i] <-  sum(as.numeric(dumdat$sex=='m' & dumdat$Ascore==0)*Nlast)
  summ$NfBB[i] <-  sum(as.numeric(dumdat$sex=='f' & dumdat$Bscore==2)*Nlast)
  summ$NfBb[i] <-  sum(as.numeric(dumdat$sex=='f' & dumdat$Bscore==1)*Nlast)
  summ$Nfbb[i] <-  sum(as.numeric(dumdat$sex=='f' & dumdat$Bscore==0)*Nlast)
  
}
sumdat <- cbind(sims_all,summ)#,data.frame(
#   Nm1 = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$sex=='m')*x*dumdat$p1)), # Nm1
#   Nm2 = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$sex=='m')*x*(1-dumdat$p1))), # Nm2
#   Nf1 = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$sex=='f')*x*dumdat$p1)), # Nf1
#   Nf2 = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$sex=='f')*x*(1-dumdat$p1))), # Nf2
#   NmAA = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Ascore==2 & dumdat$sex=='m')*x)), # NmAA
#   NmAa = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Ascore==1 & dumdat$sex=='m')*x)),
#   Nmaa = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Ascore==0 & dumdat$sex=='m')*x)),
#   NfBB = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==2 & dumdat$sex=='f')*x)), # NfBB
#   NfBb = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==1 & dumdat$sex=='f')*x)),
#   Nfbb = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==0 & dumdat$sex=='f')*x))
# ))
library(ggplot2)
library(tidyr)
sumdat$ff1 <- sumdat$Nf1/(sumdat$Nf1 + sumdat$Nf2)
sumdat$fm1 <- sumdat$Nm1/(sumdat$Nm1 + sumdat$Nm2)
sumdat$A <- sumdat$NmAA + 0.5*sumdat$NmAa
sumdat$B <- sumdat$NfBB + 0.5*sumdat$NfBb
sumdat$AA <- sumdat$NmAA/(sumdat$NmAA + sumdat$NmAa +sumdat$Nmaa)
sumdat$Aa <- sumdat$NmAa/(sumdat$NmAA + sumdat$NmAa +sumdat$Nmaa)
sumdat$aa <- sumdat$Nmaa/(sumdat$NmAA + sumdat$NmAa +sumdat$Nmaa)
sumdat$BB <- sumdat$NfBB/(sumdat$NfBB + sumdat$NfBb +sumdat$Nfbb)
sumdat$Bb <- sumdat$NfBb/(sumdat$NfBB + sumdat$NfBb +sumdat$Nfbb)
sumdat$bb <- sumdat$Nfbb/(sumdat$NfBB + sumdat$NfBb +sumdat$Nfbb)

write.csv(sumdat,"sumdat.csv")
# sumdat <- read.csv("/home/koen/Documents/projects/Project_Fur_Seals/Matrix-Model/out-kirklike-additive-geom-het/sumdat.csv")
#  
# datN <- gather(sumdat,type,val,ff1:fm1)
# datA <- gather(sumdat,type,val,AA:aa)
# datB <- gather(sumdat,type,val,BB:bb)
# 
# # central case, to which to compare things
# cent <- list(mvm=0,mvf=0,d=0.5,stype=linear,sval=0.3,wf=0.5)
# 
#cats <- expand.grid(d=unique(sims_all$d),mvm=unique(sims_all$mvm),mvf=unique(sims_all$mvf))
# 
# pdf("out-graph-2.pdf")
# for(i in nrow(cats)){
# d <- 1
# mvm <- 0
# mvf <- 0
# # 
# mplot <- function(m1,dat,xvar,facets=''){
#   tmp1 <- do.call(cbind,lapply(names(m1),FUN = function(x) m1[x] == dat[,x]))
#   ind <- apply(tmp1,1,all)
#   bla0 <- ggplot(dat[ind,],aes(y=val,x=eval(parse(text=xvar)),col=type)) +geom_hline(yintercept = 0.5) + ylab("Fraction on island 1") + xlab(xvar) + geom_boxplot() + facet_grid(eval(parse(text=facets))) +theme_bw()  + ggtitle(do.call(paste,lapply(names(m1),FUN = function(x) paste(x,m1[x]))))
#   print(bla0)  
# }
# 
# m0 <- list(mvm=0,mvf=0,stype='logistic',sval=0.3,wf=0.5,a2=3)
# m1 <- list(mvm=0,mvf=0,stype='logistic',sval=0.3,wf=0.5,a2=1.2)
# 
# mplot(m0,m1,datA,'as.factor(wf)','a2~.')
# mplot(m0,m2,datB,'as.factor(d)','.~a2')
# 
# tmp1[,c(1,8)]
# 
# sumdat2[ind,]
# mplot(list(d=1,mvm=0,mvf=0))
# 
# bla0 <- ggplot(sumdat2[sumdat2$d==d & sumdat2$mvm==mvm & sumdat2$mvf==mvf,],aes(y=val,x=factor(wf),col=type)) +geom_hline(yintercept = 0.5) + ylab("Fraction on island 1") + xlab("wf") + geom_boxplot() + facet_grid(a2~paste(stype,sval,sslope)) +theme_bw()  + ggtitle(paste("d =",d,", mvm =",mvm,", mvf =",mvf))
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
