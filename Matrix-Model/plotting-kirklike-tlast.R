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
  dumdat <- read.csv(paste(idirs[sims_all$idir[i]],"/raw/",sims_all$outfile[i],".dum",sep=""))
  Nlast <- as.numeric(dm[nrow(dm),-1])
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
sumdat <- read.csv("/home/koen/Documents/projects/Project_Fur_Seals/Matrix-Model/out-kirklike-additive-geom-het/sumdat.csv")

datN <- gather(sumdat,type,val,ff1:fm1)
datA <- gather(sumdat,type,val,AA:aa)
datB <- gather(sumdat,type,val,BB:bb)

# extracting interpretable time series data from raw csv
# ts: time points
# dat: data file with nclass columns and length(ts) rows
# dum: infomation object with nclass rows containing information on the population classes
vs <- function(ts,dat,dum){
  Nm1 <- apply(dat,1,FUN = function(x) sum(x*dum$p1*as.numeric(dum$sex=='m')))
  Nm2 <- apply(dat,1,FUN = function(x) sum(x*(1-dum$p1)*as.numeric(dum$sex=='m')))
  Nf1 <- apply(dat,1,FUN = function(x) sum(x*dum$p1*as.numeric(dum$sex=='f')))
  Nf2 <- apply(dat,1,FUN = function(x) sum(x*(1-dum$p1)*as.numeric(dum$sex=='f')))
  NmAA <- apply(dat,1,FUN = function(x) sum(x*as.numeric(dum$sex=='m' & dum$Ascore ==2)))
  NmAa <- apply(dat,1,FUN = function(x) sum(x*as.numeric(dum$sex=='m' & dum$Ascore ==1)))
  Nmaa <- apply(dat,1,FUN = function(x) sum(x*as.numeric(dum$sex=='m' & dum$Ascore ==0)))
  NfBB <- apply(dat,1,FUN = function(x) sum(x*as.numeric(dum$sex=='f' & dum$Bscore ==2)))
  NfBb <- apply(dat,1,FUN = function(x) sum(x*as.numeric(dum$sex=='f' & dum$Bscore ==1)))
  Nfbb <- apply(dat,1,FUN = function(x) sum(x*as.numeric(dum$sex=='f' & dum$Bscore ==0)))
  
  data.frame(t=ts,fm1 = Nm1/(Nm1+Nm2), ff1 = Nf1/(Nf1+Nf2), AA = NmAA/(Nmaa + NmAa + NmAA),
             Aa = NmAa/(Nmaa + NmAa + NmAA), aa = Nmaa/(Nmaa + NmAa + NmAA),
             BB = NfBB/(Nfbb + NfBb + NfBB), Bb = NfBb/(Nfbb + NfBb + NfBB),
             bb = Nfbb/(Nfbb + NfBb + NfBB))
}

refvals <- function(dat){
  out <- list
  for(i in 1:nrow(dat)){
    dm <- read.csv(paste(idirs[dat$idir[i]],"/raw/",dat$outfile[i],".csv",sep=""))[,-1]
    dumdat <- read.csv(paste(idirs[dat$idir[i]],"/raw/",dat$outfile[i],".dum",sep=""))
    out[[i]] <- vs(dm[,1],dm[,-1],dumdat)
  }
  out
}
# 
# mplot <- function(m1,dat,xvar,facets=''){
#   tmp1 <- do.call(cbind,lapply(names(m1),FUN = function(x) m1[x] == dat[,x]))
#   ind <- apply(tmp1,1,all)
#   bla0 <- ggplot(dat[ind,],aes(y=val,x=eval(parse(text=xvar)),col=type)) +geom_hline(yintercept = 0.5) + ylab("Fraction on island 1") + xlab(xvar) + geom_boxplot() + facet_grid(eval(parse(text=facets))) +theme_bw()  + ggtitle(do.call(paste,lapply(names(m1),FUN = function(x) paste(x,m1[x]))))
#   print(bla0)
# }

plotref <- function(refvals){
  par(mfrow=c(1,3))
  
  maxt <- max(unlist(lapply(refvals,FUN = function(x) max(x$ts))))

  # ff
  plot(0,0,ylim=c(0,1),xlim=c(1,maxt),main="Reference case - Island occupancy",type="n")
  for(i in 1:length(refvals)){
    lines(refvals[[i]]$ts,refvals[[i]]$ff1,col="red")
    lines(refvals[[i]]$ts,refvals[[i]]$mf1,col="blue")
    legend("topright",legend=c('f','m'),lty=1,col=c("red","blue"))
  }
  # AA, Aa, aa
  plot(0,0,ylim=c(0,1),xlim=c(1,maxt),main="Reference case - Gene frequency",type="n")
  for(i in 1:length(refvals)){
    lines(refvals[[i]]$ts,refvals[[i]]$AA,lty=2)
    lines(refvals[[i]]$ts,refvals[[i]]$Aa,lty=1)
    lines(refvals[[i]]$ts,refvals[[i]]$aa,lty=3)
    legend("topright",legend=c('AA','Aa','aa'),lty=c(2,1,3))
  }
  
  # BB, Bb, bb
  plot(0,0,ylim=c(0,1),xlim=c(1,maxt),main="Reference case - Gene frequency",type="n")
  for(i in 1:length(refvals)){
    lines(refvals[[i]]$ts,refvals[[i]]$BB,lty=2)
    lines(refvals[[i]]$ts,refvals[[i]]$Bb,lty=1)
    lines(refvals[[i]]$ts,refvals[[i]]$bb,lty=3)
    legend("topright",legend=c('BB','Bb','bb'),lty=c(2,1,3))
  }
}

caseplot <- function(sumdat,mc,curvar){
  par(mfrow=c(1,3))
  mc2 <- mc[curvar != names(mc)]
  
  tmp1 <- do.call(cbind,lapply(names(mc2),FUN = function(x) mc2[x] == sumdat[,x]))
  ind <- apply(tmp1,1,all)
  
  ddat <- sumdat[ind,]
  refcase <- ddat$curvar == mc$curvar
  
  plot(ddat[,curvar],ddat$ff1,ylim=c(0,1),pch=1+3*as.numeric(refcase),xlab=curvar,col="red")
  points(ddat[,curvar],ddat$fm1,pch=1+3*as.numeric(refcase),col="blue")
  legend("topright",legend=c('f','m'),pch = 1,col=c("red","blue"))
  
  plot(ddat[,curvar],ddat$AA,ylim=c(0,1),pch=1+3*as.numeric(refcase),xlab=curvar,col="red")
  points(ddat[,curvar],ddat$Aa,pch=1+3*as.numeric(refcase),col="blue")
  points(ddat[,curvar],ddat$aa,pch=1+3*as.numeric(refcase),col="black")
  legend("topright",legend=c('AA','Aa','aa'),pch = 1,col=c("red","blue","black"))
  
  plot(ddat[,curvar],ddat$BB,ylim=c(0,1),pch=1+3*as.numeric(refcase),xlab=curvar,col="red")
  points(ddat[,curvar],ddat$Bb,pch=1+3*as.numeric(refcase),col="blue")
  points(ddat[,curvar],ddat$bb,pch=1+3*as.numeric(refcase),col="black")
  legend("topright",legend=c('BB','Bb','bb'),pch = 1,col=c("red","blue","black"))
}

main_plot <- function(sumdat,mc){
  
  # plotting central case
  tmp1 <- do.call(cbind,lapply(names(mc),FUN = function(x) mc[x] == sumdat[,x]))
  ind <- apply(tmp1,1,all)
  if(sum(ind) != max(sumdat$replicate)){
    warning("the number of matching individual cases is not equal to the number of replicates.")
  }
  refvals <- refvals(sumdat[ind,])
  plotref(refvals)
  
  # treatments
  for(curvar in names(mc)){
    if(length(unique(sumdat[,curvar])) > 1){
      # now we have some plotting to do. 
      casplot(sumdat,mc,curvar)
    }
  }
}

m0 <- list(mvm=0,mvf=0,stype='logistic',sval=0.3,wf=0.5,a2=3,d=0.5)
pdf("sumplot.pdf")
main_plot(sumdat,m0)
dev.off()
# 
# m1 <- list(mvm=0,mvf=0,stype='logistic',sval=0.3,wf=0.5,a2=1.2)
# 
# mplot(m0,datA,'as.factor(wf)','a2~.')
# mplot(m0,m2,datB,'as.factor(d)','.~a2')
# 
# 
# mplot(list(d=1,mvm=0,mvf=0))



##### Former graphing code
# cats <- expand.grid(d=unique(sumdat$d),mvm=unique(sumdat$mvm),mvf=unique(sumdat$mvf))
# 
# library(tidyr)
# library(ggplot2)
# sumdat2 <- gather(sumdat,type,val,ff1:fm1)
# sumdat3 <- gather(sumdat,type,val,A:B)
# 
# pdf("out-graph-2.pdf")
# for(i in 1:nrow(cats)){
# d <- cats$d[i]
# mvm <- cats$mvm[i]
# mvf <- cats$mvf[i]
# bla0 <- ggplot(sumdat2[sumdat2$d==d & sumdat2$mvm==mvm & sumdat2$mvf==mvf,],aes(y=val,x=factor(wf),col=type)) +geom_hline(yintercept = 0.5) + ylab("Fraction on island 1") + xlab("wf") + geom_boxplot() + facet_grid(a2~paste(stype,sval,sslope)) +theme_bw()  + ggtitle(paste("d =",d,", mvm =",mvm,", mvf =",mvf))
# bla <- ggplot(sumdat2[sumdat2$d==d & sumdat2$mvm==mvm & sumdat2$mvf==mvf,],aes(y=val,x=paste(stype,sval,sslope),col=type)) + geom_boxplot() + ylab("Fraction on island 1") + facet_grid(a2~wf) +theme_bw()  + ggtitle(paste("d =",d,", mvm =",mvm,", mvf =",mvf))
# bla1 <- ggplot(sumdat3[sumdat3$d==d & sumdat3$mvm==mvm & sumdat3$mvf==mvf,],aes(y=val,x=factor(wf),col=type)) + geom_boxplot() + facet_grid(a2~paste(stype,sval,sslope)) +theme_bw()  + ggtitle(paste("d =",d,", mvm =",mvm,", mvf =",mvf))
# bla2 <- ggplot(sumdat3[sumdat3$d==d & sumdat3$mvm==mvm & sumdat3$mvf==mvf,],aes(y=val,x=paste(stype,sval,sslope),col=type)) + geom_boxplot() + facet_grid(a2~wf) +theme_bw() + ggtitle(paste("d =",d,", mvm =",mvm,", mvf =",mvf))
# # 
# # # -------------------------- #
# # 
# print(bla0)
# print(bla)
# print(bla1)
# print(bla2)
# }
# dev.off()
