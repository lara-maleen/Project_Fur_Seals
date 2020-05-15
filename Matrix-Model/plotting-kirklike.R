rm(list=ls())
wdir <- "~/Documents/projects/Project_Fur_Seals/Matrix-Model/out-kirklike-corrected-p01/" #"/data/home/koen/Kirk_Fur_Seals/out1/"
setwd(wdir)



equil_offs <- function(sfun,t2,a2){
  - sfun(t2)*(1 + t2*(a2 - 1))/((sfun(t2) -1 - sfun(t2)*t2)*(a2 -1))
}

sims_all <- read.csv("simruns.csv")
sims_all$series <- as.numeric(factor(with(sims_all,paste(stype,sval,sslope))))
sims_all$cats <- as.numeric(factor(with(sims_all,paste(stype,sval,sslope,a2))))
load("sfuns")
sim_cats <- sims_all[sims_all$replicate==1,]
pdf("out-graph.pdf")
for(i in unique(sims_all$series)){
  sims_series <- sims_all[sims_all$series == i,]
  sfun <- function(t2) sfuns[[as.character(sims_series$stype[i])]](t2,sims_series$sval[i],sims_series$sslope[i])
  par(mfrow=c(2,2))
  plot(seq(0,1,0.01),sfun(seq(0,1,0.01)),type="l",xlim=c(0,1),ylim=c(0,1),xlab="t2",ylab="s")
  
  # par(mfrow=c(1,3))
  
  for(j in unique(sims_series$cats)){
    # j <- 32
    # cat(j,"\t")
    files <- as.character(sims_series$outfile[sims_series$cats==j])
    cols <- rainbow(length(files))
    
    a2 <- (sims_series$a2[sims_series$cats==j])[1]
    t2s <- seq(0,1,0.01)
    p2s <- equil_offs(sfun,t2s,a2)
    p2s[p2s < 0 ] <- 0
    p2s[p2s>1] <- 1
    
    plot(t2s,p2s,xlim=c(0,1),ylim=c(0,1),xlab="Male type",ylab="Female type",type="n",main=paste("a2 = ",a2))
    abline(h=0.5,col="grey")
    abline(v=0.5,col="grey")
    for(k in 1:length(files)){
      cat(k,"\n")
      if(file.exists(paste(wdir,"raw/",files[k],".csv",sep="")))
      tmp <- read.csv(paste(wdir,"raw/",files[k],".csv",sep=""))
      dum <- read.csv(paste(wdir,"raw/",files[k],".dum",sep=""))
      Afreq <- apply(tmp[,c(-1,-2)],1,FUN = function(x) sum(x*dum$Ascore/2)/sum(x))
      Bfreq <- apply(tmp[,c(-1,-2)],1,FUN = function(x) sum(x*dum$Bscore/2)/sum(x))
      
      Ahetfreq <- apply(tmp[,c(-1,-2)],1,FUN = function(x) sum(x*(dum$Ascore==1))/sum(x))
      Bhetfreq <- apply(tmp[,c(-1,-2)],1,FUN = function(x) sum(x*(dum$Bscore==1))/sum(x))
      #tps <- c(1:50,seq(51,length(Afreq),length.out = 50))
      #Afreq <- Afreq#[tps]
      #Bfreq <- Bfreq#[tps]
      lines(Afreq,Bfreq,col=cols[k])
      points(Afreq[nrow(tmp)],Bfreq[nrow(tmp)],cex=0.5,col=cols[k])
      
      # lines(Ahetfreq,Bhetfreq,col=cols[k],lty=3)
      # points(Ahetfreq[nrow(tmp)],Bhetfreq[nrow(tmp)],cex=0.25,col=cols[k])
      
      # points(2*Afreq[nrow(tmp)]*(1-Afreq[nrow(tmp)]),2*Bfreq[nrow(tmp)]*(1-Bfreq[nrow(tmp)]),pch=3,cex=0.25,col=cols[k])
    }
  }
}
dev.off()
