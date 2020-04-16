rm(list=ls())
wdir <- "/data/home/koen/Kirk_Fur_Seals/out1/"
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
  par(mfrow=c(1,1))
  plot(seq(0,1,0.01),sfun(seq(0,1,0.01)),type="l",xlim=c(0,1),ylim=c(0,1),xlab="t2",ylab="s")
  par(mfrow=c(1,3))
  
  for(j in unique(sims_series$cats))
    files <- sims_series$outfile[sims_series$cats==j]
    cols <- rainbow(length(files))
    
    a2 <- (sims_series$a2[sims_series$cats==j])[1]
    t2s <- seq(0,1,0.01)
    p2s <- equil_offs(sfun,t2s,a2)
    p2s[p2s < 0 ] <- 0
    p2s[p2s>1] <- 1
    plot(p2s,t2s,xlim=c(0,1),ylim=c(0,1),xlab="p2",ylab="t2",type="l",main=paste("a2 = ",a2))
    
    for(k in 1:length(files)){
      tmp <- read.csv(paste(wdir,"raw/",files[k]))
      lines(tmp$x2+tmp$x4,tmp$x3+tmp$x4,col=cols[k])
      points(tmp$x2[nrow(tmp)]+tmp$x4[nrow(tmp)],temp$x3[nrow(tmp)]+tmp$x4[nrow(tmp)],cex=2,col=cols[k])
    }
  }
dev.off()