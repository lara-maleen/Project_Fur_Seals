rm(list=ls())
wdir <- "~/Documents/projects/Project_Fur_Seals/Matrix-Model/out-kirklike-additive/" #"/data/home/koen/Kirk_Fur_Seals/out1/"
setwd(wdir)

equil_offs <- function(sfun,t2,a2){
  - sfun(t2)*(1 + t2*(a2 - 1))/((sfun(t2) -1 - sfun(t2)*t2)*(a2 -1))
}
# head(sims_all)

sims_all <- read.csv("simruns.csv")
sims_all$series <- as.numeric(factor(with(sims_all,paste(stype,sval,sslope))))
sims_all$cats <- as.numeric(factor(with(sims_all,paste(stype,sval,sslope,a2,wf,wm))))
load("sfuns")

sim_cats <- sims_all[sims_all$replicate==1,]

dumdat <- read.csv("raw/out001.dum")
ress <- matrix(NA,nrow=nrow(sims_all),ncol=nrow(dumdat))
for(i in 1:nrow(sims_all)){
  dm <- read.csv(paste("raw/",sims_all$outfile[i],".csv",sep=""))[,-1]
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
  NmBB = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==2 & dumdat$sex=='f')*x)), # NfBB
  NmBb = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==1 & dumdat$sex=='f')*x)),
  Nmbb = apply(ress,1, FUN = function(x) sum(as.numeric(dumdat$Bscore==0 & dumdat$sex=='f')*x))
))
library(ggplot2)
library(tidyr)
sumdat$dnf <- sumdat$Nf1 - sumdat$Nf2
sumdat$dnm <- sumdat$Nm1 - sumdat$Nm2
sumdat2 <- gather(sumdat,type,val,dnf:dnm)
ggplot(sumdat2,aes(y=val,x=factor(wf),col=type)) + geom_boxplot() + facet_grid(a2~paste(stype,sval,sslope)) +theme_bw()
bla <- ggplot(sumdat2,aes(y=val,x=paste(stype,sval,sslope),col=type)) + geom_boxplot() + facet_grid(a2~wf) +theme_bw()

# -------------------------- #
pdf("out-graph-2.pdf")
print(bla)
dev.off()
