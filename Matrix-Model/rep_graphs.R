rm(list=ls())

source("~/Documents/projects/Project_Fur_Seals/Matrix-Model/graphing.R")
setwd("~/Documents/projects/Project_Fur_Seals/Matrix-Model/out1")

ref <- read.csv("simruns.csv")

pdf("rep_graphs.pdf")
for(i in 1:nrow(ref)){
  plotting_rep(paste("raw/",ref$outfile[i],sep=""),main = paste("Surv = ",ref$surv[i],", A.adv = ",ref$A.adv[i],", dens_reg = ",ref$dens_reg[i],sep="") )
  cat(i,"\n")
}
dev.off()
