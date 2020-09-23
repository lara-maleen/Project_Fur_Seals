rm(list=ls())
setwd("~/Documents/projects/Project_Fur_Seals/Matrix-Model/")
source("matrix_model_kirk_version.R")
source("graphing.R")

#library(profvis)
# profvis({
set.seed(20)
run_sim(filename = "test",surv_off = function(n) 0.5*plogis(20*(n-0.5)),Nt=1e4,A.adv=3,wm = 1,wf=0,min_val_m = 0,min_val_f = 0,maxfreq = 1)
# })
plotting_rep(filename = "test")
bla <- read.csv("test.csv")[,-1]
bladum <- read.csv("test.dum")[,-1]
head(bla)
tail(bla)
head(bladum)
matplot(bla[,-1],type="l")
finpop <- round(as.numeric(bla[nrow(bla),-1]),4)*100
newp <- cbind(bladum,finpop)
aggregate(newp$finpop,by=list(newp$a1,newp$a2),sum)
aggregate(newp$finpop,by=list(newp$b1,newp$b2),sum)
# The female want heterozygous offspring, which they are more likely to get by going to the 'bad' beach, b/c it is the only chance to mate with a shitty male.