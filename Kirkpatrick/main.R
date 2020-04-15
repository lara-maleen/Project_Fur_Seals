rm(list=ls())

# Requires:
# readr
# ggplot2
# gridExtra
library(foreach)
library(doMC)
cores <- detectCores()
registerDoMC(cores)

# Code dir
cdir <- "/data/home/koen/Kirk_Fur_Seals/Code/"
# cdir <- "~/Documents/projects/Project_Fur_Seals/Kirkpatrick/"
sources <- c("simfun.R","logs.R")
# Out dir
odir <- "/data/home/koen/Kirk_Fur_Seals/out1/"
# odir <- "~/Documents/projects/Project_Fur_Seals/Kirkpatrick/out1/"
odir.raw <- paste(odir,"raw",sep="")

setwd(cdir)
sapply(sources,source)

if(!dir.exists(odir)){
  dir.create(odir)
  dir.create(odir.raw)
}else{
  stop("Refusing to overwrite existing directory")
}

sfuns <- list(constant = function(t2,sval,sslope) rep(sval,length(t2)),
              linear = function(t2,sval,sslope) sval*t2,
              logistic = function(t2,sval,sslope) sval*plogis(sslope*(t2-0.5)))

a2s <- c(1.2,1.5,3)
svals <- c(0.1,0.3,0.5)
sslopes <- c(7.5,20)
start <- expand.grid(x1 = seq(0.05,0.85,0.2), x2 = seq(0.05,0.85,0.2),x3 = seq(0.05,0.85,0.2))
start$x4 <- 1 - rowSums(start)
start <- start[!apply(start,1,FUN = function(x) any(x <0 | x > 1)),]
Nrep <- nrow(start)
tmp1 <- expand.grid(a2=a2s,stype= c('constant','linear'),sval = svals,sslope=NA,replicate=1:Nrep,stringsAsFactors = FALSE)
tmp2 <- expand.grid(a2=a2s,stype=c('logistic'),sval=svals,sslope=sslopes,replicate=1:Nrep, stringsAsFactors = FALSE)
sims <- rbind(tmp1,tmp2)
sims <- cbind(sims,start[match(sims$replicate,1:nrow(start)),])

sims$outfile <- paste("out",formatC(1:nrow(sims),width=3,flag="0"),sep="")

Nt <- 1e6

# parameter values
setwd(odir)

write.csv(sims,"simruns.csv")
save(sfuns,file = "sfuns")

setwd(cdir)
construct_log(odir,c(sources,"main.R"),sims)

setwd(odir.raw)

foreach(i = 1:nrow(sims)) %dopar% {
  with(sims[i,],time_iter_offs_surv(outfile,x1,x2,x3,x4,function(t2) sfuns[[stype]](t2,sval,sslope),a2,Nt=Nt,tol=0.01))
}
