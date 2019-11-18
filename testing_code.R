rm(list=ls())

switch(Sys.info()['user'],
       koen = {setwd("/home/koen/Documents/projects/Project_Fur_Seals/")})

source("IBM_fur_seals_Lara_cluster.R")

# undebug(simulation.fun)
# undebug(simulation.fun)
newdat <- simulation.fun(p=0.1,time=100)

# simulation.fun()
str(newdat)
newdat[,'N']

plot(newdat[,'N'],type="l",ylim=c(0,500))
lines(newdat[,'N1'],type="l",col="red")
lines(newdat[,'N2'],type="l",col="blue")
