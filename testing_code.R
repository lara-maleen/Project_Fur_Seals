rm(list=ls())

switch(Sys.info()['user'],
       koen = {setwd("/home/koen/Documents/projects/Project_Fur_Seals/")})

source("IBM_fur_seals_Lara_cluster.R")

# undebug(simulation.fun)
# undebug(simulation.fun)
# newdat <- simulation.fun(p=0.1,time=100)



newdat <- simulation.fun(time = 1e4, #t  
               age = 15, 
               patches = 2, #number of Patches (two different sites: high/low density)
               territories = c(50,50), #number of territories per patch
               mutate = 0.05, #mutationfactor
               #die = 0.18, #mortality rate 
               die.fight = 0.35, #propability to die from fight
               loci.col = c(14:53), #in which columns of the pop matrix are the loci?
               p = 0.1,
               u=80,
               i = -0.8, #intercept for infanticide function
               s = 1.8,
               surv=0.9	
)
tail(newdat)
# simulation.fun()
str(newdat)
newdat[,'N']

plot(newdat[,'N'],type="l",ylim=c(0,500))
lines(newdat[,'N1'],type="l",col="red")
lines(newdat[,'N2'],type="l",col="blue")