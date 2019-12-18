rm(list=ls())

switch(Sys.info()['user'],
       koen = {setwd("/home/koen/Documents/projects/Project_Fur_Seals/")})

source("IBM_fur_seals_Lara_cluster.R")

# undebug(simulation.fun)
# undebug(simulation.fun)
# newdat <- simulation.fun(p=0.1,time=100)



newdat <- simulation.fun(time = 1e3, #t  
               age = 15, 
               patches = 2, #number of Patches (two different sites: high/low density)
               territories = c(20,20), #number of territories per patch
               mutate = 0.01, #mutationfactor
               #die = 0.18, #mortality rate 
               die.fight = 0.15, #propability to die from fight
               p = 0,
               u=80,
               i = -1.4, #intercept for infanticide function
               s = 2.8,
               surv=0.9,
               gene_file1="genes_het_ad.rds",
               gene_file2="genes2.rds",
               gene_file3="genes2.rds"
)
tail(newdat)
# simulation.fun()
str(newdat)
newdat[,'N']

par(mfrow=c(2,2))
plot(newdat[,'N'],type="l",ylim=c(0,500))
lines(newdat[,'N.females1'],type="l",col="red")
lines(newdat[,'N.females2'],type="l",col="blue")

plot(newdat[,'meantrait1'],type="l",col="red",ylim=c(0,50))
lines(newdat[,'meantrait2'],col="blue")

plot(newdat[,'meantrait.females1'],type="l",col="red",ylim=c(-0.2,0.2))
lines(newdat[,'meantrait.females2'],col="blue")

plot(newdat[,'meantrait.males1'],type="l",col="red",ylim=c(-0.2,0.2))
lines(newdat[,'meantrait.males2'],col="blue")

# odd <- newdat[,c('N1','N2')][cbind(1:nrow(newdat),rep(c(1,2),length.out = nrow(newdat)))]
# even <- newdat[,c('N1','N2')][cbind(1:nrow(newdat),rep(c(2,1),length.out = nrow(newdat)))]
# plot(newdat[,'N'],type="l",ylim=c(0,500))
# lines(odd,type="l",col="red")
# lines(even,type="l",col="blue")
