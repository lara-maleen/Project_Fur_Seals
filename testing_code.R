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
               territories = c(20,20), #number of territories per patch
               mutate = 0.01, #mutationfactor
               #die = 0.18, #mortality rate 
               die.fight = 0.15, #propability to die from fight
               p = 0,
               u=80,
               i = -1.4, #intercept for infanticide function
               s = 2.8,
               surv=0.9,
               gene_file1="genes.rds",
               gene_file2="genes2.rds",
               gene_file3="genes2.rds"
)
tail(newdat)
# simulation.fun()
str(newdat)
newdat[,'N']


## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
## END
draw_ci <- function(x,ucl,lcl,col){
  polygon(c(x,rev(x)),c(lcl,rev(ucl)),col=t_col(col,50),density=NA,border=NA)  
}

par(mfrow=c(2,2))
plot(newdat[,'N'],type="l",ylim=c(0,500))
lines(newdat[,'N.females1'],type="l",col="red")
lines(newdat[,'N.females2'],type="l",col="blue")

plot(0,0,type="n",ylim=c(0,50),xlab="Time",ylab="Trait 1",xlim=c(0,nrow(newdat)))
draw_ci(1:nrow(newdat),newdat[,'meantrait1.ucl'],newdat[,'meantrait1.lcl'],"red")
draw_ci(1:nrow(newdat),newdat[,'meantrait2.ucl'],newdat[,'meantrait2.lcl'],"blue")
lines(newdat[,'meantrait1'],col="red")
lines(newdat[,'meantrait2'],col="blue")

plot(0,0,type="n",ylim=c(-0.2,0.2),xlab="Time",ylab="Female Trait",xlim=c(0,nrow(newdat)))
draw_ci(1:nrow(newdat),newdat[,'meantrait.females1.ucl'],newdat[,'meantrait.females1.lcl'],"red")
draw_ci(1:nrow(newdat),newdat[,'meantrait.females2.ucl'],newdat[,'meantrait.females2.lcl'],"blue")
lines(newdat[,'meantrait.females1'],col="red")
lines(newdat[,'meantrait.females2'],col="blue")

plot(0,0,type="n",ylim=c(-0.2,0.2),xlab="Time",ylab="Male Trait",xlim=c(0,nrow(newdat)))
draw_ci(1:nrow(newdat),newdat[,'meantrait.males1.ucl'],newdat[,'meantrait.males1.lcl'],"red")
draw_ci(1:nrow(newdat),newdat[,'meantrait.males2.ucl'],newdat[,'meantrait.males2.lcl'],"blue")
lines(newdat[,'meantrait.males1'],col="red")
lines(newdat[,'meantrait.males2'],col="blue")

# odd <- newdat[,c('N1','N2')][cbind(1:nrow(newdat),rep(c(1,2),length.out = nrow(newdat)))]
# even <- newdat[,c('N1','N2')][cbind(1:nrow(newdat),rep(c(2,1),length.out = nrow(newdat)))]
# plot(newdat[,'N'],type="l",ylim=c(0,500))
# lines(odd,type="l",col="red")
# lines(even,type="l",col="blue")
