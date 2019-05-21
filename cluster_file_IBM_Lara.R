
##Minimal Model IBM - Fur Seal Project 
#Hypothesis 1: Marginal-Male Theory

#Parallel parameter value testing for cluster (use file: IBM_fur_seals_Lara_cluster.R)
library(foreach)
library(doMC)
cores <- detectCores()
registerDoMC(cores)

setwd("/data/home/lara")
  
  source('IBM_fur_seals_Lara_cluster.R') #load in file with the IBM model (source for simulation.fun) 
  
  #create new folder with todays date (change 'test' if necessary to number etc), also the working directory is changed to this new folder
  today <- function(prefix = "Test") {
    newdir <- paste(prefix, Sys.Date(), sep = "_")
    dir.create(newdir)
    setwd(newdir)
  }
  
  today() #new folder created and wd set 
  
  # create a list of the parameter values
  df <- expand.grid(die.fight = seq(0,1,0.25), die = seq(0,1,0.24),  replicate = 1:3) #df contains parameter matrix + replicate number
  df$filename <- paste("output",1:nrow(df),".csv",sep="") #add file names (csv files)
  set.seed(NULL)
  df$seed <- round(runif(nrow(df),0,.Machine$integer.max)) #add seed number that is fix for this run
  
  write.csv(df,file='00-param-overview.csv') #write a file with parameter values + file name to have an overview after run
  
blabla <-  foreach(i = 1:nrow(df)) %dopar% {
    set.seed(df$seed[i])
    
    #... run the simulation function with the right parameter values
    tmp <- simulation.fun(time = 100, #t  
                          age = 5, 
                          patches = 2, #number of Patches (two different sites: high/low density)
                          territories = c(50,50), #number of territories per patch
                          mutate = 0.05, #mutationfactor
                          die = df$die[i], #mortality rate 
                          die.fight = df$die.fight[i], #propability to die from fight
                          loci.col = c(12:31), #in which columns of the pop matrix are the loci?
                          i = 0.1, #intercept for infanticide function
                          s = 0.2)
    
    write.csv(tmp, file=df$filename[i])
    return(NULL)
    }
  