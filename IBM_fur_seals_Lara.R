#Minimal Model IBM - Fur Seal Project

rm(list=ls())

##### START SIMULATION.RUN-FUNCTION #####
simulation.fun <- function(replicates=1, #number of replicates
                   time=100, #number of generations
                   migrate=0.05, #migrationfactor
                   age=2, #age limit for an individual
                   patches=2, #number of Patches (two different sites: high/low density)
                   territories=20,#number of territories per patch
                   mutate=0.05, #mutationfactor
                   die=0.05, #level.vector to die
                   
                   #fecundity
                   a=0.49649467,
                   b=1.47718931,
                   c1=0.72415095,
                   c2=-0.24464625,
                   c3=0.99490196,
                   c4=-1.31337296,
                   c5=-0.06855583,
                   c6 = 0.32833236,
                   c7=-20.88383990,
                   c8=-0.66263785,
                   c9=2.39334027,
                   c10=0.11670283
){
switch(Sys.info()['user'],
       bio23user = {setwd("/home/bio23user/Documents/Projects/Hiwi-ibm/Hiwi-ibm/")},
       Leron = {setwd("C:/Users/Leron/Desktop/IBM_code/")},
       Anwender = {setwd("C:/Users/Anwender/Desktop/")})

source('Gene_generator.R')

# library(profvis)
# profvis({
  

##### FUNCTIONS #####
fitness.fun <- function(a,b,z,N,Np){ #FITNESS-FUNCTION
  y=a+b*plogis(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+c9*z*(0.5*N-Np)+c10*N*(0.5*N-Np))
  return(y)
}


ID.fun <- function(offspring.vector){ #ID-FUNCTION
  ID.offspring <-   ID.scan:(ID.scan+sum(offspring.vector)-1)
  ID.scan <<- ID.scan + sum(offspring.vector)
  return(ID.offspring)
}


trait.fun <- function(row,pop.matrix,value.matrix,loci.matrix){ #TRAIT-VALUE-FUNCTION
  value.matrix <- matrix(NA,nrow=row,ncol=10) #empty matrix for the trait values for each loci
  for(y in 1:row){ #for each individual
    for(z in 1:10){ 
      value.matrix[y,z] <- gen_phen_map[z,loci.matrix[y,z],loci.matrix[y,10+z]]
    }
    pop.matrix[y,4] <- abs(sum(value.matrix[y,]))
  }
  return(pop.matrix)
}


statistic.fun <- function(pop.matrix, Npatch){ #PATCH/STATISTIC-FUNCTION
  tmp <- aggregate(pop.matrix$trait,by=list(patch = pop.matrix$patch),mean)
  traits <- tmp$x[match(1:Npatch,tmp$patch)] 
  
  cbind(table(factor(pop.matrix$patch,levels=1:Npatch)),
        table(factor(pop.matrix[pop.matrix$gender=='male',]$patch,levels=1:Npatch)),
        table(factor(pop.matrix[pop.matrix$gender=='female',]$patch,levels=1:Npatch)),
        as.numeric(traits))
}


##### MATRICES FOR PLOTS #####
meantrait.matrix <- matrix(NA, nrow=replicates, ncol=time) #empty matrix for the mean trait value of the generations in each replicate
meanpopulation.matrix <- matrix(NA, nrow=replicates, ncol=time) #empty matrix for the mean populationsize of the generations in each replicate
meanN.patches.array <- array(NA,dim=c(patches,time,replicates)) #empty array for the mean populationsize of a patch of the generations in each replicate
meanM.patches.array <- array(NA,dim=c(patches,time,replicates)) #empty array for the mean populationsize of males of a patch of the generations in each replicate
meanF.patches.array <- array(NA,dim=c(patches,time,replicates)) #empty array for the mean populationsize of females of a patch of the generations in each replicate


##### REPLICATION LOOP START#####
for(r in 1:replicates){
  
  
  ##### INITIALISATION PATCHES #####
  population.total <- c() #emptc vector for the population matrix
  statistic.total <- array(NA,dim=c(patches,4,time)) #empty array for the statistics
  
  
  for(k in 1:patches){ #LOOP OVER PATCHES
    patchx.N <- abs(round(rnorm(1, mean=250, sd=10))) #Number of individuals in the patch 
    patchx.male <- round(runif(1,patchx.N/4,3*patchx.N/4)) #Number of males in the patch
    
    ID <- c(1:(patchx.N)) #vector ID: gives each individual an ID
    patch <- c(rep(k,patchx.N)) #vector patch: gives each individual their patch Nr.
    gender <- c(rep("male",patchx.male),rep("female",patchx.N-patchx.male)) #vector gender: is filled with males and females
    trait <- c(rep(0.5,patchx.N)) #vector trait: is for all individuals from both patches set as 0.5
    survival <- c(rep(age,patchx.N)) #vector survival: is for all new individuals of both patches the pre defined age limit 
    ID.mother <- c(rep(NA,patchx.N)) #the first generation has no mother and therefore no ID in the column for the mothers ID
    ID.father <- c(rep(NA,patchx.N)) #the first generation has no father and therefore no ID in the column for the fathers ID
    
    patchx <- data.frame(ID,patch,gender,trait,survival,ID.mother,ID.father) #the dataframe is constructed for each patch including all vectors which where defined just before
    population.total <- rbind(population.total,patchx)  #data frame including all individuals of all patches (the dataframe of a patch is included in the population matrix)
  }
  
  population.total$ID <- c(1:nrow(population.total)) #the first generation of the population becomes a new ID
  patchnumbers.vector <- c(1:patches) #vector of patchnumbers
 
  ID.scan <- nrow(population.total)+1
  
  
  ##### STATISTIC START #####
  population.N <- rep(0,time) #empty vector for the populationsize of each generation 
  population.meantrait <- rep(0,time) #empty vector for the mean traitvalue of each generation
  
  population.N[1] <- nrow(population.total) #the populationsize for the first generation is written into the vector
  population.meantrait[1] <- mean(population.total$trait) #the mean traitvalue for the first generation is written into the vector
  ########STATISTIC END  #####
  
  
    population <- nrow(population.total) #number of individuals
    loci.total <- matrix(NA,nrow=population,ncol=20+1) #empty matrix for the locis (20 numbers) and the ID of the individual (+1 number)
    
    for(x in 1:population){ #LOOP OVER THE INDIVIDUALS
      loci.total[x,] <- ceiling(runif(21,1e-16,10)) #each individual has 20 random numbers (first 10:row //last 10:column)
      loci.total[x,21] <- x #the last vector-spot is defined as x (the ID of the individual) for the first generation
    }
    
    population.total <- trait.fun(population,population.total,values.population,loci.total) #traitvalue-function: traitvalues for the population are included and overwrite the population matrix
    
  
    ##### GENERATION LOOP START #####  
    for(t in 1:time){
      N <- nrow(population.total) #number of individuals in total (all patches included)

      
      if(N>0) { #START IS ANYBODY THERE-LOOP: if there are any individuals and the population is not extinct 
        N.local <- c() #empty vector for local populationsize
        N.female <- subset(population.total,population.total$gender=="female") #number of female individuals in total
        N.male <- subset(population.total,population.total$gender=="male") #number of male individuals in total
        level.vector <- c() #empty vector
        
        for(pls in 1:patches){ #START IS OFFSPRING POSSIBLE?
          level.vector <- c(level.vector,nlevels(subset(population.total,population.total$patch==pls)$gender)) #create a Vector which shows how many different arguments(levels) are in a vector 
        }
        
        if(max(level.vector)==2){ #if one patch contains both genders then it has a level of 2
          N.0 <- N/500
          
          for(j in 1:patches){ #loop over patches
            N.local <- c(N.local,nrow(subset(population.total,population.total$patch==j))/500) #vector of local population sizes
          }
          
          if(nrow(N.female)>0){ #number of offspring per femal
            offspring.vector <- 2*rpois(nrow(N.female),fitness.fun(a,b,N.female$trait,N.0,N.local[N.female$patch])) #each female gets a random number of offspring based on the fitness-function
          }
        
          ID.offspring <- c() #empty vector for the ID of the offspring
          patch.offspring <- c() #empty vector for the patch of the offspring
          gender.offspring <- c() #empty vector for the gender of the offspring
          trait.offspring <- c() #empty vector for the trait of the offspring
          survival.offspring <- c() #each offspring gets the survival of the maximum age
          ID.mother.offspring <- c() #empty vector for the mothers ID of the offspring
          ID.father.offspring <- c() #empty vector for the fathers ID of the offspring
          
          loci.offspring <- matrix(NA,nrow=sum(offspring.vector),ncol=21) #empty vector for the locis of the offspring
          
          
          #### START LOOP PARTNERFINDING #####
          patchbook <- c() #empty vector for the patchnumber of the offspring
          genderbook <- c() #empty vector for the gender of the offspring
          
          N.female.patch <- table(factor(N.female$patch,levels = 1:patches)) #number of females in each patch (as a vector)
          N.male.patch <- table(factor(N.male$patch,levels = 1:patches)) #number of males in each patch (as a vector)
          current.offspring <- 1 #counter that keeps track of how much offspring have emerged so far during the loop below
          
          if(nrow(N.female)>0){ #START ANY FEMALES?: loop starts if there is at least one female individual
            for(u in 1:nrow(N.female)){ #START LOOP PARTNERFINDING/mother 
              if(offspring.vector[u]>0){ #START GETS THE MOTHER OFFSPRING?
                mother <- N.female$ID[u] #gives the ID of the mother
                ID.mother.offspring <- c(ID.mother.offspring, rep(mother,offspring.vector[u])) #ID of the mother is written into the vector for all her offspring
                
              
                ###FATHER####
                if(N.male.patch[N.female$patch[u]]>0){ #START ANY MALES IN THE PATCH OF THE MOTHER?: loop starts if there is at least one male individual in the mothers patch
                  father <- sample(N.male$ID[N.male$patch==N.female$patch[u]],1) #sample the ID of one male which patchnumber is the same as the patchnumber of the mother
                  ID.father.offspring <- c(ID.father.offspring,rep(father,offspring.vector[u])) #ID of the father is written into the vector as often as he becomes offspring with the mother
                  
                  #GENETICS:
                  loci.mother <- loci.total[loci.total[,21]==mother,] #vector of locis of the mother
                  loci.father <- loci.total[loci.total[,21]==father,] #vector of locis of the father
                  loci.child <- rep(0,ncol(loci.total)) #empty vector with fixed length for the locis of the offspring
                  
                  
                    for(o in 1:offspring.vector[u]){ #START LOOP NUMBER CHILDREN per female
                      loci.child[1:10] <- loci.mother[(1:10) +sample(c(0,10),10,replace=TRUE)] #the offspring becomes 10 locis sampled from the mother
                      loci.child[11:20] <- loci.father[(1:10) +sample(c(0,10),10,replace=TRUE)] #the offspring becomes 10 locis sampled from the father
                      
                      #MUTATION
                      if(runif(1,0,1) < mutate){ #if a random number is lower than the mutationrate the offspring becomes a random distributed loci
                        loci.child[round(runif(1,1,20))] <- round(runif(1,1,10))
                      }
                      
                      loci.offspring[current.offspring,] <-  loci.child #connects loci of the offspring to the matrix of the other offspring in this generation
                      current.offspring <- current.offspring + 1
            
                      if(runif(1,0,1)>0.5){ #if random number is higher as 0.5, the offspring is female
                        genderbook <- c(genderbook,"female") #the gender is written in the gender vector for the offspring
                      } else{ #otherwise the offspring is male
                        genderbook <- c(genderbook,"male") #the gender is written in the gender vector for the offspring  
                      }
                      
                    } #END LOOP NUMBER CHILDREN
              } #END ANY MALES IN THE PATCH OF THE MOTHER?
            } #END GETS THE MOTHER OFFSPRING?
          } #END LOOP PARTNERFINDING/mother
          
          patchbook <- rep(N.female$patch,offspring.vector) #each offspring becomes the patchnumber of the mother
          ID.offspring <- ID.fun(offspring.vector) #the ID of the offspring is calculated by the ID-function and written into the vector for their ID
          trait.offspring <- c(rep(0,length(patchbook))) #the traitvalue of the offspring is set to 0 for the moment
          survival.offspring <- c(rep(age,length(patchbook))) #each offspring gets the survival of the age limit pre defined
          gender.offspring <- genderbook #genders of the offspring are written into the matrix
          patch.offspring <- patchbook #patches of offspring are written into the matrix
          population.offspring <- data.frame(ID.offspring,patch.offspring,gender.offspring,trait.offspring,survival.offspring,ID.mother.offspring,ID.father.offspring) #a new dataframe is made for the offspring of this generation
          colnames(population.offspring) <- c("ID","patch","gender","trait","survival","ID.mother","ID.father") #column names of the dataframe
          loci.offspring[,21] <- ID.offspring #the ID of the offspring is written into the matrix of the locis of the offspring
          
          population.offspring <- trait.fun(sum(offspring.vector),population.offspring,values.offspring,loci.offspring) #the offspring matrix is overwritten including the traitvalues calculated by the traitvalue-function
          
          population.total <- rbind(population.total,population.offspring) #the offspring population matrix is added to the general population matrix
          rownames(population.total) <- 1:nrow(population.total) #rownames are overwritten
          loci.total <- rbind(loci.total,loci.offspring) #the offspring loci matrix is added to the general loci matrix
          
        }#END ANY FEMALES?
      }#END IS OFFSPRING POSSIBLE?
        
        
      ##### DEATH START #####
      #death by age:
      population.total$survival[1:N] <- population.total$survival[1:N]-1 #every adult loses one survival counter per generation
      loci.total[population.total$survival==0,1] <- -2  #if the survival is 0, it replaces the first loci with -2
       
      #random Death:
      dying.individuals <- runif(nrow(population.total),0,1) < die #for each individual is a random number distributed. if the number is belo the deathrate the individual is written into a vector
      population.total$survival[dying.individuals] <- 0 #the individuals that where written into the vactor below, become a 0 in their survival
      loci.total[dying.individuals,1] <- -2 #the individuals that where written into the vactor below, become a -2 in the first space of their loci row
        
      #erasing dead individuals:
      loci.total <- subset(loci.total,loci.total[,1]>(-2 )) #loci matrix: all rows with a -2 in the beginning are deleted
      population.total <-subset(population.total,population.total$survival>0) #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
      ##### END DEATH #####
        
        
      ##### MIGRATION START #####
      if(nrow(population.total)>0){ #loop starts if there is at least one individual in the population
        migrating.individuals <- runif(nrow(population.total),0,1) < migrate #draws one uniformly distribued number for every individual deciding if and where the individual migrates to
        population.total$patch[migrating.individuals] <- (population.total$patch[migrating.individuals] - 1 + floor(runif(sum(migrating.individuals),1,patches)))%%patches + 1 #migration
          
        rownames(population.total) <- 1:nrow(population.total) #re-indexing the population
      }
      ##### MIGRATION END #####
        
      
      ###Statistic 2##
      population.N[t] <-nrow(population.total) #overwrites the populationsizes for each generation in the empty vector
      population.meantrait[t] <- mean(population.total$trait)#trait.N1.vector[t] <- mean(pop$trait[pop$patch==1]) #overwrites the average trait-value for each generation in the empty vector
        
      statistic.total[,,t] <- statistic.fun(population.total,patches) #fills the arry with the statistic: N-pop, m-pop, w-pop, mean trait
      ##### End Statistic 2#############
        
      
    }#END IS ANYBODY THERE? 
    print(t)
  }##### END GENERATION LOOP #####
  meantrait.matrix[r,] <- population.meantrait #writes the mean traitvalues of the replicate into the matrix for all replicates
  meanpopulation.matrix[r,] <- population.N #writes the mean populationsize of the replicate into the matrix for all replicates
  
  for(pat in 1:patches){
    meanN.patches.array[pat,,r] <- statistic.total[pat,1,] #the populationsize of all generations of patch pat is written into an array for each replicate
    meanM.patches.array[pat,,r] <- statistic.total[pat,2,] #the populationsize of males of all generations of patch pat is written into an array for each replicate
    meanF.patches.array[pat,,r] <- statistic.total[pat,3,] #the populationsize of females of all generations of patch pat is written into an array for each replicate
  }
  
    
}##### END REPLICATION LOOP #####
meantrait.replicates <- colMeans(meantrait.matrix) #calculates the mean traitvalue of a generation over all replicates
meanpopulationsize.replicates <- colMeans(meanpopulation.matrix) #calculates the mean populationsize of a generation over all replicates
meanN.patches.replicates <- rowMeans(meanN.patches.array, dims=2) #calculates the mean populationsize of a generation of a patch over all replicates
meanM.patches.replicates <- rowMeans(meanM.patches.array, dims=2) #calculates the mean populationsize of male of a generation of a patch over all replicates
meanF.patches.replicates <- rowMeans(meanF.patches.array, dims=2) #calculates the mean populationsize of female of a generation of a patch over all replicates
#})


##### PLOTS #####
#PLOT 1
colours <- c("turquoise","violet","orange","blue","red4","seagreen4")
plot(meanpopulationsize.replicates,main="Population over time", xlab="generations",ylab="population",type="l",col="black",ylim =c(0,max(population.N))) #plot traitvalue
for(ink in 1:patches){
  lines(meanN.patches.replicates[ink,],type="l",col=colours[ink])
}

#PLOT 2
plot(meantrait.replicates,main="mean trait value over time", xlab="generations",ylab="mean trait value",type="l",col="turquoise",ylim = c(min(population.meantrait),max(population.meantrait)))

#PLOT3
par(mfrow=c(2,2))
for(paint in 1:patches){
  plot(meanM.patches.replicates[paint,],main="Frequency of the Sexes Patch x", xlab="generations",ylab="frequency",type="l",col="blue")
  lines(meanF.patches.replicates[paint,],type="l",col="red")
}

}#END SIMULATION.RUN
simulation.fun()
