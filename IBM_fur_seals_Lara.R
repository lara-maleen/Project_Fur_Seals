#Minimal Model IBM - Fur Seal Project
#Hypothesis 1: Marginal-Male Theory

rm(list=ls())

##### START SIMULATION.RUN-FUNCTION #####
simulation.fun <- function(replicates=1, #number of replicates
                           time=20, #number of generations
                           migrate=0.05, #migrationfactor
                           age=2, #age limit for an individual
                           patches=2, #number of Patches (two different sites: high/low density)
                           territories=c(50,50), #number of territories per patch
                           mutate=0.05, #mutationfactor
                           die=0.05, #level.vector to die
                           die.fight=0.25, #propability to die from fight
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
                           c10=0.11670283,
                           i=0, #intercept for infanticide function
                           s=0.5 #slope for infanticide function
){
switch(Sys.info()['user'],
         Lara = {setwd("C:/Users/Lara/Documents/Studium/WHK/WHK Bielefeld Meike/Project_Fur_Seals/")},
        koen = {setwd("/home/koen/Documents/projects/Project_Fur_Seals/")})
  
source('Gene_generator.R')
  
##### FUNCTIONS #####

fitness.fun <- function(a,b,z,N,Np){ #FITNESS-FUNCTION (a,b = coefficients to change function, z = trait value, N= total pop size, Np = pop size of patch)
    y=a+b*plogis(c1+c2*N+c3*z+c4*(0.5*N-Np)+c5*N^2+c6*z^2+c7*(0.5*N-Np)^2+c8*z*N+c9*z*(0.5*N-Np)+c10*N*(0.5*N-Np))
    return(y)
}
  
  
ID.fun <- function(offspring.vector){ #ID-FUNCTION
    ID.offspring <-   ID.scan:(ID.scan+sum(offspring.vector)-1)
    ID.scan <<- ID.scan + sum(offspring.vector)
    return(ID.offspring)
}
  
trait.fun <- function(row,pop.matrix,value.matrix,loci.matrix){ #TRAIT-VALUE-FUNCTION - used for male quality 
    value.matrix <- matrix(NA,nrow=row,ncol=10) #empty matrix for the trait values for each loci
    for(y in 1:row){ #for each individual
      for(z in 1:10){ 
        value.matrix[y,z] <- 1000*(gen_phen_map[z,loci.matrix[y,z],loci.matrix[y,10+z]]) #Intercept=1, Slope=2, slope multiplied with output from genmap
      }
      pop.matrix[y,4] <- abs(sum(value.matrix[y,]))+1
    }
    return(pop.matrix)
}

choice.fun <- function(population.males, population.total, N.male){
  average.success <- mean(population.total[population.total$gender=='male',]$no.offspring) #calculate the overall average male reproductive success of previous year
  decision.matrix <- matrix(NA,nrow=population.males,ncol=1)
  for(y2 in 1:population.males){ #for each male
    if (N.male$no.offspring[y2]>=average.success){ #If reproductive success is equal or greater than total pop average than patch stays the same (from last year)
      decision.matrix[y2] <- N.male$patch.last.year[y2]
    }
    else{ #Otherwise the patch is changed in contrary patch 
      if(N.male$patch.last.year[y2]==1){
        decision.matrix[y2] <- 2
      }
      else{
        decision.matrix[y2] <- 1
      }
    }
    N.male[y2,2] <- decision.matrix[y2]
  }
  return(N.male) #New patch is written into patch column in male matrix
}


competition.fun <- function(N.male, patches, population.males, territories){ #LET MALES COMPETE FOR TERRITORIES, DEPENDING ON THEIR QUALITY TRAIT
  ### 1.) Males choose their territory in this patch
  if(nrow(N.male)>0){ #Are their any males?
    terr.matrix <- c() #for storing the territory number
    for(p2 in 1:patches){ #Going through previous determined patches of males (at first Patch I than Patch II)
      if(nrow(N.male[N.male$patch==p2,])>0){ #Are their any males in the patch
        terr.matrix.patch <- matrix(0,nrow(N.male[N.male$patch==p2,]),ncol=1) #matrix for territories obtained (just for males in these patch)
        for(i2 in 1:length(terr.matrix.patch)){ #go through all males 
          terr.matrix.patch[i2] <- sample(territories[p2], 1) #randomly decide which territory male goes to
          N.male[N.male$patch==p2,]$terr[i2] <- terr.matrix.patch[i2] #write the territory number in matrix of males in this patch 
        }#End individual's loop
      }
    }#End 1.) patch loop
  }#End 1) Step
  
  ### 2) Males compete for their territory - the one with highest quality obtains it
  
  male.matrix <- c()
  for(p3 in 1:patches){ #Go again trough all patches
    male.matrix2 <- c()
    for(t in 1:territories[p3]){ #loop over all territory numbers (1-50)
      matrix.terr <- N.male[which(N.male[,"terr"]==t&N.male[,"patch"]==p3),] #Choose all males in this particular territory (as matrix)
      if(nrow(matrix.terr)>=2){ #If there are at least two in the territory...
        winner <- matrix.terr[matrix.terr$trait==(max(matrix.terr[,"trait"])),] #That's the WINNER in this territory
        matrix.terr <- matrix.terr[which(matrix.terr$ID!=winner$ID),] #remove winner from matrix
        for (i4 in 1:nrow(matrix.terr)){ #For the looser(s) change territory to NA
          matrix.terr$terr[i4] <- 0
        }
        male.matrix2 <- rbind(male.matrix2, winner, matrix.terr) #Safe new info in matrix 
      }
      else{ #What happens when there is just one male (or zero) with this territory number? 
        winner <- N.male[which(N.male[,"terr"]==t&N.male[,"patch"]==p3),] #He "wins" and is added to matrix
        male.matrix2 <- rbind(male.matrix2, winner) 
      }
    }#End territory loop
    male.matrix <- rbind(male.matrix,male.matrix2)
  }#End 2) step
  
  N.male <- male.matrix
  return(N.male)
  
}


statistic.fun <- function(pop.matrix, Npatch){ #PATCH/STATISTIC-FUNCTION
    tmp <- aggregate(pop.matrix$trait,by=list(patch = pop.matrix$patch),mean)
    traits <- tmp$x[match(1:Npatch,tmp$patch)] 
    
    cbind(table(factor(pop.matrix$patch,levels=1:Npatch)),
          table(factor(pop.matrix[pop.matrix$gender=='male',]$patch,levels=1:Npatch)),
          table(factor(pop.matrix[pop.matrix$gender=='female',]$patch,levels=1:Npatch)),
          as.numeric(traits))
}
  
##### MATRICES FOR PLOTS ##### with 100 spaces for mean values (1 replicat, 100 years (time))
  meantrait.matrix <- matrix(NA, nrow=replicates, ncol=time) #empty matrix for the mean trait value of the generations in each replicate
  meanpopulation.matrix <- matrix(NA, nrow=replicates, ncol=time) #empty matrix for the mean populationsize of the generations in each replicate
  meantrait.matrix.males <- matrix(NA, nrow=replicates, ncol=time) #empty matrix for the mean male quality trait value of the generations in each replicate
  meanN.patches.array <- array(NA,dim=c(patches,time,replicates)) #empty array for the populationsize of a patch of the generations in each replicate
  meanM.patches.array <- array(NA,dim=c(patches,time,replicates)) #empty array for the populationsize of males of a patch of the generations in each replicate
  meanF.patches.array <- array(NA,dim=c(patches,time,replicates)) #empty array for the populationsize of females of a patch of the generations in each replicate
  meanquality.patches.array <- array(NA,dim=c(patches,time,replicates)) #empty array for the quality of males of a patch of the generations in each replicate
  
##### REPLICATION LOOP START#####
  
for(r in 1:replicates){
    
    
    ##### INITIALISATION PATCHES #####
    population.total <- c() #empty vector for the population matrix
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
      patch.last.year <- ceiling(runif(patchx.N, min=0, max=2)) #generates randomly ID of last years patch for each individual (patch 1 or 2)
      no.offspring <- c(rep(0,patchx.N)) #no offspring in first generation, will be filled with males success/offspring from last year
      terr <- c(rep(NA, patchx.N)) #here the obtained territory is stored, emptied every t
      
      patchx <- data.frame(ID,patch,gender,trait,survival,ID.mother,ID.father, patch.last.year, no.offspring, terr) #the dataframe is constructed for each patch including all vectors which where defined just before
      population.total <- rbind(population.total,patchx)  #data frame including all individuals of all patches (the dataframe of a patch is included in the population matrix)
    }
    
    population.total$ID <- c(1:nrow(population.total)) #the first generation of the population becomes a new ID
    patchnumbers.vector <- c(1:patches) #vector of patchnumbers
    
    ID.scan <- nrow(population.total)+1
    
    ##### STATISTIC START #####
    population.N <- rep(0,time) #empty vector for the populationsize of each generation 
    population.meantrait <- rep(0,time) #empty vector for the mean traitvalue of each generation
    population.meantrait.males <- rep(0,time) #empty vector for the mean male quality traitvalue of each generation
    
    population.N[1] <- nrow(population.total) #the populationsize for the first generation is written into the vector
    population.meantrait[1] <- mean(population.total$trait) #the mean traitvalue for the first generation is written into the vector
    population.meantrait.males [1] <- mean(population.total$trait[gender=="male"]) #the mean male quality traitvalue for the first generation is written into the vector
    
    ########STATISTIC END  #####
   
    population <- nrow(population.total) #number of individuals
    loci.total <- matrix(NA,nrow=population,ncol=20+1) #empty matrix for the locis (20 numbers) and the ID of the individual (+1 number)
    #loser.matrix <- c() #empty loser matrix 
    
    for(x in 1:population){ #LOOP OVER THE INDIVIDUALS
      loci.total[x,] <- ceiling(runif(21,1e-16,10)) #each individual has 20 random numbers (first 10:row //last 10:column)
      loci.total[x,21] <- x #the last vector-spot is defined as x (the ID of the individual) for the first generation
    }
    
    population.total <- trait.fun(population,population.total,values.population,loci.total) #traitvalue-function: traitvalues for the population are included and overwrite the population matrix
    
    ##### GENERATION LOOP START #####  
    for(t in 1:time){
      #population.total <- rbind(population.total, loser.matrix) #include losers from last year again in pop matrix
      N <- nrow(population.total) #number of individuals in total (all patches included)
      
      if(N>0) { #START IS ANYBODY THERE-LOOP: if there are any individuals and the population is not extinct 
        N.local <- c() #empty vector for local populationsize
        N.male <- subset(population.total,population.total$gender=="male") #number of male individuals in total
        population.males <- nrow(N.male) #number of male individuals
        level.vector <- c() #empty vector
        
        ##### MALE PATCH CHOICE - WHICH PATCH THEY GO #####
        patchbook_males <- c()
        N.male <- choice.fun(population.males, population.total, N.male) #Males decide where to go this year depending on last years success
        patchbook_males <- N.male$patch #overwrite patch from previous year
        population.total[population.total$gender=='male',]$patch <- patchbook_males
        
        ##### MALE CHOICE END #####
        
        ##### MALE COMPETITION - HOW MANY TERRITORIES #####
        population.total$terr <- c(rep(0, nrow(population.total))) #empty territory vector for all indivduals, every t
        terrbook_males <- c()
        N.male <- competition.fun(N.male, patches, population.males, territories) #territories are obtained after competition of males 
        terrbook_males <- N.male$terr
        population.total[population.total$gender=='male',]$terr <- terrbook_males #obtained territories of "winners" are written into pop.matrix
        
        #All males that lost territory competition have certain mortality:
        dying.males <- matrix(NA,nrow(population.total[population.total$gender=="male"&population.total$terr==0,]),ncol=2) #empty matrix for all losers
        dying.males[,2] <- runif(nrow(population.total[population.total$gender=="male"&population.total$terr==0,]),0,1) < die.fight #for each individual is a random number distributed. if the number is below the deathrate a true is written into the vector + ID
        dying.males[,1] <- population.total[population.total$gender=="male"&population.total$terr==0,]$ID #IDS of the males are written here 
        dying.males.ID <- c()
        dying.males.ID <- dying.males[,1][dying.males[,2]==1] #IDs of the males that died are stored
        for(d2 in dying.males.ID){ #go trough the died males and change survival number and the loci matrix 
          population.total[population.total$ID==d2,]$survival <- 0
          loci.total[loci.total[,21]==d2][1] <- -2 
        }
        
        #Update all population info after males died 
        loci.total <- subset(loci.total,loci.total[,1]>(-2 )) #loci matrix: all rows with a -2 in the beginning are deleted
        population.total <-subset(population.total,population.total$survival>0) #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
        N.male <- subset(population.total,population.total$gender=="male") 
        N <- nrow(population.total)
        N.male.patch <- table(factor(N.male$patch,levels = 1:patches)) #number of males in each patch (as a vector)
        
        ##### MALE COMPETITION - FIGHT FOR TERRITORIES II - Remaining males #####
        
        #Let males that lost in previous fight switch to other patch
        #Alternative: N.male[N.male$terr==0,]$patch <- round(runif(1,1,patches)) #how the exclude the existing patch?
        
        decision.matrix <- matrix(NA,nrow(N.male),ncol=1)
        for(i3 in 1:nrow(N.male)){ #for each male
          if (N.male$terr[i3]>0){ #If reproductive success is equal or greater than total pop average than patch stays the same (from last year)
            decision.matrix[i3] <- N.male$patch[i3]
          }
          else{ #Otherwise the patch is changed in contrary patch 
            if(N.male$patch[i3]==1){
              decision.matrix[i3] <- 2
            }
            else{
              decision.matrix[i3] <- 1
            }
          }
          N.male[i3,2] <- decision.matrix[i3]
        }
        
        #Males choose their territory again, fight again (also with existing male) !NOT WORKING YET
        
        #terrbook_males <- c()
        #N.male <- competition.fun(N.male, patches, population.males, territories) 
        #terrbook_males <- N.male$terr
        #population.total[population.total$gender=='male',]$terr <- terrbook_males
        
        #All males that lost territory competition have certain mortality:
        #dying.males <- matrix(NA,nrow(population.total[population.total$gender=="male"&population.total$terr==0,]),ncol=2)
        #dying.males[,2] <- runif(nrow(population.total[population.total$gender=="male"&population.total$terr==0,]),0,1) < die.fight #for each individual is a random number distributed. if the number is below the deathrate the individual is written into a vector
        #dying.males[,1] <- population.total[population.total$gender=="male"&population.total$terr==0,]$ID
        #dying.males.ID <- c()
        #dying.males.ID <- dying.males[,1][dying.males[,2]==1]
        #for(d2 in dying.males.ID){
         # population.total[population.total$ID==d2,]$survival <- 0
          #loci.total[loci.total[,21]==d2][1] <- -2 
        #}
        
        #Update the population matrix
         #population.total <-subset(population.total,population.total$survival>0) #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
       
        #loser.matrix <- subset(population.total,population.total$terr==0&population.total$gender=="male") #store the loser males in a matrix, they will be included in pop matrix next t again
        #population.total <- rbind(subset(population.total,population.total$terr>0),subset(population.total,population.total$gender=="female")) #Exclude losers from population.matrix, as they don't reproduce
        #N.male <- subset(population.total,population.total$gender=="male") 
        #N <- nrow(population.total)
        #N.male.patch <- table(factor(N.male$patch,levels = 1:patches)) #number of males in each patch (as a vector)
        #loci.total <- subset(loci.total,loci.total[,1]>(-2)) #loci matrix: all rows with a -2 in the beginning are deleted
        
        ##### COMPETITION END #####
        
        ### FEMALE PATCH CHOICE ### After male patch and territory decision, because females arrive afterwards in nature
        #patchbook_females <- c() #empty vector for patch females go - WRITE IT IN LATER!
        N.female <- subset(population.total,population.total$gender=="female") #number of female individuals in total
        N.female.patch <- table(factor(N.female$patch,levels = 1:patches)) #number of females in each patch (as a vector)
        #depending on previous density (from last year) and the fitness of existing males?
        
        ### FEMALE CHOICE END ###

        for(pls in 1:patches){ #START IS OFFSPRING POSSIBLE?
          level.vector <- c(level.vector,nlevels(subset(population.total,population.total$patch==pls)$gender)) #create a Vector which shows how many different arguments(levels) are in a vector 
        }
        
        if(max(level.vector)==2){ #if one patch contains both genders then it has a level of 2
          N.0 <- N/500
          
          for(j in 1:patches){ #loop over patches
            N.local <- c(N.local,nrow(subset(population.total,population.total$patch==j))/500) #vector of local population sizes
          }
          
          if(nrow(N.female)>0){ #number of offspring per female
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
          
          current.offspring <- 1 #counter that keeps track of how much offspring have emerged so far during the loop below
          
          if(nrow(N.female)>0){ #START ANY FEMALES?: loop starts if there is at least one female individual
            for(u in 1:nrow(N.female)){ #START LOOP PARTNERFINDING/mother 
              if(offspring.vector[u]>0){ #START GETS THE MOTHER OFFSPRING?
                mother <- N.female$ID[u] #gives the ID of the mother
                ID.mother.offspring <- c(ID.mother.offspring, rep(mother,offspring.vector[u])) #ID of the mother is written into the vector for all her offspring
                
                ###FATHER####
                no_offspring_vector <- c() #empty vector for number of offspring per father 
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
                    
                    #OFFSPRING SURVIVAL
                    
                    
                  } #END LOOP NUMBER CHILDREN
                } #END ANY MALES IN THE PATCH OF THE MOTHER?
              } #END GETS THE MOTHER OFFSPRING?
            } #END LOOP PARTNERFINDING/mother
            
            population.total$no.offspring <- table(factor(ID.father.offspring,levels=population.total$ID)) #writing the number of offspring into the no_offspring columns, stored for one t
            patchbook <- rep(N.female$patch,offspring.vector) #each offspring becomes the patchnumber of the mother
            ID.offspring <- ID.fun(offspring.vector) #the ID of the offspring is calculated by the ID-function and written into the vector for their ID
            trait.offspring <- c(rep(0,length(patchbook))) #the traitvalue of the offspring is set to 0 for the moment
            survival.offspring <- c(rep(age,length(patchbook))) #each offspring gets the survival of the age limit pre defined
            gender.offspring <- genderbook #genders of the offspring are written into the matrix
            patch.offspring <- patchbook #patches of offspring are written into the matrix
            no.offspring.offspring <-  c(rep(0,length(patchbook))) #empty column for the subsequent offspring they will get
            terr.offspring <- c(rep(NA, length(patchbook))) #empty column for subsequent territory
            population.offspring <- data.frame(ID.offspring,patch.offspring,gender.offspring,trait.offspring,survival.offspring,ID.mother.offspring,ID.father.offspring, patch.offspring, no.offspring.offspring, terr.offspring) #a new dataframe is made for the offspring of this generation
            colnames(population.offspring) <- c("ID","patch","gender","trait","survival","ID.mother","ID.father", "patch.last.year", "no.offspring","terr") #column names of the dataframe
            loci.offspring[,21] <- ID.offspring #the ID of the offspring is written into the matrix of the locis of the offspring
            
            population.offspring <- trait.fun(sum(offspring.vector),population.offspring,values.offspring,loci.offspring) #the offspring matrix is overwritten including the traitvalues calculated by the traitvalue-function
            
            #INFATICIDE: Let offspring die with mortality depending on patch density
            infanticide.vector <- c(rep(NA,patches))
            for(p4 in 1:patches){ #for each patch a specific mortality/infanticide rate, depending on density on the patch
              y=i+s*plogis((nrow(population.offspring[population.offspring$patch==p4,])/territories[p4])/25) #mortality is created 
              infanticide.vector[p4] <- y #safed in vector 
            }
            
            dying.offspring <- matrix(NA,nrow(population.offspring),ncol=2)
            dying.offspring[,2] <- runif(nrow(population.offspring),0,1) < infanticide.vector[population.offspring$patch] #for each individual is a random number distributed. if the number is below the deathrate the individual is written into a vector
            dying.offspring[,1] <- population.offspring$ID
            dying.offspring.ID <- c()
            dying.offspring.ID <- dying.offspring[,1][dying.offspring[,2]==1]
            for(d3 in dying.offspring.ID){
              population.offspring[population.offspring$ID==d3,]$survival <- 0
              loci.offspring[loci.offspring[,21]==d3][1] <- -2 
            }
            
            population.offspring <-subset(population.offspring,population.offspring$survival>0) #remove all dead offspring
            loci.offspring <- subset(loci.offspring,loci.offspring[,1]>(-2 )) #loci matrix: all rows with a -2 in the beginning are deleted
            
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
        
        ###Statistic 2##
        population.N[t] <-nrow(population.total) #overwrites the populationsizes for each generation in the empty vector
        population.meantrait[t] <- mean(population.total$trait)#trait.N1.vector[t] <- mean(pop$trait[pop$patch==1]) #overwrites the average trait-value for each generation in the empty vector
        population.meantrait.males[t] <- mean(population.total$trait[gender=="male"])#trait.N1.vector[t] <- mean(pop$trait[pop$patch==1]) #overwrites the average male quality trait-value for each generation in the empty vector
        
        statistic.total[,,t] <- statistic.fun(population.total,patches) #fills the arry with the statistic: N-pop, m-pop, w-pop, mean trait
        ##### End Statistic 2#############
        
      }#END IS ANYBODY THERE? 
      print(t)
    }##### END GENERATION LOOP #####
    meantrait.matrix[r,] <- population.meantrait #writes the mean traitvalues of the replicate into the matrix for all replicates
    meantrait.matrix.males[r,] <- population.meantrait.males #writes the mean male quality traitvalues of the replicate into the matrix for all replicates
    meanpopulation.matrix[r,] <- population.N #writes the mean populationsize of the replicate into the matrix for all replicates
    
    for(pat in 1:patches){
      meanN.patches.array[pat,,r] <- statistic.total[pat,1,] #the populationsize of all generations of patch pat is written into an array for each replicate
      meanM.patches.array[pat,,r] <- statistic.total[pat,2,] #the populationsize of males of all generations of patch pat is written into an array for each replicate
      meanF.patches.array[pat,,r] <- statistic.total[pat,3,] #the populationsize of females of all generations of patch pat is written into an array for each replicate
      meanquality.patches.array[pat,,r] <- statistic.total[pat,4,] #the quality of males of all generations of patch pat is written into an array for each replicate
      }
    
    
}##### END REPLICATION LOOP #####
  meantrait.replicates <- colMeans(meantrait.matrix) #calculates the mean traitvalue of a generation over all replicates
  meantrait.males.replicates <- colMeans(meantrait.matrix.males) #calculates the mean male quality traitvalue of a generation over all replicates
  meanpopulationsize.replicates <- colMeans(meanpopulation.matrix) #calculates the mean populationsize of a generation over all replicates
  meanN.patches.replicates <- rowMeans(meanN.patches.array, dims=2) #calculates the mean populationsize of a generation of a patch over all replicates
  meanM.patches.replicates <- rowMeans(meanM.patches.array, dims=2) #calculates the mean populationsize of male of a generation of a patch over all replicates
  meanF.patches.replicates <- rowMeans(meanF.patches.array, dims=2) #calculates the mean populationsize of female of a generation of a patch over all replicates
  meanquality.patches.replicates <- rowMeans(meanquality.patches.array, dims=2) #calculates the mean quality of males in a generation of a patch over all replicates
  #})
  
##### PLOTS #####
layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4,1))
par(mai=rep(0.8,4.5))

colours <- c("goldenrod1","deepskyblue") #first = patch 1, second = patch 2
  
#First plot: DENSITY - mean population sizes over time/per patch 
i <- c(1,2)
plot(meanpopulationsize.replicates/min(territories),main="Density over time", xlab="Time",ylab="Density",type="l",col="white",ylim =c(0,max(meanN.patches.replicates/min(territories)))) 
  for(ink in 1:patches){
    lines(meanN.patches.replicates[ink,]/territories[ink],type="l",col=colours[ink])
  #}
}
  
#Second plot: MALE QUALITY - average trait value over time/per patch 
plot(meantrait.males.replicates,main="Average male quality over time", xlab="Time",ylab="Mean Quality Trait",col="white",type="l",ylim = c(min(meanquality.patches.replicates),max(meanquality.patches.replicates)))
  for(ink2 in 1:patches){
    lines(meanquality.patches.replicates[ink2,],type="l",col=colours[ink2])
  }
  
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=2,legend=c("Patch I","Patch II"), fill=c("goldenrod1","deepskyblue")) 

}#END SIMULATION.RUN

#Run function 
simulation.fun()
