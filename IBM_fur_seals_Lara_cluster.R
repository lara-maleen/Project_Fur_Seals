#Minimal Model IBM - Fur Seal Project
#Hypothesis 1: Marginal-Male Theory
#Code for running in cluster (without replicates and storing information in matrix every t)
rm(list=ls())

setwd("/home/koen/Documents/projects/Project_Fur_Seals/")

sample.vec <- function(x, ...) x[sample.int(length(x), ...)] 

##### START SIMULATION.RUN-FUNCTION #####
simulation.fun <- function(
                           time=1e4, #number of generations
                           age=15, #age limit for an individual; life span for A. gazella. 15-25 years --> literature!?
                           patches=2, #number of Patches (two different sites: high/low density)
                           territories=c(50,50), #number of territories per patch
                           mutate=0.05, #mutation factor
                           die.fight=0.35, #propability to die from fight/competition
                           loci.col=c(14:53), #loci column numbers of the pop matrix 
                           p= 0.25, #parameter for philopatry function (female patch choice) -> the higher p is, the more intense the philopatric (side-fidelity) influence
                           u = 80, #assumed normal average density (for each patch), used for female patch choice function
                           i=-0.8, #intercept for infanticide function
                           s=1.8, #slope for infanticide function
                           surv=0.9, #survival for total population 
                           gene_file1="genes.rds",
                           gene_file2="genes2.rds"
 ){

#gen_phen_map <- readRDS('/data/home/lara/genes.rds') #load the gene array (10 loci, 10 alleles) #gene map used in cluster
#gen_phen_map2 <- readRDS('/data/home/lara/genes2.rds') #load the gene array (10 loci, 10 alleles) #gene map used in cluster
gen_phen_map <- readRDS(gene_file1) #load the gene array (10 loci, 10 alleles), used for male trait values
gen_phen_map2 <- readRDS(gene_file2) #second gene map for female trait value (10 loci, 10 alleles). Phenotype of -0.2 and +0.2 initially

##### FUNCTIONS #####


ID.fun <- function(offspring.vector){ #ID-FUNCTION: for each individual a new ID is created
    ID.offspring <-   ID.scan:(ID.scan+sum(offspring.vector)-1)
    ID.scan <<- ID.scan + sum(offspring.vector)
    return(ID.offspring)
}
  

trait.fun <- function(population.total, loci.matrix, gen_phen_map, gen_phen_map2){ #TRAIT-VALUE-FUNCTION - used for male quality + female philopatry trait 

  #Male Trait Value
  
  # for the male trait value
  lc1 <- as.numeric(t(loci.matrix[,1:10]))
  lc2 <- as.numeric(t(loci.matrix[,11:20]))
  zs <- rep(1:10,nrow(population.total))
  phen <- gen_phen_map[cbind(lc1,lc2,zs)]
  population.total$trait <- colSums(matrix(phen,nrow=10))

  #Female Trait Value:
  lc1 <- as.numeric(t(loci.matrix[,21:30]))
  lc2 <- as.numeric(t(loci.matrix[,31:40]))
  zs <- rep(1:10,nrow(population.total))
  phen <- gen_phen_map2[cbind(lc1,lc2,zs)]
  population.total$female.trait <- colSums(matrix(phen,nrow=10))
  
  return(population.total) 
}


choice.fun <- function(N.male, patches){ #MALE PATCH CHOICE: decide where each adult male goes to this t
  
  # Males that had offspring the year before, stick to the same patch, those that didn't switch patch
  N.male$patch <- (N.male$patch - 1 + as.numeric(N.male$nr.offspring == 0)*floor(runif(nrow(N.male),1,patches)))%%patches + 1

  return(N.male) #New patch is written into patch column in male matrix
}


choice.fun.females <- function(N.female,p,u,N.last1,N.last2, patches){ #FEMALE PATCH CHOICE: determines the patch for females, depending on last years N (from last years patch) as well as density-preference trait & philopatry
  
  N.last <- c(N.last1,N.last2) #get the population size from the previous year per patch
  
  #Add philopatric decision (parameter set at the beginning)
  p.patch <- p > runif(nrow(N.female),0,1) #decide wether female is philopatric or not (stays at birth patch = TRUE), if not then the density-dependent choice takes place
  
  #on the positions where p.patch is TRUE, the patch number is the birth patch:
  N.female$patch[p.patch] <- N.female$patch.born[p.patch]
  
  # for the other patches, it is determined whether or not the patch changes from the current patch
  patch.u  <- plogis(N.female$female.trait[!p.patch]*(N.last[N.female$patch[!p.patch]] - u)) > runif(sum(!p.patch),0,1)

  if(any(is.na(patch.u))){
    print(summary(patch.u))
  }
  if(any(patch.u)){
    N.female$patch[!p.patch][patch.u] <- (N.female$patch[!p.patch][patch.u] - 1 + floor(runif(sum(patch.u),1,patches)))%%patches + 1
  }
  return(N.female)
}


competition.fun <- function(N.male, patches, population.males, territories){ #LET MALES COMPETE FOR TERRITORIES, DEPENDING ON THEIR QUALITY TRAIT
  
  ### 1.) Males choose their territory in this patch
  for(p in 1:patches){ #Going through previous determined patches of males (at first Patch I than Patch II)
    if(nrow(N.male[N.male$patch==p&N.male$terr==0,])>0){ #Are their any males in the patch (with no territory yet)
      ID.terr.males <- matrix(NA, nrow=nrow(N.male[N.male$patch==p&N.male$terr==0,]), ncol=2) #new matrix for storing IDs
      ID.terr.males[,1] <- N.male[N.male$patch==p & N.male$terr==0,]$ID #get IDs of males that have no territory yet
      ID.terr.males[,2] <- sample(territories[p],nrow(ID.terr.males), replace=TRUE)
      N.male$terr[match(ID.terr.males[,1],N.male$ID)] <- ID.terr.males[,2] #write the territory number in matrix of males in this patch 
    }
  }#End 1.) patch loop
  
  ### 2) Males compete for their territory - the one with highest quality trait obtains it

  max_per_terr <- aggregate(N.male$trait,by=list(patch = N.male$patch,terr = N.male$terr),max)
  max_per_terr$dum <- paste(max_per_terr$patch,max_per_terr$terr,sep="-")
  dum <- paste(N.male$patch,N.male$terr,sep="-")
  matched_max <- max_per_terr$x[match(dum,max_per_terr$dum)]
  N.male$terr[N.male$trait < matched_max] <- 0 
  N.male$terr[N.male$terr!=0 & duplicated(N.male$terr)] <- 0 # what if there are multiple males with the same, maximum, trait value
  
  return(N.male)
  
}



mortality <- function(N, surv){ #Calculate density-dependent mortality rate. Dependent on total population size in this t and initially specified surviving rate 
  1-(plogis(qlogis(surv)-(N-600)*0.005)) #carying capacity with ~600 individuals total 
}



##### INITIALISATION #####

  population.total <- c() #empty vector for the population matrix
  statistic.matrix <- matrix(ncol=15, nrow=time) #empty array for the statistics
    
    
    for(k in 1:patches){ #LOOP OVER PATCHES
      patchx.N <- abs(round(rnorm(1, mean=300, sd=5))) #Number of individuals in the patch 
      patchx.male <- round(runif(1,patchx.N/4,3*patchx.N/4)) #Number of males in the patch
      
      ID <- c(1:(patchx.N)) #vector ID: gives each individual an ID
      patch <- c(rep(k,patchx.N)) #vector patch: gives each individual their patch Nr.
      gender <- c(rep("male",patchx.male),rep("female",patchx.N-patchx.male)) #vector gender: is filled with males and females
      trait <- c(rep(0.5,patchx.N)) #vector trait: is for all individuals from both patches set as 0.5
      female.trait <- c(rep(0.5,patchx.N))
      survival <- ceiling(runif(patchx.N, min=0, max=age)) #vector survival: randomly distributed between 1 and age limit 
      ID.mother <- c(rep(NA,patchx.N)) #the first generation has no mother and therefore no ID in the column for the mothers ID
      ID.father <- c(rep(NA,patchx.N)) #the first generation has no father and therefore no ID in the column for the fathers ID
      patch.last.year <- ceiling(runif(patchx.N, min=0, max=2)) #generates randomly ID of last years patch for each individual (patch 1 or 2)
      nr.offspring <- c(rep(0,patchx.N)) #number of offspring in first generation, will be filled with males success/offspring from last year
      terr <- c(rep(0, patchx.N)) #here the obtained territory is stored, emptied every t
      repro <- c(rep(0, patchx.N)) #decision stored if male can reproduce this t or not (1=True, 0=False)
      patch.born <- patch 
      loci <- c(1:40) #empty space for loci (nr of loci=40) 20 for male quality trait + 20 for female philopatry trait value
      
      patchx <- data.frame(ID,patch,gender,trait,female.trait,survival,ID.mother,ID.father, patch.last.year, nr.offspring, terr, repro, patch.born) #the dataframe is constructed for each patch including all vectors which where defined just before
      loci.matrix.pop <- matrix(ncol=length(loci.col), nrow=patchx.N) 
      patchx <- cbind(patchx, loci.matrix.pop)
      population.total <- rbind(population.total,patchx)  #data frame including all individuals of all patches (the dataframe of a patch is included in the population matrix)
    }
    
    population.total$ID <- c(1:nrow(population.total)) #the first generation of the population becomes a new ID
    ID.scan <- nrow(population.total)+1
    
    ##### STATISTIC START #####
    population.N <- rep(0,time) #empty vector for the populationsize of each generation (includes also pending males...)
    population.N1 <- rep(0,time) #empty vector for the pop size in patch 1 of each generation
    population.N2 <- rep(0,time) #empty vector for the pop size in patch 2 of each generation
    population.meantrait1.males <- rep(0,time) #empty vector for the mean trait in patch 1 of each generation
    population.meantrait2.males <- rep(0,time) #empty vector for the mean trait in patch 2 of each generation
    population.meantrait1.females <- rep(0,time) #empty vector for the mean trait in patch 1 of each generation
    population.meantrait2.females <- rep(0,time) #empty vector for the mean trait in patch 2 of each generation
    population.males1 <- rep(0,time) #empty vector for the number of males in patch 1 of each generation
    population.males2 <- rep(0,time) #empty vector for the number of males in patch 2 of each generation
    population.females1 <- rep(0,time) #empty vector for the number of females in patch 1 of each generation
    population.females2 <- rep(0,time) #empty vector for the number of females in patch 2 of each generation
    offspring.produced1 <- rep(0, time) #empty vector for number of offspring produced in patch 1
    offspring.produced2 <- rep(0, time) #empty vector for number of offspring produced in patch 2
    cov.males1 <- rep(0,time) #empty vector for covariance of number of offspring and male quality in patch 1
    cov.males2 <- rep(0,time) #empty vector for covariance of number of offspring and male quality in patch 2
    
    ########STATISTIC END  #####
    
    population <- nrow(population.total) #number of individuals
    
    for(x in 1:population){ #LOOP OVER THE INDIVIDUALS
      population.total[x,loci.col] <- ceiling(runif(40,1e-16,10)) #each individual has 40 random numbers (first 10:row //last 10:column), the first 20 are for male trait, the last 20 for female trait
    }
    
    loci.matrix <- population.total[,loci.col] #get all loci from current pop matrix 
    population.total <- trait.fun(population.total,loci.matrix, gen_phen_map, gen_phen_map2) #traitvalue-function: traitvalues for the population are included and overwrite the population matrix
    #population.total <- female.trait.fun(population.total,values.population,loci.matrix, gen_phen_map2) #traitvalue-function: traitvalues for the population are included and overwrite the population matrix
    
    ##### GENERATION LOOP START #####  
    for(t in 1:time){
      
      N.last1 <- nrow(population.total[population.total$patch==1,]) #storing patch 1 N for previous year
      N.last2 <- nrow(population.total[population.total$patch==2,]) #storing patch 2 N for previous year
      N <- nrow(population.total) #number of individuals in total (all patches included)
      
      if(N>0) { #START IS ANYBODY THERE-LOOP: if there are any individuals and the population is not extinct 
    
        ##### WHICH MALES ARE READY TO COMPETE? ####
        population.total$repro <- 0
        population.total$repro[population.total$gender=="female"] <- 1 #all females are able to reproduce
        population.total$repro[population.total$survival<(age-3)&population.total$gender=="male"] <- 1 #males that are old enough get a 1 to make sure they can compete and reproduce afterwards, will be changed when they loose fight (dont obtain a territory)
        
        if(any(population.total$gender=="male" & population.total$repro>0)){ #are there any reproductively active males?
          
        ##### 
        
        N.male <- population.total[population.total$gender=="male" & population.total$repro==1,] #get all male individuals as new matrix
        population.males <- nrow(N.male) #number of male individuals
        
        ##### MALE PATCH CHOICE #####
        
        N.male <- choice.fun(N.male, patches) #Males decide where to go this year depending on last years success
        population.total[population.total$gender=='male'&population.total$repro==1,]$patch <- N.male$patch #add info to population matrix
        
        population.total[population.total$gender=='male'&population.total$repro==1,]$nr.offspring <- 0 #set number of offspring to zero, so that number of offspring in this t can be added after reproduction again
        
        ##### MALE PATCH CHOICE END #####
        
        ##### MALE COMPETITION - WHICH TERRITORY MALE ESTABLISH/OBTAIN #####
        
        population.total$terr <- 0 #empty the territory vector for all indivduals
        N.male <- competition.fun(N.male, patches, population.males, territories) #territories are obtained after competition of males 
        N.male <- N.male[order(N.male$ID),] #order ID's because in the comp. function the individuals are reorderd and not the same order as in male matrix before. Ordering ID's gets it back in previous order
        
        population.total[population.total$gender=='male' & population.total$repro==1,]$terr <- N.male$terr #obtained territories of "winners" are written into pop.matrix
        
        #All males that lost territory competition have certain mortality:
        population.total$survival[population.total$gender=="male" & population.total$terr==0 & population.total$repro==1 & runif(nrow(population.total)) < die.fight] <- 0
        
        #Update all population info after males died 
        population.total <- population.total[population.total$survival>0,] #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
        N.male <- population.total[population.total$gender=="male" & population.total$repro==1,]
        N <- nrow(population.total)
        
        ##### MALE COMPETITION - FIGHT FOR TERRITORIES II  #####
        
        #Let males that lost in previous fight switch to random other patch
        N.male$patch[N.male$terr==0] <- (N.male$patch[N.male$terr==0] - 1 + floor(runif(1,1,patches)))%%patches + 1
        population.total$patch[population.total$gender=='male' & population.total$repro==1] <- N.male$patch  #overwrite patch choice from before 
        
        #Males choose their territory again, fight again 
        N.male <- competition.fun(N.male, patches, population.males, territories) 
        N.male <- N.male[order(N.male$ID),] #order ID's because in the comp. function the individuals are reorderd and not the same order as in male matrix before. Ordering ID's gets it back in previous order
        population.total[population.total$gender=='male'&population.total$repro==1,]$terr <- N.male$terr
        
        #All males that lost territory competition have certain mortality:
        population.total$survival[runif(nrow(population.total)) < die.fight & population.total$gender=="male" & population.total$terr==0 & population.total$repro==1] <- 0 # die.fight #for each individual is a random number distributed. if the number is below the deathrate a true is written into the vector + ID
        
        #Update all population info after males died 
        population.total <- population.total[population.total$survival > 0,] #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
        N <- nrow(population.total)
       
        ##### CHOOSE MALES FOR REPRODUCTION ####
        
        population.total$repro[population.total$terr>0&population.total$gender=="male"] <- 1 #males that obtained territory during competition are able to reproduce
        population.total$repro[population.total$terr==0&population.total$gender=="male"] <- 0 #males that did not obtain territory during competition are not able to reproduce
        N.male <- population.total[population.total$gender=="male" & population.total$repro==1,]
        
        } else{
            N.male <- subset(population.total,population.total$gender=="male"&population.total$repro==1) #this happens when there are no males 
        }
        
        
        ##### COMPETITION END #####
       
        #### FEMALE PATCH CHOICE #### 
        
        N.female <- subset(population.total,population.total$gender=="female") #matrix with just female individuals
        if(nrow(N.female)>0){ #are there any females?
          N.female <- choice.fun.females(N.female,p,u,N.last1,N.last2, patches) #patch choice this year, depending on philopatry trait and last years density on birth patch
          patchbook_females <- N.female$patch
          population.total[population.total$gender=='female',]$patch <- patchbook_females #overwrite patch choice from before       
        }
        
        ### FEMALE CHOICE END ###
        
        #Check if male and female are in same patch for mating:
        tryst <- c(rep(0,patches))
        N.female <- c()
        
        for(pls in 1:patches){#in which patches are both males and females
          if(
            nrow(subset(population.total,population.total$patch==pls&population.total$gender=="male"&population.total$repro==1))>0 &
            nrow(subset(population.total,population.total$patch==pls&population.total$gender=="female"))>0
          ){
            tryst[pls]<-2
          }
        }
        
        for(neko in 1:patches){#all females which have males in their patches
          if(tryst[neko]>0){
            N.female<-rbind(N.female,subset(population.total,population.total$gender=="female"&population.total$patch==neko))
          }
        }
        
        
        # if(max(tryst)==2){ #IS OFFSPRING POSSIBLE? If one patch contains both genders then tyst has a level of 2
          N.0 <- N/500
          offspring.vector <- 0
          
          if(nrow(N.female)>0){ #number of offspring per female
            offspring.vector <- rep(1,nrow(N.female)) #each female gets one pup 
          }
          
          ID.offspring <- c() #empty vector for the ID of the offspring
          patch.offspring <- c() #empty vector for the patch of the offspring
          gender.offspring <- c() #empty vector for the gender of the offspring
          trait.offspring <- c() #empty vector for the trait of the offspring
          female.trait.offspring <- c()
          survival.offspring <- c() #each offspring gets the survival of the maximum age
          ID.mother.offspring <- c() #empty vector for the mothers ID of the offspring
          ID.father.offspring <- c() #empty vector for the fathers ID of the offspring
          
          loci.offspring <- matrix(NA,nrow=sum(offspring.vector),ncol=40) #empty vector for the locis of the offspring
          
          #### START LOOP PARTNERFINDING #####
          population.offspring <- c() #empty vector for the patchnumber of the offspring
          genderbook <- c() #empty vector for the gender of the offspring
          
          N.female.patch <- table(factor(N.female$patch,levels = 1:patches)) #number of females in each patch (as a vector)
          N.male.patch <- table(factor(N.male$patch,levels = 1:patches))#number of males in each patch (as a vector)
          
          current.offspring <- 1 #counter that keeps track of how much offspring have emerged so far during the loop below
          
          # restructure, perform actions per patch, rahter then per female (might allow for more vectorisation)
          for(pat in 1:patches){
            if(N.female.patch[[pat]] > 0 & N.male.patch[[pat]] > 0){ # no need to separately check whether the number of offspring > 0, since every female obtains 1 pup
              N.fat <- sum(offspring.vector[N.female$patch==pat] > 0) # the number of required fathers is equal to the number of nonzero litters
              fat.row <- sample.vec(1:nrow(N.male),N.fat,replace=TRUE,prob=as.numeric(N.male$patch==pat)) # rownumbers of the selected fathers
              fat.row <- rep(fat.row,times=offspring.vector[N.female$patch==pat])
              mot.row <- rep(1:nrow(N.female),times=as.numeric(N.female$patch==pat)*offspring.vector)
              ID.father.offspring <- c(ID.father.offspring,N.male$ID[fat.row])
              offs <- N.female[mot.row,]  # initially, offspring are basically their mother
              offs$ID.father <- N.male$ID[fat.row]
              offs$ID.mother <- offs$ID
              offs$gender <- sample(c('male','female'),nrow(offs),replace=TRUE)
              offs[,c('patch.born','patch','patch.last.year')] <- pat
              offs[,c('trait','female.trait','terr','repro','nr.offspring')] <- 0
              offs$survival <- age
              offs$ID <- ID.scan:(ID.scan+nrow(offs)-1) # resetting some of the columns
              # cat("t = ",t,"\t","ID.scan = ",ID.scan,"\n")
              ID.scan <- ID.scan + nrow(offs)
              # cat("t = ",t,"\t","ID.scan = ",ID.scan,"\n")
              # flush.console()
              fat.dum <- expand.grid(loc=c(1:10,21:30), fat = fat.row)
              mot.dum <- expand.grid(loc=c(1:10,21:30), mot = mot.row)
              
              fat.dum$loc <- loci.col[fat.dum$loc + sample(c(0,10),nrow(fat.dum),replace=TRUE)]
              mot.dum$loc <- loci.col[mot.dum$loc + sample(c(0,10),nrow(mot.dum),replace=TRUE)]
              
              off.loc.f <- matrix(as.numeric(N.male[as.matrix(fat.dum[,2:1])]),ncol=20,byrow = TRUE)
              off.loc.m <- matrix(as.numeric(N.female[as.matrix(mot.dum[,2:1])]),ncol=20,byrow = TRUE)
              
              offs[,loci.col] <- cbind(off.loc.m[,1:10],off.loc.f[,1:10],off.loc.m[,11:20],off.loc.f[,11:20])

              # mutate
              mut_loc <- matrix(runif(nrow(offs)*length(loci.col)) < mutate,ncol=length(loci.col))
              if(any(mut_loc)){
                # print(which(mut_inds))
                offs[,loci.col][mut_loc] <- sample(1:10,sum(mut_loc),replace=TRUE)
              }
              
              test_fun <- function(i){
                gen.own <- offs[i,loci.col]
                fat.source <- N.male[fat.row[i],loci.col]
                mot.source <- N.female[mot.row[i],loci.col]
                gen.source <- cbind(rbind(as.numeric(mot.source[,1:10]),as.numeric(mot.source[,11:20])),
                rbind(as.numeric(fat.source[,1:10]),as.numeric(fat.source[,11:20])),
                rbind(as.numeric(mot.source[,21:30]),as.numeric(mot.source[,31:40])),
                rbind(as.numeric(fat.source[,21:30]),as.numeric(fat.source[,31:40])))
                all(gen.own == gen.source[1,] | gen.own == gen.source[2,])
              }
              # test_fun(13)
              # loci.col
              population.offspring <- rbind(population.offspring,offs)
            }
          }
          population.offspring <- trait.fun(population.offspring, population.offspring[,loci.col], gen_phen_map, gen_phen_map2) #the offspring matrix is overwritten including the traitvalues calculated by the traitvalue-function
                
          #INFATICIDE: Let offspring die with mortality depending on patch density
          densities <- N.female.patch + N.male.patch
          population.offspring <- population.offspring[!runif(nrow(population.offspring)) < i+s*plogis(0.01*densities[population.offspring$patch]),]
              
          population.total <- rbind(population.total,population.offspring) #the offspring population matrix is added to the general population matrix
          rownames(population.total) <- 1:nrow(population.total) #rownames are overwritten
        
        population.total$nr.offspring <- table(factor(ID.father.offspring,levels=population.total$ID)) #writing the number of offspring into the nr_offspring columns, stored for one t
        
        ##### DEATH START #####
        #death by age:
        N <- nrow(population.total)
        
        population.total$survival <- population.total$survival-1 #every adult loses one survival counter per generation
        
        #random Death:
        die <- mortality(N,surv) #get density dependend mortality rate (=die)
        
        dying.individuals <- runif(nrow(population.total),0,1) < die #for each individual is a random number distributed. if the number is below the deathrate the individual is written into a vector
        population.total$survival[dying.individuals] <- 0 #the individuals that where written into the vector below, become a 0 in their survival
        
        #erasing dead individuals:
        if(nrow(population.total)>0){
        population.total <-subset(population.total,population.total$survival>0) #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
        }
        
        
        ##### END DEATH #####   
        ###Statistic 2##
        population.N[t] <- nrow(population.total) #overwrites the populationsizes for each generation in the empty vector
        population.N1[t] <- nrow(population.total[population.total$patch==1&population.total$repro==1,]) #get population size from patch 1 for all individuals that reproduced
        population.N2[t] <- nrow(population.total[population.total$patch==2&population.total$repro==1,]) #get population size from patch 2  for all ind. that reproduced 
        population.meantrait1.males[t] <- mean(population.total[population.total$gender=="male"&population.total$patch==1&population.total$repro==1,]$trait)  #average trait-value from males for patch 1  for first generation
        population.meantrait2.males[t] <- mean(population.total[population.total$gender=="male"&population.total$patch==2&population.total$repro==1,]$trait) #average trait-value from males for patch 2  for first generation
        population.meantrait1.females[t] <- mean(population.total[population.total$gender=="female"&population.total$patch==1&population.total$repro==1,]$female.trait)  #average trait-value from females for patch 1  for first generation
        population.meantrait2.females[t] <- mean(population.total[population.total$gender=="female"&population.total$patch==2&population.total$repro==1,]$female.trait) #average trait-value from females for patch 2  for first generation
        population.males1[t] <- nrow(population.total[population.total$gender=="male"&population.total$patch==1&population.total$repro==1,]) #Number of males in patch 1  for first generation
        population.males2[t] <- nrow(population.total[population.total$gender=="male"&population.total$patch==2&population.total$repro==1,]) #Number of males in patch 2  for first generation
        population.females1[t] <- nrow(population.total[population.total$gender=="female"&population.total$patch==1,]) #Number of females in patch 1  for first generation
        population.females2[t] <- nrow(population.total[population.total$gender=="female"&population.total$patch==2,]) #Number of females in patch 2  for first generation
        offspring.produced1[t] <- nrow(population.total[population.total$survival==(age-1)&population.total$patch==1,])#number of new offspring in patch 1 
        offspring.produced2[t] <- nrow(population.total[population.total$survival==(age-1)&population.total$patch==2,])#number of new offspring in patch 2
        cov.males1[t] <- cov((population.total[population.total$gender=="male"&population.total$patch==1&population.total$repro==1,]$nr.offspring),(population.total[population.total$gender=="male"&population.total$patch==1&population.total$repro==1,]$trait)) #covariance of number of offspring and male quality in patch 1
        cov.males2[t] <- cov((population.total[population.total$gender=="male"&population.total$patch==2&population.total$repro==1,]$nr.offspring),(population.total[population.total$gender=="male"&population.total$patch==2&population.total$repro==1,]$trait)) #covariance of number of offspring and male quality in patch 2
        
        statistic.matrix[t,] <- cbind(population.N[t],population.N1[t],population.N2[t],population.meantrait1.males[t], population.meantrait2.males[t],population.meantrait1.females[t], population.meantrait2.females[t], population.males1[t], population.males2[t], population.females1[t], population.females2[t], offspring.produced1[t], offspring.produced2[t],  cov.males1[t],  cov.males2[t])
        
        ##### End Statistic 2#############
        
      }#END IS ANYBODY THERE? 
      #print(t)
    }##### END GENERATION LOOP #####
    
    #Stored summary statistic formatted for output data
    statistic.matrix[is.na(statistic.matrix)] <- 0 #NaN can be produced when trait values are not existing (remove these and call them 0)
    colnames(statistic.matrix) <- c("N","N1","N2","meantrait.males1","meantrait.males2","meantrait.females1","meantrait.females2","N.males1","N.males2", "N.females1", "N.females2", "offspring.produced1", "offspring.produced2", "cov.males1", "cov.males2") #column names of statistic store matrix
    return(statistic.matrix)
}#END SIMULATION.RUN

#Run function
#debug(simulation.fun)
# library(profvis)
# # profvis({
# statistic <- simulation.fun()
# })




