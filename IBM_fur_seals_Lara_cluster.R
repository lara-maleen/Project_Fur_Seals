#Minimal Model IBM - Fur Seal Project
#Hypothesis 1: Marginal-Male Theory
#Code for running in cluster (without replicates and storing information in matrix every t)
rm(list=ls())


##### START SIMULATION.RUN-FUNCTION #####
# simulation.fun <- function(
  time=100#, #number of generations
                           age=15#, #age limit for an individual; life span for A. gazella. 15-25 years --> literature!?
                           patches=2#, #number of Patches (two different sites: high/low density)
                           territories=c(50,50)#, #number of territories per patch
                           mutate=0.05#, #mutation factor
                           die.fight=0.35#, #propability to die from fight/competition
                           loci.col=c(14:53)#, #loci column numbers of the pop matrix 
                           p= 0.25#, #parameter for philopatry function (female patch choice) -> the higher p is, the more intense the philopatric (side-fidelity) influence
                           u = 100#, #assumed normal average density (for each patch), used for female patch choice function
                           i=-0.8#, #intercept for infanticide function
                           s=1.8#, #slope for infanticide function
                           surv=0.75#, #survival for total population 
                           gene_file1="genes.rds"#,
                           gene_file2="genes2.rds"
# ){

profvis({
  
#gen_phen_map <- readRDS('/data/home/lara/genes.rds') #load the gene array (10 loci, 10 alleles) #gene map used in cluster
#gen_phen_map2 <- readRDS('/data/home/lara/genes2.rds') #load the gene array (10 loci, 10 alleles) #gene map used in cluster
gen_phen_map <- readRDS(gene_file1) #load the gene array (10 loci, 10 alleles), used for male trait values
gen_phen_map2 <- readRDS(gene_file2) #second gene map for female trait value (10 loci, 10 alleles). Phenotype of -0.2 and +0.2 initially
gen_phen_map2 <- array(round(runif(1e3),2),dim=c(10,10,10))
##### FUNCTIONS #####


ID.fun <- function(offspring.vector){ #ID-FUNCTION: for each individual a new ID is created
    ID.offspring <-   ID.scan:(ID.scan+sum(offspring.vector)-1)
    ID.scan <<- ID.scan + sum(offspring.vector)
    return(ID.offspring)
}
  

trait.fun <- function(population.total,value.matrix, loci.matrix, gen_phen_map, gen_phen_map2){ #TRAIT-VALUE-FUNCTION - used for male quality + female philopatry trait 

  #Male Trait Value
  # value.matrix <- matrix(NA,nrow(population.total),ncol=10) #empty matrix for the trait values for each loci
  # for(y in 1:nrow(population.total)){ #for each individual
  #   for(z in 1:10){ #for number of loci
  #     value.matrix[y,z] <- (gen_phen_map[loci.matrix[y,z],loci.matrix[y,10+z],z]) #get value from gene map 1 (this is the male trait gene map), go through all loci and see what alleles individual have
  #   }
  #   population.total[y,4] <- abs(sum(value.matrix[y,])) #calculate additive phenotypic trait value, stored in column number 4 (male trait value)
  # }
  # 
  
  # for the male trait value
  lc1 <- as.numeric(t(loci.matrix[,1:10]))
  lc2 <- as.numeric(t(loci.matrix[,11:20]))
  zs <- rep(1:10,nrow(population.total))
  phen <- gen_phen_map[cbind(lc1,lc2,zs)]
  population.total[,4] <- colSums(matrix(phen,nrow=10))
  #                                 
  # if(!all(colSums(matrix(phen,nrow=10)) -population.total[,4] == 0)){
  # stop("help")
  # }
  
  #Female Trait Value:
  lc1 <- as.numeric(t(loci.matrix[,21:30]))
  lc2 <- as.numeric(t(loci.matrix[,31:40]))
  zs <- rep(1:10,nrow(population.total))
  phen <- gen_phen_map2[cbind(lc1,lc2,zs)]
  population.total[,5] <- colSums(matrix(phen,nrow=10))
  
  # value.matrix <- matrix(NA,nrow(population.total),ncol=10) #empty matrix for the trait values for each loci
  # for(y2 in 1:nrow(population.total)){ #for each individual
  #   for(z2 in 1:10){ #for each loci 
  #     value.matrix[y2,z2] <- gen_phen_map2[loci.matrix[y2,z2+20],loci.matrix[y2,10+z2+20],z2] #get value from gene map 2 (female trait gene map), loci columns 21-40 in pop matrix (i.e. loci matrix 21-40)
  #   }
  #   population.total[y2,5] <- (sum(value.matrix[y2,])) #calculate additive phenotypic trait value, stored in column number 5 (female trait value)
  # }
  # 
  # print(all(colSums(matrix(phen,nrow=10)) -population.total[,5] == 0))
  # stop("help")
  # }
  return(population.total) 
}


choice.fun <- function(N.male, patches){ #MALE PATCH CHOICE: decide where each adult male goes to this t
  for(i in 1:nrow(N.male)){ #for each male 
    if (N.male$nr.offspring[i]>0){ #If reproductive success (offspring number) is greater than 0 patch stays the same (from last year)
      N.male$patch[i] <- N.male$patch.last.year[i]
    }
    else{ #Otherwise the patch is changed in contrary patch 
       N.male$patch[i] <- (N.male$patch[i] - 1 + floor(runif(1,1,patches)))%%patches + 1
      }
    }
  return(N.male) #New patch is written into patch column in male matrix
}


choice.fun.females <- function(N.female,p,u,N.last1,N.last2, patches){ #FEMALE PATCH CHOICE: determines the patch for females, depending on last years N (from last years patch) as well as density-preference trait & philopatry
  
  N.last <- c(N.last1,N.last2) #get the population size from the previous year per patch
  
  #Add philopatric decision (parameter set at the beginning)
  p.patch <- p > runif(nrow(N.female),0,1) #decide wether female is philopatric or not (stays at birth patch = TRUE), if not then the density-dependent choice takes place
  
  #on the positions where p.patch is TRUE, the patch number is the birth patch:
  for (i in 1:length(p.patch)){
    
    if (p.patch[i]){ #if this is true, than female go to the patch it was born
      N.female$patch[i] <- N.female$patch.born[i]
    }
    
    else{ #otherwise the female gets a new TRUE or FAlSE depending on the density of last years patch 
      patch.u  <- plogis(N.female$female.trait[i]*(N.last[N.female$patch[i]] - u)) > runif(1,0,1)
      if ((patch.u)){ #if that is true, the patch is changed to the other patch  
        N.female$patch[i] <- (N.female$patch[i] - 1 + floor(runif(1,1,patches)))%%patches + 1
      }
    }
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
    population.total <- trait.fun(population.total,values.population,loci.matrix, gen_phen_map, gen_phen_map2) #traitvalue-function: traitvalues for the population are included and overwrite the population matrix
    #population.total <- female.trait.fun(population.total,values.population,loci.matrix, gen_phen_map2) #traitvalue-function: traitvalues for the population are included and overwrite the population matrix
    
    ##### GENERATION LOOP START #####  
    for(t in 1:time){
      
      N.last1 <- nrow(population.total[population.total$patch==1,]) #storing patch 1 N for previous year
      N.last2 <- nrow(population.total[population.total$patch==2,]) #storing patch 2 N for previous year
      N <- nrow(population.total) #number of individuals in total (all patches included)
      
      if(N>0) { #START IS ANYBODY THERE-LOOP: if there are any individuals and the population is not extinct 
    
        ##### WHICH MALES ARE READY TO COMPETE? ####
        
        if(nrow(population.total[population.total$gender=="female",])>0){ #are there any females?
        population.total[population.total$gender=="female",]$repro <- 1 #all females are able to reproduce
        }
        
        if(nrow(population.total[population.total$gender=="male",])>0){ #are there any males?
          
          if(nrow(population.total[population.total$gender=="male"&population.total$survival<(age-3),])>0){ #are there any males that are over 3 years old --> Hoffman 2003 'MALE REPRODUCTIVE STRATEGY AND THE IMPORTANCE OF MATERNAL STATUS IN THE ANTARCTIC FUR SEAL ARCTOCEPHALUS GAZELLA'
          population.total[population.total$survival<(age-3)&population.total$gender=="male",]$repro <- 1 #males that are old enough get a 1 to make sure they can compete and reproduce afterwards, will be changed when they loose fight (dont obtain a territory)
          if(nrow(population.total[population.total$gender=="male"&population.total$survival>=(age-3),])>0){
          population.total[population.total$survival>=(age-3)&population.total$gender=="male",]$repro <- 0 #males that are not old enough get a 0 to make sure they cannot compete and reproduce afterwards, will be changed when they loose fight (dont obtain a territory)
        
        ##### 
        
        N.male <- subset(population.total,population.total$gender=="male"&population.total$repro==1) #get all male individuals as new matrix
        population.males <- nrow(N.male) #number of male individuals
        
        ##### MALE PATCH CHOICE #####
        
        patchbook_males <- c() #vector for storing the patch choice of males
        N.male <- choice.fun(N.male, patches) #Males decide where to go this year depending on last years success
        patchbook_males <- N.male$patch #overwrite patch from previous year
        population.total[population.total$gender=='male'&population.total$repro==1,]$patch <- patchbook_males #add info to population matrix
        
        population.total[population.total$gender=='male'&population.total$repro==1,]$nr.offspring <- rep(0,nrow(N.male)) #set number of offspring to zero, so that number of offspring in this t can be added after reproduction again
        
        ##### MALE PATCH CHOICE END #####
        
        ##### MALE COMPETITION - WHICH TERRITORY MALE ESTABLISH/OBTAIN #####
        
        population.total$terr <- c(rep(0, nrow(population.total))) #empty the territory vector for all indivduals
        terrbook_males <- c() #vector for storing the territory choice for each male
        N.male <- competition.fun(N.male, patches, population.males, territories) #territories are obtained after competition of males 
        N.male <- N.male[order(N.male$ID),] #order ID's because in the comp. function the individuals are reorderd and not the same order as in male matrix before. Ordering ID's gets it back in previous order
        terrbook_males <- N.male$terr
        population.total[population.total$gender=='male'&population.total$repro==1,]$terr <- terrbook_males #obtained territories of "winners" are written into pop.matrix
        
        #All males that lost territory competition have certain mortality:
        dying.males <- matrix(NA,nrow(population.total[population.total$gender=="male"&population.total$terr==0&population.total$repro==1,]),ncol=2) #empty matrix for all losers
        dying.males[,2] <- runif(nrow(population.total[population.total$gender=="male"&population.total$terr==0&population.total$repro==1,]),0,1) < die.fight #for each individual is a random number distributed. if the number is below the deathrate a true is written into the vector + ID
        dying.males[,1] <- population.total[population.total$gender=="male"&population.total$terr==0&population.total$repro==1,]$ID #IDS of the males are written here 
        dying.males.ID <- c()
        dying.males.ID <- dying.males[,1][dying.males[,2]==1] #IDs of the males that died are stored
        population.total[population.total$ID %in% dying.males,]$survival <- 0

        
        #Update all population info after males died 
        population.total <-subset(population.total,population.total$survival>0) #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
        N.male <- subset(population.total,population.total$gender=="male"&population.total$repro==1)
        N <- nrow(population.total)
        
        ##### MALE COMPETITION - FIGHT FOR TERRITORIES II  #####
        
        #Let males that lost in previous fight switch to other patch
        males.patch.shift <- N.male[N.male$terr==0,]$ID #get the males that didnt obtain territory, they shift patches (ID is safed) 
        N.male$patch[N.male$terr==0] <- (N.male$patch[N.male$terr==0] - 1 + floor(runif(1,1,patches)))%%patches + 1
        
        patchbook_males <- c()
        patchbook_males <- N.male$patch 
        population.total[population.total$gender=='male'&population.total$repro==1,]$patch <- patchbook_males #overwrite patch choice from before 
        
        #Males choose their territory again, fight again 
        
        terrbook_males <- c()
        N.male <- competition.fun(N.male, patches, population.males, territories) 
        terrbook_males <- N.male$terr
        population.total[population.total$gender=='male'&population.total$repro==1,]$terr <- terrbook_males
        
        #All males that lost territory competition have certain mortality:
        dying.males <- matrix(NA,nrow(population.total[population.total$gender=="male"&population.total$terr==0&population.total$repro==1,]),ncol=2) #empty matrix for all losers
        dying.males[,2] <- runif(nrow(population.total[population.total$gender=="male"&population.total$terr==0&population.total$repro==1,]),0,1) < die.fight #for each individual is a random number distributed. if the number is below the deathrate a true is written into the vector + ID
        dying.males[,1] <- population.total[population.total$gender=="male"&population.total$terr==0&population.total$repro==1,]$ID #IDS of the males are written here 
        dying.males.ID <- c()
        dying.males.ID <- dying.males[,1][dying.males[,2]==1] #IDs of the males that died are stored
        for(d2 in dying.males.ID){ #go trough the died males and change survival number and the loci matrix 
          population.total[population.total$ID==d2,]$survival <- 0
        }
        
        #Update all population info after males died 
        population.total <-subset(population.total,population.total$survival>0) #population matrix: Individuals which have a survival higher then 0 stay alive in the dataframe. the others are deleted
        N.male <- subset(population.total,population.total$gender=="male"&population.total$repro==1)
        N <- nrow(population.total)
       
        ##### CHOOSE MALES FOR REPRODUCTION ####
        
        population.total[population.total$terr>0&population.total$gender=="male",]$repro <- 1 #males that obtained territory during competition are able to reproduce
        population.total[population.total$terr==0&population.total$gender=="male",]$repro <- 0 #males that did not obtain territory during competition are not able to reproduce
        N.male <- subset(population.total,population.total$gender=="male"&population.total$repro==1) #change male matrix 
        
        }
          else{
            N.male <- subset(population.total,population.total$gender=="male"&population.total$repro==1) #this happens when there are no males 
          }
        }
          else{
            N.male <- subset(population.total,population.total$gender=="male"&population.total$repro==1) #this happens when there are no males over 3 years
          }
        
        }  #End are there any males for fight?
        else{
          N.male <- subset(population.total,population.total$gender=="male"&population.total$repro==1) #this happens when there are no males under 3 years
        }
       
        
        ##### COMPETITION END #####
       
        #### FEMALE PATCH CHOICE #### 
        
        N.female <- subset(population.total,population.total$gender=="female") #matrix with just female individuals
        if(nrow(N.female)>0){ #are there any females?
        patchbook_females <- c()
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
        
        
        if(max(tryst)==2){ #IS OFFSPRING POSSIBLE? If one patch contains both genders then tyst has a level of 2
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
          patchbook <- c() #empty vector for the patchnumber of the offspring
          genderbook <- c() #empty vector for the gender of the offspring
          
          N.female.patch <- table(factor(N.female$patch,levels = 1:patches)) #number of females in each patch (as a vector)
          N.male.patch <- table(factor(N.male$patch,levels = 1:patches))#number of males in each patch (as a vector)
          
          current.offspring <- 1 #counter that keeps track of how much offspring have emerged so far during the loop below
          
          if(nrow(N.female)>0){ #START ANY FEMALES?: loop starts if there is at least one female individual
            if(nrow(N.male)>0){
              for(u in 1:nrow(N.female)){ #START LOOP PARTNERFINDING/mother 
                if(offspring.vector[u]>0){ #START GETS THE MOTHER OFFSPRING?
                  if(N.male.patch[N.female$patch[u]]>0){ #START ANY MALES IN THE PATCH OF THE MOTHER?: loop starts if there is at least one male individual in the mothers patch
                    mother <- N.female$ID[u] #gives the ID of the mother
                    ID.mother.offspring <- c(ID.mother.offspring, rep(mother,offspring.vector[u])) #ID of the mother is written into the vector for all her offspring
                  
                    ###FATHER####
                    
                    potfather <- N.male$ID[N.male$patch==N.female$patch[u]] # storing the id's of potential fathers
                    if(length(potfather) > 1){
                      father <- sample(N.male$ID[N.male$patch==N.female$patch[u]],1) #sample the ID of one male which patchnumber is the same as the patchnumber of the mother
                    }else{
                      father <- potfather
                    }
                    
                    ID.father.offspring <- c(ID.father.offspring,rep(father,offspring.vector[u])) #ID of the father is written into the vector as often as he becomes offspring with the mother
                    
                    #GENETICS:
                    loci.mother <- population.total[population.total$ID==mother,loci.col] #vector of locis of the mother
                    loci.father <- population.total[population.total$ID==father,loci.col] #vector of locis of the father
                    loci.child <- rep(0,length(loci.col)) #empty vector with fixed length for the locis of the offspring
                    
                    for(o in 1:offspring.vector[u]){ #START LOOP NUMBER CHILDREN per female
                      loci.child[1:10] <- loci.mother[(1:10) +sample(c(0,10),10,replace=TRUE)] #the offspring becomes 10 locis for male trait sampled from the mother (this is the allele row)
                      loci.child[11:20] <- loci.father[(1:10) +sample(c(0,10),10,replace=TRUE)] #the offspring becomes 10 locis sampled for male trait from the father (this is allele column)
                      loci.child[21:30] <- loci.mother[(1:10) +sample(c(0,10),10,replace=TRUE)] #the offspring becomes 10 locis sampled for female trait from the mother (this is the allele row)
                      loci.child[31:40] <- loci.father[(1:10) +sample(c(0,10),10,replace=TRUE)] #the offspring becomes 10 locis sampled for female trait from the father (this is allele column)
                      #total sum of loci = 40
                      loci.child <- unlist(loci.child)
                      
                      #MUTATION
                      if(runif(1,0,1) < mutate){ #if a random number is lower than the mutationrate the offspring becomes a random distributed loci
                        loci.child[round(runif(1,1,40))] <- round(runif(1,1,10)) #the first runif selects the mutated loci and the second runif determines the new loci value (1-10 alleles)
                      }
                      
                      loci.child <- unlist(loci.child)
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
              
                loci.offspring <- as.data.frame(loci.offspring)
                patchbook <- rep(N.female$patch,offspring.vector) #each offspring becomes the patchnumber of the mother
                ID.offspring <- ID.fun(offspring.vector) #the ID of the offspring is calculated by the ID-function and written into the vector for their ID
                trait.offspring <- c(rep(0,length(patchbook))) #the traitvalue of the offspring is set to 0 for the moment
                female.trait.offspring <- c(rep(0,length(patchbook))) #the traitvalue of the offspring is set to 0 for the moment
                survival.offspring <- c(rep(age,length(patchbook))) #each offspring gets the survival of the age limit pre defined
                gender.offspring <- genderbook #genders of the offspring are written into the matrix
                patch.offspring <- patchbook #patches of offspring are written into the matrix
                nr.offspring.offspring <-  c(rep(0,length(patchbook))) #empty column for the subsequent offspring they will get
                terr.offspring <- c(rep(0, length(patchbook))) #empty column for subsequent territory
                repro.offspring <- c(rep(0, length(patchbook))) #empty column for subsequent decision for reproduction in t
                patch.born.offspring <- rep(N.female$patch,offspring.vector) #this stores info where on individuals is born (not changed over time)
                population.offspring <- data.frame(ID.offspring,patch.offspring,gender.offspring,trait.offspring, female.trait.offspring, survival.offspring,ID.mother.offspring,ID.father.offspring, patch.offspring, nr.offspring.offspring, terr.offspring, repro.offspring, patch.born.offspring) #a new dataframe is made for the offspring of this generation
                colnames(population.offspring) <- c("ID","patch","gender","trait","female.trait","survival","ID.mother","ID.father", "patch.last.year", "nr.offspring","terr","repro", "patch.born") #column names of the dataframe
                
                population.offspring <- cbind(population.offspring, loci.offspring)
                colnames(population.offspring) <- c("ID","patch","gender","trait","female.trait","survival","ID.mother","ID.father", "patch.last.year", "nr.offspring","terr","repro", "patch.born",1:40) #column names of the dataframe
                
                population.offspring <- trait.fun(population.offspring,values.offspring, loci.offspring, gen_phen_map, gen_phen_map2) #the offspring matrix is overwritten including the traitvalues calculated by the traitvalue-function
                #population.offspring <- female.trait.fun(population.offspring,values.offspring, loci.offspring, gen_phen_map2) #the offspring matrix is overwritten including the traitvalues calculated by the traitvalue-function
                
                #INFATICIDE: Let offspring die with mortality depending on patch density
                
                infanticide.vector <- c(rep(NA,patches))
                for(p4 in 1:patches){ #for each patch a specific mortality/infanticide rate, depending on density on the patch
                  curr_N <- sum(population.total$patch==p4)
                  y=i+s*plogis(0.01*(curr_N)) #mortality is created 
                  infanticide.vector[p4] <- y #safed in vector 
                }
                
                dying.offspring <- matrix(NA,nrow(population.offspring),ncol=2)
                dying.offspring[,2] <- runif(nrow(population.offspring),0,1) < infanticide.vector[population.offspring$patch] #for each individual is a random number distributed. if the number is below the deathrate the individual is written into a vector
                dying.offspring[,1] <- population.offspring$ID
                dying.offspring.ID <- c()
                dying.offspring.ID <- dying.offspring[,1][dying.offspring[,2]==1]
                population.offspring[population.offspring$ID %in% dying.offspring.ID,]$survival <- 0

                
                population.offspring <-subset(population.offspring,population.offspring$survival>0) #remove all dead offspring
                
                population.total <- rbind(population.total,population.offspring) #the offspring population matrix is added to the general population matrix
                rownames(population.total) <- 1:nrow(population.total) #rownames are overwritten
              
            }#END ANY MALES?
          }#END ANY FEMALES?
        }#END IS OFFSPRING POSSIBLE?
        
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
    # return(statistic.matrix)
 })  
# }#END SIMULATION.RUN

#Run function
#debug(simulation.fun)
# library(profvis)
# # profvis({
# statistic <- simulation.fun()
# })




