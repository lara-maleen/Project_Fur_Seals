# Faster function to calculate the offspring distribution for a female of class 'fem' when fathers are described by male.vals (the probability that a father passes on certain alleles)
# Input:
# fem: a single row from the 'dum' object, that describes the properties of the focal female (i.e. mother)
# dum: a data frame with information on the different classes in the population. Every row describes one class [same order as the popvect], furthermore it has columns:
#       - a1 and a2: allele at the a loci. Value of each can be either a or A
#       - b1 and b2: allele at the b loci. Value of each can be either b or B
# male.vals: a list with values mAB, mAb, maB and mab, each representing the probability that the offspring obtains a combination of (a/A) and (b/B) allele from the father
#
# output:
# a vector with the offspring distribution for a female of type 'fem'
calc_off_dist_alt <- function(fem,dum,male.vals){
  
  as.f <- 0.5*(as.numeric(fem$a1 == 'A') + as.numeric(fem$a2 == 'A')) # prob. of inheriting A from the mother
  bs.f <- 0.5*(as.numeric(fem$b1 == 'B') + as.numeric(fem$b2 == 'B')) # prob. of inheriting B from the mother
  

  # creating a 3x3 matrix with all potential combinations of scores and how likely they are to occur
  # Ascore 0,1,2 (with: 0 = aa, 1 = aA or Aa and 2 = AA) and Bscore 0,1,2 
  
  #            aa         Aa             AA
  #    bb
  #    Bb
  #    BB
  offs.score <- with(male.vals,{
    offs.score <- matrix(NA,nrow=3,ncol=3)
    offs.score[1,1] <- mab*(1-as.f)*(1-bs.f) # aabb
    offs.score[1,2] <- mAb*(1-as.f)*(1-bs.f) + mab*as.f*(1-bs.f) # Aabb
    offs.score[1,3] <- mAb*(1-bs.f)*as.f# AAbb
    offs.score[2,1] <- maB*(1-as.f)*(1-bs.f) + mab*(1-as.f)*bs.f # aabB
    offs.score[2,2] <- mAB*(1-as.f)*(1-bs.f) + mAb*(1-as.f)*bs.f + maB*as.f*(1-bs.f) + mab*as.f*bs.f # aAbB
    offs.score[2,3] <- mAB*as.f*(1-bs.f) + mAb*as.f*bs.f# AAbB
    offs.score[3,1] <- maB*(1-as.f)*bs.f # aaBB
    offs.score[3,2] <- mAB*(1-as.f)*bs.f + maB*as.f*bs.f# aABB
    offs.score[3,3] <- mAB*as.f*bs.f # AABB 
    offs.score
  })
  
  # create an offs vector that describes the resulting offspring distribution as a vector that has the same order as dum
  offs <- offs.score[cbind(dum$Bscore + 1,dum$Ascore + 1)] 
  
  # divide by two, to accout for the fact that the above procedure counts 1 male and 1 female offspring every time
  # instead of 0.5 male and 0.5 female offspring
  offs/2
}

# (old and slower) function to calculate the offspring distribution for a female of class 'fem'
# explicitly evaluates all combinations of fathers and mothers
calc_off_dist <- function(fem,dum,male.dist,A.adv){
  # print(fem)
  male.dist[dum$Ascore == 2] <- A.adv*male.dist[dum$Ascore == 2]
  male.dist <- male.dist/sum(male.dist)
  
  df.male <- rbind(data.frame(a1=dum$a1,b1=dum$b1,pfat=0.25*male.dist,stringsAsFactors = FALSE),
                   data.frame(a1=dum$a2,b1=dum$b1,pfat=0.25*male.dist,stringsAsFactors = FALSE),
                   data.frame(a1=dum$a1,b1=dum$b2,pfat=0.25*male.dist,stringsAsFactors = FALSE),
                   data.frame(a1=dum$a2,b1=dum$b2,pfat=0.25*male.dist,stringsAsFactors = FALSE))
  
  df.offs <- as.data.frame(rbind(df.male,df.male,df.male,df.male),stringsAsFactors=FALSE)
  
  colnames(df.offs) <- c('a1','b1','pfat')
  
  df.offs$a2 <- rep(rep(c(fem$a1,fem$a2),each=nrow(df.male)),2)
  df.offs$b2 <- rep(c(fem$b1,fem$b2),each=2*nrow(df.male))
  
  df.offs$pmat <- 0.25
  
  offs <- rep(0,length(male.dist))
  
  df.offs$Ascore <- as.numeric(df.offs$a1 == 'A') + as.numeric(df.offs$a2 == 'A')
  df.offs$Bscore <- as.numeric(df.offs$b1 == 'B') + as.numeric(df.offs$b2 == 'B')
  
  for(i in 1:nrow(df.offs)){
    loc <- which(df.offs$Ascore[i] == dum$Ascore & df.offs$Bscore[i] == dum$Bscore)
    offs[loc] <- offs[loc] + df.offs$pmat[i]*df.offs$pfat[i]/2
  }
  offs
}

female.dist <- function(dum,popvect,normalize=TRUE){
  # male dist isle 1
  female.dist.1 <- dum$p1*popvect*as.numeric(dum$sex=='f')
  
  # female dist isle 2
  female.dist.2 <- (1-dum$p1)*popvect*as.numeric(dum$sex=='f')
  
  
  if(sum(female.dist.1) > 0 & normalize){
    female.dist.1 <- female.dist.1 / sum(female.dist.1)
  }
  
  if(sum(female.dist.2) > 0 & normalize){
    female.dist.2 <- female.dist.2 / sum(female.dist.2)
  }
  
  
  return(list(female.dist.1,female.dist.2))
}

# Function for determining the male distribution on each island based on the population vector.
# input
# dum: a data frame with information on the different classes in the population. Every row describes one class [same order as the popvect], furthermore it has columns:
#       - sex: either 'm' or 'f' for male and female
#       - p1: the chance that an inidividual from a class goes to island 1
# popvect: vector of length nrow(dum), each entry describing the population size of each class
# maxfreq: the maximum fraction of the population that is allowed to be on one island (hard density constrained) should be between 0.5 and 1.
#          when one island has more individuals than allowed by maxfreq, the difference automatically migrates to the other island.
# normalize: if TRUE, the output distribution of males is normalized within each island.
#
# output:
# list of length 2. Each item of the list is the male distribution on one of the islands (i.e. a vector of length nrow(dum)). Possibly normalized within each island.
male.dist <- function(dum,popvect,maxfreq,normalize=TRUE){
  # male dist isle 1
  male.dist.1 <- dum$p1*popvect*as.numeric(dum$sex=='m')

  # male dist isle 2
  male.dist.2 <- (1-dum$p1)*popvect*as.numeric(dum$sex=='m')
  
  freq <- sum(male.dist.1) / sum(male.dist.1 + male.dist.2)
  if(maxfreq < 0.5){
    stop("the maximum frequency for males in either should be between 0.5 and 1")
  }
  
  if(freq > maxfreq){
    # a=freq-maxfreq: fraction of the male population that is too much in patch 1
    # b=sum(male.dist.1 + male.dist.2): total number of males
    # a*b: total number of males that should move
    # c=a*b/sum(male.dist.1): the fraction of males from patch 1 that has to move
    # c*male.dist.1: distribution of moving males
    diff <- sum(male.dist.1 + male.dist.2)*(freq-maxfreq)/sum(male.dist.1)*male.dist.1
    male.dist.1 <- male.dist.1 - diff
    male.dist.2 <- male.dist.2 + diff
  }else if(freq < (1-maxfreq)){
    diff <- sum(male.dist.1 + male.dist.2)*((1-freq)-maxfreq)/sum(male.dist.2)*male.dist.2
    male.dist.1 <- male.dist.1 + diff
    male.dist.2 <- male.dist.2 - diff
  }
  
  if(sum(male.dist.1) > 0 & normalize){
    male.dist.1 <- male.dist.1 / sum(male.dist.1)
  }
  
  if(sum(male.dist.2) > 0 & normalize){
    male.dist.2 <- male.dist.2 / sum(male.dist.2)
  }

    
  return(list(male.dist.1,male.dist.2))
}

# Function determining the male contributions to offspring, based on a given male distribution.
# input
# dum: a data frame with information on the different classes in the population. Every row describes one class [same order as the male.dist], furthermore it has columns:
#       - Ascore: integer, 0: aa, 1: aA or Aa, 2: AA
#       - Bscore: integer, 0: bb, 2: bB or Bb, 3: BB
# male.dist: a vector with the number of males per class
# A.adv: an AA male has an A.adv x higher chance of siring an offspring
# d: [0,1] dominance factor. The fitness of aan AA male is A.adv, of an Aa male A.adv^d and of an aa male 1.

# output
# a list with entriy mAB, mAb, maB and mab.
# each entry describes the chance for an offspring of obtaining these alleles when a father is sampled from the input distribution
male.vals <- function(dum,male.dist,A.adv,d){
  
  male.dist[dum$Ascore == 2] <- A.adv*male.dist[dum$Ascore == 2] # scale the AA males up by a factor A.adv to represent their higher quality
  male.dist[dum$Ascore == 1] <- A.adv^d*male.dist[dum$Ascore == 1]
  male.dist <- male.dist/sum(male.dist) # renormalize
  
  mAB <- sum(male.dist*0.5*dum$Ascore*0.5*dum$Bscore) # prob. of inheriting AB from father
  mAb <- sum(male.dist*0.5*dum$Ascore*(1-0.5*dum$Bscore)) # prob. of inheriting Ab from father
  maB <- sum(male.dist*(1-0.5*dum$Ascore)*0.5*dum$Bscore) # prob. of inheriting aB from father
  mab <- sum(male.dist*(1-0.5*dum$Ascore)*(1-0.5*dum$Bscore)) # prob. of inheriting ab from father
  
  list(mAB=mAB,mAb=mAb,maB=maB,mab=mab)
}

surv_adult <- function(dum,male.dists,female.dists){
  
  # on island 1, male survival depends on the amount of AA males
  # on island 2, all males survive equally
  
  # calculate sex and genotype specific survival on each island, based on density/frequency
  # N.m.1.AA <- sum(male.dists[[1]])
  #N.m.2 <- sum(male.dists[[2]])
  
  #N.f.1 <- sum(female.dists[[1]])
  #N.f.2 <- sum(female.dists[[2]])
  
  surv.f.1 <- 1
  surv.f.2 <- 1
  surv.m.1 <- dum$surv # density independent, but genotype determined survival
  surv.m.2 <- 1
  # print(surv.m.1)
  # finally: average survival based on actual distributions. i.e. 
  pre.surv <- ifelse(dum$sex=='f',(female.dists[[1]]*surv.f.1 + female.dists[[2]]*surv.f.2)/(female.dists[[1]]+female.dists[[2]]),(male.dists[[1]]*surv.m.1 + male.dists[[2]]*surv.m.2)/(male.dists[[1]]+male.dists[[2]]))
  
  pre.surv[male.dists[[1]] == 0 & male.dists[[2]] == 0 & dum$sex=='m'] <- 1
  pre.surv[female.dists[[1]] == 0 & female.dists[[2]] == 0 & dum$sex=='f'] <- 1
  # account for double zeroes in denominator!
  
  surv <- pre.surv # combine surv.f and surv.m

  # print(cbind(dum,surv))
  # stop("err")
  return(surv)
}

# Function that generates the matrix to advance the population vector to the next time step
# input:
# surv: survival of adults
# surv_off: function with one argument, density that returns [0,1] 1 - survival for offspring.
# popvect: the population distribution
# dum: a dataframe describing the classes of the population vector. It has length(popvect) rows and each column describes a different property of the classes, e.g.:
#       - sex: either 'm' or 'f' for male and female
#       - p1: the chance that an inidividual from a class goes to island 1
# A.adv: an AA male has an A.adv x higher chance of siring an offspring
# wm: relative weight of male abundance for offspring survival
# wf: relative weight of female abundance for offspring survival
# maxfreq: the maximum fraction of the population that is allowed to be on one island (hard density constrained) should be between 0.5 and 1.
#          when one island has more individuals than allowed by maxfreq, the difference automatically migrates to the other island.
# d: [0,1] dominance factor. The fitness of aan AA male is A.adv, of an Aa male A.adv^d and of an aa male 1.
#
# output:
# a length(popvect) x length(popvect) matrix that describes the population dynamics from t to t+1. Offspring is attributed to females only,
# but the distribution of the offspring depends on the male distribution on the islands.
make_mat <- function(surv_off,popvect,dum,A.adv,wm,wf,maxfreq=1,d,random_father=FALSE,rel_dens=TRUE){
  
  male.dists <- male.dist(dum,popvect,maxfreq,normalize = FALSE)
  female.dists <- female.dist(dum,popvect,normalize=FALSE) # TODO: combine dist functions 
  
  A <- matrix(0,nrow=18,ncol=18)
  
  #diag(A) <- surv
  
  N.f.1 <- sum(female.dists[[1]])
  N.f.2 <- sum(female.dists[[2]])

  N.m.1 <- sum(male.dists[[1]])
  N.m.2 <- sum(male.dists[[2]])
  
  # in order to generate patch specific survival
  # surv.1 and surv.2 vectors of length nrow(dum) that describe the survival probability of different genotypes
  # on the different beaches. Could be density dependent. Could be different for males and females
  diag(A) <- surv_adult(dum,male.dists,female.dists) #//dum$p1*surv.1 + (1-dum$p1)*surv.2
  
  if(rel_dens){
    surv.1 <- 1-surv_off(wm/(wm+wf)*N.m.1/(N.m.1+N.m.2) + wf/(wm+wf)*N.f.1/(N.f.1+N.f.2)) #1-dens_reg+dens_reg*plogis(5*(0.25-N.f.1))
    surv.2 <- 1-surv_off(wm/(wm+wf)*N.m.2/(N.m.1+N.m.2) + wf/(wm+wf)*N.f.2/(N.f.1+N.f.2)) #1-dens_reg+dens_reg*plogis(5*(0.25-N.f.2))
  }else{
    surv.1 <- 1-surv_off(wm/(wm+wf)*N.m.1 + wf/(wm+wf)*N.f.1) #1-dens_reg+dens_reg*plogis(5*(0.25-N.f.1))
    surv.2 <- 1-surv_off(wm/(wm+wf)*N.m.2 + wf/(wm+wf)*N.f.2) #1-dens_reg+dens_reg*plogis(5*(0.25-N.f.2))
    
  }
  A.adv.1 <- A.adv #2*plogis(5*(ml1))
  A.adv.2 <- A.adv #2*plogis(5*(ml2))
  
  male.vals.1 <- male.vals(dum,male.dists[[1]]/N.m.1,A.adv.1,d)
  male.vals.2 <- male.vals(dum,male.dists[[2]]/N.m.2,A.adv.2,d)
  
  if(random_father){
    male.vals.1 <- male.vals(dum,(male.dists[[1]]+male.dists[[2]])/(N.m.1+N.m.2),A.adv.1,d)
    male.vals.2 <- male.vals.1
  }
  
  for(i in which(dum$sex == 'f')){
    if(N.m.1 > 1e-32){
    offs.dist.1 <- calc_off_dist_alt(dum[i,],dum,male.vals.1)
    }else{
      offs.dist.1 <- rep(0,18)
    }
    
    if(N.m.2 > 1e-32){
    offs.dist.2 <- calc_off_dist_alt(dum[i,],dum,male.vals.2)
    }else{
      offs.dist.2 <- rep(0,18)
    }
    # offs.dist.1.fake <- calc_off_dist(dum[i,],dum,male.dist.1,A.adv)
    # offs.dist.2.fake <- calc_off_dist(dum[i,],dum,male.dist.2,1) # male advantage of having genotype AA only counts on island 1
    # 
    # cat(all.equal(offs.dist.1.fake,offs.dist.1),"\t",
    # all.equal(offs.dist.2.fake,offs.dist.2),"\n")

    # change these p1's to change female distribution across the islands
    A[,i] <- A[,i] + dum$p1[i]*offs.dist.1*surv.1 + (1-dum$p1[i])*offs.dist.2*surv.2
    
  }
  
  return(A)
}

# function that calculates euclidian distance.
# input:
# x,y: numeric vectors
#
# output:
# the euclidian distance between x and y
distance <- function(x,y){
  sqrt(sum((x-y)^2))
}

## Matrix-like model for the Fur Seals
# main function that runs the simulations.
#
# Input:
# filename: character, file where output is stored. If filename="bla", the code creates "bla.csv" and "bla.dum" [unless dumgen or test is TRUE]
# surv: survival of adults between two time steps. This acts to scale between offspring and adults. (b/c overall density is constant and 1)
# surv_off: function with 1 argument that returns 1 - offspring survival
# A.adv: an AA male has an A.adv x higher chance of siring an offspring
# wm: relative weight of male abundance for offspring survival
# wf: relative weight of female abundance for offspring survival
# Nt: number of time steps
# min_val_m: the probability that an aa male goes to island 2 and that an AA male goes to island 1
# min_val_f: the probability that a bb female goes to island 2 and that a BB female goes to island 1
# N0: starting population vector
# tol: minimum euclidian difference between N(t) and N(t+1) for a new data point to be stored
# maxfreq: the maximum fraction of the population that is allowed to be on one island (hard density constrained) should be between 0.5 and 1.
#          when one island has more individuals than allowed by maxfreq, the difference automatically migrates to the other island.
# d: [0,1] dominance factor. The fitness of aan AA male is A.adv, of an Aa male A.adv^d and of an aa male 1.
# dumgen: TRUE/FALSE if TRUE, only a dum object is generated and no actual simulations are run
# test: TRUE/FALSE if TRUE, no files are generated and instead both the dum and population trajectories are directly returned as a list
#
# Output:
# depending on 'test' and 'dumgen' one of the following:
# TRUE/FALSE a list of ts (timeseries) and the descriptive dum object returned directly
# FALSE/TRUE or TRUE/TRUE only the descriptive dum object
# FALSE/FALSE timeseries and descriptive dum object are stored to two files, filename.csv and filename.dum
#
run_sim <- function(filename,surv=0,surv_off=function(n) 0,A.adv=1.5,wm=1,wf=0,Nt=1e3, min_val_m=0.3, min_val_f=0.1,N0,tol=1e-4,maxfreq=1,d=0.5,d2=0.5,dumgen=FALSE,test=FALSE,Apenalty=0.1,random_father=FALSE,rel_dens=TRUE){
  if(missing(N0)){
    N0 <- runif(18)
  }
  
  # descriptive object that contains all possible genotypes and sexes
  dum <- expand.grid(a1 = c('a','A'),a2 = c('a','A'),b1 = c('b','B'),b2 = c('b','B'),sex = c('m','f'),stringsAsFactors = FALSE)
  dum <- dum[!(dum$a1 == 'A' & dum$a2 == 'a') & !(dum$b1 =='B' & dum$b2 == 'b'),]

  dum2 <- dum
  
  # translating genotypes into integer values (aa = 0, aA = 1, AA = 2, same for b)
  dum2$Ascore <- as.numeric(dum2$a1 == 'A') + as.numeric(dum2$a2 == 'A')
  dum2$Bscore <- as.numeric(dum2$b1 == 'B') + as.numeric(dum2$b2 == 'B')
  
  # adding the chances of going to island 1 to each class in the dum object:
  dum2$p1 <- as.numeric(dum2$sex=='m')*(min_val_m+(as.numeric(dum2$Ascore==1)*d2+as.numeric(dum2$Ascore==2))*(1-2*min_val_m)) + 
              as.numeric(dum2$sex == 'f')*(min_val_f+(as.numeric(dum2$Bscore==1)*d2+as.numeric(dum2$Bscore==2))*(1-2*min_val_f))
  dum2$surv <- (1 - Apenalty)^(as.numeric(dum2$sex=='m')*(as.numeric(dum2$Ascore==2) + d*as.numeric(dum2$Ascore==1)))
    
  if(dumgen){return(dum2)}
  if(!test){write.csv(dum2,file=paste(filename,".dum",sep=""))}

  # objects to store timeseries
  store <- matrix(NA,nrow=length(N0),ncol=Nt)
  
  # object to store timeseries, but first row is set to be time units
  final_store <- matrix(NA,nrow=1+length(N0),ncol=Nt)
  
  # normalization and initializatoin
  store[,1] <- N0 / sum(N0)
  final_store[,1] <- c(1,store[,1])
  max_store <- 1 # highest index currently used in the final_store matrix
  for(t in 2:Nt){
    A <- make_mat(surv_off,store[,t-1],dum2,A.adv,wm,wf,maxfreq = maxfreq,d,random_father = random_father,rel_dens=rel_dens)
    store[,t] <- A %*% store[,t-1] 
    store[,t] <- store[,t]/sum(store[,t]) # keeping population size constant

    if(distance(store[,t],final_store[-1,max_store]) > tol | t == Nt){ # only store new data point if difference with currently stored value exceeds tol
      max_store <- max_store + 1
      final_store[,max_store] <- c(t,store[,t])
    }
  }
  
  if(!test){
      write.csv(t(final_store[,1:max_store]),file=paste(filename,".csv",sep=""))
    }else{
      return(list(ts = final_store[,1:max_store],dum=dum2))
  }
}
