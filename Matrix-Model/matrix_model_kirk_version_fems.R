calc_off_dist_alt <- function(fem,dum,male.vals){
  
  as.f <- 0.5*(as.numeric(fem$a1 == 'A') + as.numeric(fem$a2 == 'A')) # prob. of inheriting A from mother
  bs.f <- 0.5*(as.numeric(fem$b1 == 'B') + as.numeric(fem$b2 == 'B')) # prob. of inheriting B from mother
  

  # creating a 3x3 matrix with all potential combinations of scores
  # Ascore 0,1,2 and Bscore 0,1,2
  
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
  offs <- offs.score[cbind(dum$Bscore + 1,dum$Ascore + 1)] 
  
  offs/2
}


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



male.dist <- function(dum,popvect,maxfreq,normalize=TRUE){
  # male dist isle 1
  male.dist.1 <- dum$p1*popvect*as.numeric(dum$sex=='m')

  # male dist isle 2
  male.dist.2 <- (1-dum$p1)*popvect*as.numeric(dum$sex=='m')
  
  freq <- sum(male.dist.1) / sum(male.dist.1 + male.dist.2)

  if(freq > maxfreq){
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

male.vals <- function(dum,male.dist,A.adv){
  
  male.dist[dum$Ascore == 2] <- A.adv*male.dist[dum$Ascore == 2]
  male.dist <- male.dist/sum(male.dist)
  
  mAB <- sum(male.dist*0.5*dum$Ascore*0.5*dum$Bscore) # prob. of inheriting AB from father
  mAb <- sum(male.dist*0.5*dum$Ascore*(1-0.5*dum$Bscore)) # prob. of inheriting Ab from father
  maB <- sum(male.dist*(1-0.5*dum$Ascore)*0.5*dum$Bscore) # prob. of inheriting aB from father
  mab <- sum(male.dist*(1-0.5*dum$Ascore)*(1-0.5*dum$Bscore)) # prob. of inheriting ab from father
  
  list(mAB=mAB,mAb=mAb,maB=maB,mab=mab)
}

make_mat <- function(surv,surv_off,popvect,dum,A.adv,maxfreq=1){
  
  male.dists <- male.dist(dum,popvect,maxfreq,normalize = FALSE)
  
  A <- matrix(0,nrow=18,ncol=18)
  diag(A) <- surv
  N.f.1 <- sum((dum$p1*popvect)[dum$sex == 'f'])
  N.f.2 <- sum(((1-dum$p1)*popvect)[dum$sex == 'f'])
  N.m.1 <- sum(male.dists[[1]])
  N.m.2 <- sum(male.dists[[2]])
  
  surv.1 <- 1-surv_off(N.f.1/(N.f.1+N.f.2)) #1-dens_reg+dens_reg*plogis(5*(0.25-N.f.1))
  surv.2 <- 1-surv_off(N.f.2/(N.f.1+N.f.2)) #1-dens_reg+dens_reg*plogis(5*(0.25-N.f.2))

  A.adv.1 <- A.adv #2*plogis(5*(ml1))
  A.adv.2 <- A.adv #2*plogis(5*(ml2))
  
  male.vals.1 <- male.vals(dum,male.dists[[1]]/N.m.1,A.adv.1)
  male.vals.2 <- male.vals(dum,male.dists[[2]]/N.m.2,A.adv.2)
  
  for(i in which(dum$sex == 'f')){
    if(N.m.1 > 1e-32){
    offs.dist.1 <- calc_off_dist_alt(dum[i,],dum,male.vals.1)
    }else{
      offs.dist.1 <- rep(0,18)
    }
    
    if(N.m.2 > 1e-32){
    offs.dist.2 <- calc_off_dist_alt(dum[i,],dum,male.vals.2) # male advantage of having genotype AA only counts on island 1
    }else{
      offs.dist.2 <- rep(0,18)
    }
    # offs.dist.1.fake <- calc_off_dist(dum[i,],dum,male.dist.1,A.adv)
    # offs.dist.2.fake <- calc_off_dist(dum[i,],dum,male.dist.2,1) # male advantage of having genotype AA only counts on island 1
    # 
    # cat(all.equal(offs.dist.1.fake,offs.dist.1),"\t",
    # all.equal(offs.dist.2.fake,offs.dist.2),"\n")

    
    A[,i] <- A[,i] + dum$p1[i]*offs.dist.1*surv.1 + (1-dum$p1[i])*offs.dist.2*surv.2
    
  }
  
  return(A)
}

distance <- function(x,y){
  sqrt(sum((x-y)^2))
}

## Matrix-like model for the Fur Seals
run_sim <- function(filename,surv=0,surv_off=function(n) 1,A.adv=1.5,Nt=1e3, min_val_m=0.3, min_val_f=0.1,N0,tol=1e-4,maxfreq=1){
  if(missing(N0)){
    N0 <- runif(18)
  }
  
  dum <- expand.grid(a1 = c('a','A'),a2 = c('a','A'),b1 = c('b','B'),b2 = c('b','B'),sex = c('m','f'),stringsAsFactors = FALSE)
  dum <- dum[!(dum$a1 == 'A' & dum$a2 == 'a') & !(dum$b1 =='B' & dum$b2 == 'b'),]

  dum2 <- dum
  dum2$p1 <- min_val_m*as.numeric(dum2$sex=='m') + min_val_f*as.numeric(dum2$sex == 'f') + 
    0.5*(1-2*min_val_m)*(as.numeric(dum2$sex=='m')*(as.numeric(dum2$a1=='A') + as.numeric(dum2$a2=='A'))) + 0.5*(1-2*min_val_f)*(as.numeric(dum2$sex=='f')*(as.numeric(dum2$b1=='B') + as.numeric(dum2$b2=='B'))) # probability of going to island 1 for each indiv
  dum2$Ascore <- as.numeric(dum2$a1 == 'A') + as.numeric(dum2$a2 == 'A')
  dum2$Bscore <- as.numeric(dum2$b1 == 'B') + as.numeric(dum2$b2 == 'B')
  
  write.csv(dum2,file=paste(filename,".dum",sep=""))

  store <- matrix(NA,nrow=length(N0),ncol=Nt)
  final_store <- matrix(NA,nrow=1+length(N0),ncol=Nt)
  
  store[,1] <- N0 / sum(N0)
  final_store[,1] <- c(1,store[,1])
  max_store <- 1 # highest index currently used in the final_store matrix
  for(t in 2:Nt){
    A <- make_mat(surv,surv_off,store[,t-1],dum2,A.adv,maxfreq = maxfreq)
    store[,t] <- A %*% store[,t-1] 
    store[,t] <- store[,t]/sum(store[,t])

    if(distance(store[,t],final_store[-1,max_store]) > tol | t == Nt){
      max_store <- max_store + 1
      final_store[,max_store] <- c(t,store[,t])
    }
  }
  
  write.csv(t(final_store[,1:max_store]),file=paste(filename,".csv",sep=""))
}
