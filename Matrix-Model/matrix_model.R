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

male.vals <- function(dum,male.dist,A.adv){
  
  male.dist[dum$Ascore == 2] <- A.adv*male.dist[dum$Ascore == 2]
  male.dist <- male.dist/sum(male.dist)
  
  mAB <- sum(male.dist*0.5*dum$Ascore*0.5*dum$Bscore) # prob. of inheriting AB from father
  mAb <- sum(male.dist*0.5*dum$Ascore*(1-0.5*dum$Bscore)) # prob. of inheriting Ab from father
  maB <- sum(male.dist*(1-0.5*dum$Ascore)*0.5*dum$Bscore) # prob. of inheriting aB from father
  mab <- sum(male.dist*(1-0.5*dum$Ascore)*(1-0.5*dum$Bscore)) # prob. of inheriting ab from father
  
  list(mAB=mAB,mAb=mAb,maB=maB,mab=mab)
}

make_mat <- function(surv,popvect,dum,A.adv,dens_reg){
  
  # male dist isle 1
  male.dist.1 <- dum$p1*popvect*as.numeric(dum$sex=='m')
  if(sum(male.dist.1) > 0){
    male.dist.1 <- male.dist.1 / sum(male.dist.1)
  }
  # male dist isle 2
  male.dist.2 <- (1-dum$p1)*popvect*as.numeric(dum$sex=='m')
  if(sum(male.dist.2) > 0){
    male.dist.2 <- male.dist.2 / sum(male.dist.2)
  }
  
  A <- matrix(0,nrow=18,ncol=18)
  diag(A) <- surv
  N.f.1 <- sum((dum$p1*popvect)[dum$sex == 'f'])
  N.f.2 <- sum(((1-dum$p1)*popvect)[dum$sex == 'f'])
  surv.1 <- 1-dens_reg+dens_reg*plogis(5*(0.25-N.f.1))
  surv.2 <- 1-dens_reg+dens_reg*plogis(5*(0.25-N.f.2))
  
  male.vals.1 <- male.vals(dum,male.dist.1,A.adv)
  male.vals.2 <- male.vals(dum,male.dist.2,1)
  
  for(i in which(dum$sex == 'f')){
    offs.dist.1 <- calc_off_dist_alt(dum[i,],dum,male.vals.1)
    offs.dist.2 <- calc_off_dist_alt(dum[i,],dum,male.vals.2) # male advantage of having genotype AA only counts on island 1
    
    # offs.dist.1.fake <- calc_off_dist(dum[i,],dum,male.dist.1,A.adv)
    # offs.dist.2.fake <- calc_off_dist(dum[i,],dum,male.dist.2,1) # male advantage of having genotype AA only counts on island 1
    # 
    # cat(all.equal(offs.dist.1.fake,offs.dist.1),"\t",
    # all.equal(offs.dist.2.fake,offs.dist.2),"\n")

    
    A[,i] <- A[,i] + dum$p1[i]*offs.dist.1*surv.1 + (1-dum$p1[i])*offs.dist.2*surv.2
    
  }
  
  return(A)
}


## Matrix-like model for the Fur Seals
run_sim <- function(filename,surv=0,A.adv=1.5,Nt=1e3, min_val_m=0.3, min_val_f=0.1,dens_reg=0,N0){
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
  store[,1] <- N0 / sum(N0)
  for(t in 2:Nt){
    A <- make_mat(surv,store[,t-1],dum2,A.adv,dens_reg)
    store[,t] <- A %*% store[,t-1] 
    store[,t] <- store[,t]/sum(store[,t])
  }
  
  write.csv(store,file=paste(filename,".csv",sep=""))
}