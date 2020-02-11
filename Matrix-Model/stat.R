
dist <- function(x,y){
  sqrt(sum((x-y)^2))
}

statistics <- function(filename,surv,A.adv,dens_reg,Tp=10,empty=FALSE){
  if(empty){
    return(data.frame(stability=NA,N1=NA,N2=NA,N.1.m=NA,N.2.m=NA,N.1.f=NA,N.2.f=NA,
                      N.1.AA.m = NA, N.1.Aa.m = NA, N.1.aa.m = NA,
                      N.2.AA.m = NA, N.2.Aa.m = NA, N.2.aa.m = NA,
                      N.1.BB.f = NA, N.1.Bb.f = NA, N.1.bb.f = NA,
                      N.2.BB.f = NA, N.2.Bb.f = NA, N.2.bb.f = NA))
  }
  # stability (using perturbations?!)
  dum <- read.csv(paste(filename,".dum",sep=""),stringsAsFactors = FALSE)
  dat <- t(read.csv(paste(filename,".csv",sep=""),stringsAsFactors = FALSE))
  
  delta <- 0.01 #fractional perturbation
  base_vec <- as.numeric(dat[nrow(dat),])

  up_pert <- logical(length(base_vec))
  down_pert <- logical(length(base_vec))
  pert_ind <- 0
  max_dist <- 0
  for(i in 1:length(base_vec)){
    atzero <- round(base_vec[i],4)==0
    pert_up <- base_vec
    pert_down <- base_vec
    
    pert_up[i] <- pert_up[i] + delta
    pert_down[i] <- pert_down[i] - delta
    
    pert_down[pert_down < 0] <- 0
    
    pert_up <- pert_up/sum(pert_up)
    pert_down <- pert_down/sum(pert_down)
    
    dist_prior_up <- dist(pert_up,base_vec)
    dist_prior_down <- dist(pert_down,base_vec)
    
    for(t in 1:Tp){
      A_up <- make_mat(surv,pert_up,dum,A.adv,dens_reg=dens_reg)
      pert_up <- A_up %*% pert_up
      pert_up <- pert_up/sum(pert_up)
      
      if(!atzero){
        A_down <- make_mat(surv,pert_down,dum,A.adv,dens_reg=dens_reg)
        pert_down <- A_down %*% pert_down
        pert_down <- pert_down/sum(pert_down)
      }
      
      # if(i == 8){
      #   pert_up <- 0.9*pert_up + 0.1*rep(1/length(base_vec),length(base_vec))
      # }
      # 
      # if(i == 6){
      #   pert_down <- 0.0*pert_down + 1.0*rep(1/length(base_vec),length(base_vec))
      # }
      
    }
    
    dist_post_up <- dist(pert_up,base_vec)
    dist_post_down <- dist(pert_down,base_vec)
    
    # cat(dist_post_up,"\t",dist_prior_up,"\n")
    
    up_pert <- dist_post_up - dist_prior_up # negative values indicate that the distance has decreased over time
    # and that the system is stable.
    # positive values indicate an system that's still moving forward.
    # more positive values mean larger differences over time
    down_pert <- dist_post_down - dist_prior_down
    
    
    if(up_pert > 0 & up_pert > max_dist){
      pert_ind <- i
      max_dist <- up_pert
    }
    
    if(!atzero & down_pert > 0 & down_pert > max_dist){
      pert_ind <- -i
      max_dist <- down_pert
    }
    
  }
  
  
  # differences in abundances between the islands
  
  sumdat <- dum
  
  # add columns for island 1 and 2
  sumdat$N.1 <- sumdat$p1*base_vec
  sumdat$N.2 <- (1-sumdat$p1)*base_vec
  
  N1 <- sum(sumdat$N.1)
  N2 <- sum(sumdat$N.2)
  
  N.1.m <- sum(sumdat$N.1[sumdat$sex=='m'])
  N.2.m <- sum(sumdat$N.2[sumdat$sex=='m'])
  N.1.f <- sum(sumdat$N.1[sumdat$sex=='f'])
  N.2.f <- sum(sumdat$N.2[sumdat$sex=='f'])
  
  N.1.AA.m <- sum(sumdat$N.1[sumdat$sex=='m' & sumdat$Ascore == 2])
  N.1.Aa.m <- sum(sumdat$N.1[sumdat$sex=='m' & sumdat$Ascore == 1])
  N.1.aa.m <- sum(sumdat$N.1[sumdat$sex=='m' & sumdat$Ascore == 0])
    
  N.2.AA.m <- sum(sumdat$N.2[sumdat$sex=='m' & sumdat$Ascore == 2])
  N.2.Aa.m <- sum(sumdat$N.2[sumdat$sex=='m' & sumdat$Ascore == 1])
  N.2.aa.m <- sum(sumdat$N.2[sumdat$sex=='m' & sumdat$Ascore == 0])
  
  N.1.BB.f <- sum(sumdat$N.1[sumdat$sex=='f' & sumdat$Bscore == 2])
  N.1.Bb.f <- sum(sumdat$N.1[sumdat$sex=='f' & sumdat$Bscore == 1])
  N.1.bb.f <- sum(sumdat$N.1[sumdat$sex=='f' & sumdat$Bscore == 0])
  
  N.2.BB.f <- sum(sumdat$N.2[sumdat$sex=='f' & sumdat$Bscore == 2])
  N.2.Bb.f <- sum(sumdat$N.2[sumdat$sex=='f' & sumdat$Bscore == 1])
  N.2.bb.f <- sum(sumdat$N.2[sumdat$sex=='f' & sumdat$Bscore == 0])
  
  # summarize these numbers in relevant statistics
  
  # differences in allele frequencies
  
  # write time state file, if requested
  # write.csv(sumdat,file = paste(filename,".stat",sep=""))
  
  # return last data point stats
  return(data.frame(stability = pert_ind, N1 = N1, N2 = N2, N.1.m = N.1.m, N.2.m = N.2.m, N.1.f = N.1.f, N.2.f = N.2.f,
                    N.1.AA.m = N.1.AA.m, N.1.Aa.m = N.1.Aa.m, N.1.aa.m = N.1.aa.m,
                    N.2.AA.m = N.2.AA.m, N.2.Aa.m = N.2.Aa.m, N.2.aa.m = N.2.aa.m,
                    N.1.BB.f = N.1.BB.f, N.1.Bb.f = N.1.Bb.f, N.1.bb.f = N.1.bb.f,
                    N.2.BB.f = N.2.BB.f, N.2.Bb.f = N.2.Bb.f, N.2.bb.f = N.2.bb.f))
}