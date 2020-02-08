
dist <- function(x,y){
  sqrt(sum((x-y)^2))
}

statistics <- function(filename,surv,A.adv,empty=FALSE){
  if(empty){
    return(data.frame(stability=NA,N1=NA,N2=NA))
  }
  # stability (using perturbations?!)
  dum <- read.csv(paste(filename,".dum",sep=""),stringsAsFactors = FALSE)
  dat <- t(read.csv(paste(filename,".csv",sep=""),stringsAsFactors = FALSE))
  
  delta <- 0.01 #fractional perturbation
  base_vec <- as.numeric(dat[nrow(dat),])
  Tp <- 10
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
      A_up <- make_mat(surv,pert_up,dum,A.adv)
      pert_up <- A_up %*% pert_up
      pert_up <- pert_up/sum(pert_up)
      
      if(!atzero){
        A_down <- make_mat(surv,pert_down,dum,A.adv)
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
  
  # summarize these numbers in relevant statistics
  
  # differences in allele frequencies
  
  # write time state file, if requested
  write.csv(sumdat,filename = paste(filename,".stat",sep=""))
  
  # return last data point stats
  return(data.frame(stability = pert_ind, N1 = sumdat$N.1[nrow(sumdat)], N2 = sumdat$N.2[nrow(sumdat)]))
}