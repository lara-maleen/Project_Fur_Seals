statistics <- function(filename,empty=FALSE){
  if(empty){
    return(data.frame(N1=NA,N2=NA,N.1.m=NA,N.2.m=NA,N.1.f=NA,N.2.f=NA,
                      N.1.AA.m = NA, N.1.Aa.m = NA, N.1.aa.m = NA,
                      N.2.AA.m = NA, N.2.Aa.m = NA, N.2.aa.m = NA,
                      N.1.BB.f = NA, N.1.Bb.f = NA, N.1.bb.f = NA,
                      N.2.BB.f = NA, N.2.Bb.f = NA, N.2.bb.f = NA))
  }

  dum <- read.csv(paste(filename,".dum",sep=""),stringsAsFactors = FALSE)
  dat <- read.csv(paste(filename,".csv",sep=""),stringsAsFactors = FALSE)
  
  # differences in abundances between the islands
  
  sumdat <- dum
  popvec <- as.numeric(dat[nrow(dat),-2:-1])
  
  # add columns for island 1 and 2
  sumdat$N.1 <- sumdat$p1*popvec
  sumdat$N.2 <- (1-sumdat$p1)*popvec
  
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
  return(data.frame(N1 = N1, N2 = N2, N.1.m = N.1.m, N.2.m = N.2.m, N.1.f = N.1.f, N.2.f = N.2.f,
                    N.1.AA.m = N.1.AA.m, N.1.Aa.m = N.1.Aa.m, N.1.aa.m = N.1.aa.m,
                    N.2.AA.m = N.2.AA.m, N.2.Aa.m = N.2.Aa.m, N.2.aa.m = N.2.aa.m,
                    N.1.BB.f = N.1.BB.f, N.1.Bb.f = N.1.Bb.f, N.1.bb.f = N.1.bb.f,
                    N.2.BB.f = N.2.BB.f, N.2.Bb.f = N.2.Bb.f, N.2.bb.f = N.2.bb.f))
}

get_stats <- function(dir){
  all_files <- read.csv(file.path(dir,"simruns.csv"))
  
  new_inf <- lapply(all_files$outfile, FUN = function(x) statistics(file.path(dir,"raw",x)))
  
  new_inf2 <- do.call(rbind,new_inf)
  
  write.csv(new_inf2,file.path(dir,"simruns.stat"))
}