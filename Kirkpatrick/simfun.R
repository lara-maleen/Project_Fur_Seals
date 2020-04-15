z <- function(s,t2,a2){
  1 + (1-s)*(a2-1)*t2/(1-s*t2)
}

B <- function(s,t2,p2,a2){
  (1-s)/(1-s*t2)*(1+p2*(a2/z(s,t2,a2)-1))
}

delta_t2 <- function(s,t2,p2,a2){
  0.5*t2*(B(s,t2,p2,a2) - 1)
}

delta_p2 <- function(s,t2,p2,a2,D){
  delta_t2(s,t2,p2,a2)*D/(t2*(1-t2))  
}

D <- function(x1,x2,x3,x4){
  x1*x4-x2*x3  
}

equil <- function(s,p2,a2){
  (1/s+1/(a2*(1-s)-1))*p2 - 1/(a2*(1-s)-1)
}

equil_f <- function(s,p2,a2){
  (1/s+1/(a2*(1-s)-1))*p2/2 - 1/(a2*(1-s)-1)
}

equil_offs <- function(sfun,t2,a2){
  - sfun(t2)*(1 + t2*(a2 - 1))/((sfun(t2) -1 - sfun(t2)*t2)*(a2 -1))
}

t2 <- function(x1,x2,x3,x4){
  x3+x4  
}

p2 <- function(x1,x2,x3,x4){
  x2+x4
}

time_iter_offs_surv <- function(filename,x1,x2,x3,x4,sfun,a2,Nt=1,verb=FALSE,tol=0){
  x1mat <- matrix(c(
    1,   0.5,   0.5,  0.25,
    0.5, 0,     0.25, 0,
    0.5, 0.25,  0,    0,
    0.25,0,     0,    0
  ),ncol=4,byrow=TRUE)
  
  
  # T1P2
  x2mat <- matrix(c(
    0,   0.5,  0,     0.25,
    0.5, 1,    0.25,  0.5,
    0,   0.25, 0,     0,
    0.25,0.5,  0,     0
  ),ncol=4,byrow=TRUE)
  
  # T2P1
  x3mat <- matrix(c(
    0, 0, 0.5, 0.25,
    0, 0, 0.25,0,
    0.5,0.25,1,0.5,
    0.25,0,0.5,0
  ),ncol=4,byrow=TRUE)
  
  # T2P2
  x4mat <- matrix(c(
    0, 0, 0, 0.25,
    0, 0, 0.25,0.5,
    0, 0.25,0,0.5,
    0.25,0.5,0.5,1
  ),ncol=4,byrow=TRUE)
  
  out <- matrix(NA,ncol=Nt+1, nrow=4)
  tpoints <- numeric(Nt+1)
  tpoints[1] <- 1
  out[,1] <- c(x1,x2,x3,x4)
  
  
  s <- 0
  maxt <- 1
  for(t in 1:Nt){
    t2h <- t2(x1,x2,x3,x4)
    p2h <- p2(x1,x2,x3,x4)
    Dh <- D(x1,x2,x3,x4)
    
    x1p <- x1/(1-s*t2h)
    x2p <- x2/(1-s*t2h)
    x3p <- (1-s)*x3/(1-s*t2h)
    x4p <- (1-s)*x4/(1-s*t2h)
    
    x1pf <- x1/(1-s*t2h)
    x2pf <- (1-s)*x2/(1-s*t2h)
    x3pf <- x3/(1-s*t2h)
    x4pf <- (1-s)*x4/(1-s*t2h)
    
    zval <- z(s,t2h,a2)
    
    tab1 <- matrix(c(
      x1pf*x1p,        x1pf*x2p,       x1pf*x3p,          x1pf*x4p,
      x2pf*x1p/zval,   x2pf*x2p/zval,  a2*x2pf*x3p/zval,  a2*x2pf*x4p/zval,
      x3pf*x1p,        x3pf*x2p,       x3pf*x3p,          x3pf*x4p,
      x4pf*x1p/zval,   x4pf*x2p/zval,  a2*x4pf*x3p/zval,  a2*x4pf*x4p/zval
    ),ncol=4,byrow=TRUE)
    # cat("tab1-sum: ", sum(tab1),"\n")
    
    # applying survival penalty on the offspring
    # cat(sum(tab1),"\n")
    s_offs <- sfun(t2h)
    tab1[,3:4] <- (1-s_offs)*tab1[,3:4]
    # cat(sum(tab1),"\n")
    tab1 <- tab1/sum(tab1)
    
    x1n <- sum(tab1 * x1mat)
    x2n <- sum(tab1 * x2mat)
    x3n <- sum(tab1 * x3mat)
    x4n <- sum(tab1 * x4mat)
    
    if(verb){
      cat("Delta t2:", t2(x1n,x2n,x3n,x4n)-t2h,"\t Delta p2:",p2(x1n,x2n,x3n,x4n)- p2h,"\n")
    }
    tot <- x1n + x2n + x3n + x4n
    x1 <- x1n/tot
    x2 <- x2n/tot
    x3 <- x3n/tot
    x4 <- x4n/tot
    
    distance <- sqrt(sum((c(x1,x2,x3,x4)-out[,maxt])^2))
    if(t == Nt | distance >= tol){
      out[,maxt+1] <- c(x1,x2,x3,x4)
      tpoints[maxt+1] <- t
      maxt <- maxt + 1
    }
  }
  output <- data.frame(t = tpoints[1:maxt],x1 = out[1,1:maxt], x2 = out[2,1:maxt], x3 = out[3,1:maxt],x4 = out[4,1:maxt])
  write.csv(output,file = filename)
}