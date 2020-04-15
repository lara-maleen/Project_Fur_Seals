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

time_iter <- function(x1,x2,x3,x4,s,a2,Nt=1,verb=FALSE){
  
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
  out[,1] <- c(x1,x2,x3,x4)
  
  for(t in 1:Nt){
    t2h <- t2(x1,x2,x3,x4)
    p2h <- p2(x1,x2,x3,x4)
    Dh <- D(x1,x2,x3,x4)
    
    x1p <- x1/(1-s*t2h)
    x2p <- x2/(1-s*t2h)
    x3p <- (1-s)*x3/(1-s*t2h)
    x4p <- (1-s)*x4/(1-s*t2h)
    
    zval <- z(s,t2h,a2)
    
    tab1 <- matrix(c(
      x1*x1p,     x1*x2p,    x1*x3p,       x1*x4p,
      x2*x1p/zval,   x2*x2p/zval,  a2*x2*x3p/zval,  a2*x2*x4p/zval,
      x3*x1p,     x3*x2p,    x3*x3p,       x3*x4p,
      x4*x1p/zval,   x4*x2p/zval,  a2*x4*x3p/zval,  a2*x4*x4p/zval
    ),ncol=4,byrow=TRUE)
    # cat("tab1-sum: ", sum(tab1),"\n")
  
    x1n <- sum(tab1 * x1mat)
    x2n <- sum(tab1 * x2mat)
    x3n <- sum(tab1 * x3mat)
    x4n <- sum(tab1 * x4mat)
  
    if(verb){
      cat("Delta t2:", t2(x1n,x2n,x3n,x4n)-t2h," (",delta_t2(s,t2h,p2h,a2),")\t Delta p2",p2(x1n,x2n,x3n,x4n)- p2h,"(",delta_p2(s,t2h,p2h,a2,Dh),")\n")
    }
    tot <- x1n + x2n + x3n + x4n
    x1 <- x1n/tot
    x2 <- x2n/tot
    x3 <- x3n/tot
    x4 <- x4n/tot
    
    out[,t+1] <- c(x1,x2,x3,x4)
  }
  
  out
}

time_iter_fem_surv <- function(x1,x2,x3,x4,s,a2,Nt=1,verb=FALSE){
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
  out[,1] <- c(x1,x2,x3,x4)
  
  for(t in 1:Nt){
    t2h <- t2(x1,x2,x3,x4)
    p2h <- p2(x1,x2,x3,x4)
    Dh <- D(x1,x2,x3,x4)
  
    x1p <- x1/(1-s*t2h)
    x2p <- x2/(1-s*t2h)
    x3p <- (1-s)*x3/(1-s*t2h)
    x4p <- (1-s)*x4/(1-s*t2h)
  
    zval <- z(s,t2h,a2)
  
    tab1 <- matrix(c(
      x1p*x1p,        x1p*x2p,       x1p*x3p,          x1p*x4p,
      x2p*x1p/zval,   x2p*x2p/zval,  a2*x2p*x3p/zval,  a2*x2p*x4p/zval,
      x3p*x1p,        x3p*x2p,       x3p*x3p,          x3p*x4p,
      x4p*x1p/zval,   x4p*x2p/zval,  a2*x4p*x3p/zval,  a2*x4p*x4p/zval
    ),ncol=4,byrow=TRUE)
  # cat("tab1-sum: ", sum(tab1),"\n")
  
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
    
    out[,t+1] <- c(x1,x2,x3,x4)
  }
  out  
}

time_iter_fem_surv_p <- function(x1,x2,x3,x4,s,a2,Nt=1,verb=FALSE){
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
  out[,1] <- c(x1,x2,x3,x4)
  
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
    
    out[,t+1] <- c(x1,x2,x3,x4)
  }
  out  
}



time_iter_offs_surv <- function(x1,x2,x3,x4,sfun,a2,Nt=1,verb=FALSE,tol=0){
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
  rbind(tpoints[1:maxt],out[,1:maxt])
}
par(mfrow=c(2,2))

# testing theoretical functions function(s,t2,p2,a2)
# delta_t2(0.3,t2(x1,x2,x3,x4),p2(x1,x2,x3,x4),2)
# delta_p2(0.3,t2(x1,x2,x3,x4),p2(x1,x2,x3,x4),2,D(x1,x2,x3,x4))
a2 <- 2.5
s <- 0.25
Nt <- 1e3
Nrep <- 50
cols <- rainbow(Nrep)

p2s <- seq(0,1,0.01)
t2s <- equil(s,p2s,a2)
t2s[t2s < 0 ] <- 0
t2s[t2s>1] <- 1
plot(p2s,t2s,xlim=c(0,1),ylim=c(0,1),xlab="p2",ylab="t2",type="l",main="T affects males")

for(i in 1:Nrep){
nums <- runif(4)
nums <- nums/sum(nums)
x1 <- nums[1]
x2 <- nums[2]
x3 <- nums[3]
x4 <- nums[4]
tmp <- time_iter(x1,x2,x3,x4,s,a2,Nt=Nt)
p2n <-apply(tmp[c(2,4),],2,sum)
t2n <- apply(tmp[c(3,4),],2,sum)
lines(p2n,t2n,col=cols[i])
points(p2n[length(p2n)],t2n[length(t2n)],cex=2,col=cols[i])
}

########### T affects males and females
p2s <- seq(0,1,0.01)
t2s <- equil_f(s,p2s,a2)
t2s[t2s < 0 ] <- 0
t2s[t2s>1] <- 1
plot(p2s,t2s,xlim=c(0,1),ylim=c(0,1),xlab="p2",ylab="t2",type="l",main="T affects males and females")

for(i in 1:Nrep){
  nums <- runif(4)
  nums <- nums/sum(nums)
  x1 <- nums[1]
  x2 <- nums[2]
  x3 <- nums[3]
  x4 <- nums[4]
  tmp <- time_iter_fem_surv(x1,x2,x3,x4,s,a2,Nt=Nt)
  p2n <-apply(tmp[c(2,4),],2,sum)
  t2n <- apply(tmp[c(3,4),],2,sum)
  lines(p2n,t2n,col=cols[i])
  points(p2n[length(p2n)],t2n[length(t2n)],cex=2,col=cols[i])
}


########### T affects males and P affects females

p2s <- seq(0,1,0.01)
t2s <- equil_f(s,p2s,a2)
t2s[t2s < 0 ] <- 0
t2s[t2s>1] <- 1
plot(p2s,t2s,xlim=c(0,1),ylim=c(0,1),xlab="p2",ylab="t2",type="l",main="T affects males and P affects females")

for(i in 1:Nrep){
  nums <- runif(4)
  nums <- nums/sum(nums)
  x1 <- nums[1]
  x2 <- nums[2]
  x3 <- nums[3]
  x4 <- nums[4]
  tmp <- time_iter_fem_surv_p(x1,x2,x3,x4,s,a2,Nt=Nt)
  p2n <-apply(tmp[c(2,4),],2,sum)
  t2n <- apply(tmp[c(3,4),],2,sum)
  lines(p2n,t2n,col=cols[i])
  points(p2n[length(p2n)],t2n[length(t2n)],cex=2,col=cols[i])
}


########### T affects offspring survival
# par(mfrow=c(1,1))
# s <- 0.22
# Nrep <- 1
# Nt <- 1e3
# a2 <- 6000
sfuns <- list(constant = function(t2,sval,sslope) sval,
            linear = function(t2,sval,sslope) sval*t2,
            logistic.max.s = function(t2,sval,sslope) sval*plogis(sslope*(t2-0.5)))
matplot(seq(0,1,0.01),do.call(cbind,lapply(sfuns,function(x) x(seq(0,1,0.01),0.5,20))),type="l",xlab="t2",ylab="survival penalty")

sfuns[['logistic']](20)
a2s <- c(1.2,1.5,3)
svals <- c(0.1,0.3,0.5)
sslopes <- c(7.5,20)
start <- expand.grid(x1 = seq(0,1,0.25), x2 = seq(0,1,0.25),x3 = seq(0,1,0.25))
start$x4 <- 1 - rowSums(start)
start <- start[!apply(start,1,FUN = function(x) any(x <0 | x > 1)),]
Nrep <- nrow(start)
tmp1 <- expand.grid(a2=a2s,stype= c('constant','linear'),sval = svals,sslope=NA,replicate=1:Nrep)
tmp2 <- expand.grid(a2=a2s,stype=c('logistic'),sval=svals,sslope=sslopes,replicate=1:Nrep)
overview <- rbind(tmp1,tmp2)

overview$filename <- paste("file",1:nrow(overview),".csv",sep="")
Nt <- 1e6
sfun <- function(t2){
  0.5*plogis(10*(t2-0.8))
}
a2 <- 1.2
t2s <- seq(0,1,0.01)
p2s <- equil_offs(sfun,t2s,a2)
p2s[p2s < 0 ] <- 0
p2s[p2s>1] <- 1
plot(p2s,t2s,xlim=c(0,1),ylim=c(0,1),xlab="p2",ylab="t2",type="l",main="T affects offspring")

Nt <- 1e6
for(i in 1:Nrep){
  nums <- runif(4)
  nums <- nums/sum(nums)
  x1 <- nums[1]
  x2 <- nums[2]
  x3 <- nums[3]
  x4 <- nums[4]
  tmp <- time_iter_offs_surv(x1,x2,x3,x4,sfun,a2,Nt=Nt,tol=0.01)
  p2n <-apply(tmp[c(3,5),],2,sum)
  t2n <- apply(tmp[c(4,5),],2,sum)
  lines(p2n,t2n,col=cols[i])
  points(p2n[length(p2n)],t2n[length(t2n)],cex=2,col=cols[i])
}

# p2n <- 0.55
# 
# eqnorm <- equil(0.3,p2n,2)
# x4 <- runif(1,eqnorm+p2n-1,min(c(p2n,eqnorm)))
# x2 <- p2n - x4
# x3 <- eqnorm - x4
# x1 <- 1 - x2 - x3 - x4
# time_iter(x1,x2,x3,x4,0.3,2)
# 
# p2n <- 0.9
# eqfem <- equil_f(0.3,p2n,2)
# x4 <- runif(1,eqfem+p2n-1,min(c(p2n,eqfem)))
# x2 <- p2n - x4
# x3 <- eqfem - x4
# x1 <- 1- x2 - x3 - x4
# # x1 <- 1 - p2n -eqfem + x4
# time_iter_fem_surv(x1,x2,x3,x4,0.3,2)
# 
