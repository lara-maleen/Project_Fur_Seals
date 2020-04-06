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

t2 <- function(x1,x2,x3,x4){
  x3+x4  
}

p2 <- function(x1,x2,x3,x4){
  x2+x4
}

time_iter <- function(x1,x2,x3,x4,s,a2){
  t2h <- t2(x1,x2,x3,x4)
  p2h <- p2(x1,x2,x3,x4)
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
  cat("tab1-sum: ", sum(tab1),"\n")

  x1mat <- matrix(c(
    1,   0.5,   0.5,  0.25,
    0.5, 0,     0.25, 0,
    0.5, 0.25,  0,    0,
    0.25,0,     0,    0
  ),ncol=4,byrow=TRUE)
  
  x1n <- sum(tab1 * x1mat)
  
  # T1P2
  x2mat <- matrix(c(
    0,   0.5,  0,     0.25,
    0.5, 1,    0.25,  0.5,
    0,   0.25, 0,     0,
    0.25,0.5,  0,     0
  ),ncol=4,byrow=TRUE)
  
  x2n <- sum(tab1 * x2mat)
  
  # T2P1
  x3mat <- matrix(c(
    0, 0, 0.5, 0.25,
    0, 0, 0.25,0,
    0.5,0.25,1,0.5,
    0.25,0,0.5,0
  ),ncol=4,byrow=TRUE)
    
  x3n <- sum(tab1 * x3mat)

  # T2P2
  x4mat <- matrix(c(
    0, 0, 0, 0.25,
    0, 0, 0.25,0.5,
    0, 0.25,0,0.5,
    0.25,0.5,0.5,1
  ),ncol=4,byrow=TRUE)
  
  x4n <- sum(tab1 * x4mat)
  Dh <- D(x1,x2,x3,x4)
  cat("Delta t2:", t2(x1n,x2n,x3n,x4n)-t2h," (",delta_t2(s,t2h,p2h,a2),")\t Delta p2",p2(x1n,x2n,x3n,x4n)- p2h,"(",delta_p2(s,t2h,p2h,a2,Dh),")\n")
}

time_iter_fem_surv <- function(x1,x2,x3,x4,s,a2){
  t2h <- t2(x1,x2,x3,x4)
  p2h <- p2(x1,x2,x3,x4)
  
  x1p <- x1/(1-s*t2h)
  x2p <- x2/(1-s*t2h)
  x3p <- (1-s)*x3/(1-s*t2h)
  x4p <- (1-s)*x4/(1-s*t2h)
  
  zval <- z(s,t2h,a2)
  
  tab1 <- matrix(c(
    x1p*x1p,     x1p*x2p,    x1p*x3p,       x1p*x4p,
    x2p*x1p/zval,   x2p*x2p/zval,  a2*x2p*x3p/zval,  a2*x2p*x4p/zval,
    x3p*x1p,     x3p*x2p,    x3p*x3p,       x3p*x4p,
    x4p*x1p/zval,   x4p*x2p/zval,  a2*x4p*x3p/zval,  a2*x4p*x4p/zval
  ),ncol=4,byrow=TRUE)
  cat("tab1-sum: ", sum(tab1),"\n")
  
  x1mat <- matrix(c(
    1,   0.5,   0.5,  0.25,
    0.5, 0,     0.25, 0,
    0.5, 0.25,  0,    0,
    0.25,0,     0,    0
  ),ncol=4,byrow=TRUE)
  
  x1n <- sum(tab1 * x1mat)
  
  # T1P2
  x2mat <- matrix(c(
    0,   0.5,  0,     0.25,
    0.5, 1,    0.25,  0.5,
    0,   0.25, 0,     0,
    0.25,0.5,  0,     0
  ),ncol=4,byrow=TRUE)
  
  x2n <- sum(tab1 * x2mat)
  
  # T2P1
  x3mat <- matrix(c(
    0, 0, 0.5, 0.25,
    0, 0, 0.25,0,
    0.5,0.25,1,0.5,
    0.25,0,0.5,0
  ),ncol=4,byrow=TRUE)
  
  x3n <- sum(tab1 * x3mat)
  
  # T2P2
  x4mat <- matrix(c(
    0, 0, 0, 0.25,
    0, 0, 0.25,0.5,
    0, 0.25,0,0.5,
    0.25,0.5,0.5,1
  ),ncol=4,byrow=TRUE)
  
  x4n <- sum(tab1 * x4mat)
  Dh <- D(x1,x2,x3,x4)
  cat("Delta t2:", t2(x1n,x2n,x3n,x4n)-t2h,"\t Delta p2:",p2(x1n,x2n,x3n,x4n)- p2h,"\n")
}


# testing theoretical functions function(s,t2,p2,a2)
# delta_t2(0.3,t2(x1,x2,x3,x4),p2(x1,x2,x3,x4),2)
# delta_p2(0.3,t2(x1,x2,x3,x4),p2(x1,x2,x3,x4),2,D(x1,x2,x3,x4))

p2n <- 0.55

eqnorm <- equil(0.3,p2n,2)
x4 <- runif(1,eqnorm+p2n-1,min(c(p2n,eqnorm)))
x2 <- p2n - x4
x3 <- eqnorm - x4
x1 <- 1 - x2 - x3 - x4
time_iter(x1,x2,x3,x4,0.3,2)

p2n <- 0.9
eqfem <- equil_f(0.3,p2n,2)
x4 <- runif(1,eqfem+p2n-1,min(c(p2n,eqfem)))
x2 <- p2n - x4
x3 <- eqfem - x4
x1 <- 1- x2 - x3 - x4
# x1 <- 1 - p2n -eqfem + x4
time_iter_fem_surv(x1,x2,x3,x4,0.3,2)
