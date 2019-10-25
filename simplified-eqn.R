f <- function(N,Na,s){
  N/ ((1-s)*Na + N - Na)
}

N_granddaughters <- function(N,Na, s){
  Nb <- N - Na
  f1 <- f(N,Na,s)
  
  Na2 <- (1-s)*f1*Na 
  Nb2 <- f1*Nb 
  
  f2 <- f(Nb2+Na2,Na2,s)
  
  c((1-s)*f2*Na2/Na + N/Na,f2*Nb2/Nb)
  # (1-s)*f*Na*f*(1-s) + f*Nb*f
  
  #(N-Na)*((1-s)^2*f^2 +  N/Na) + f^2*(Na)
}

N_granddaughters(17,15,0.4)
