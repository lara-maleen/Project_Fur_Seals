# an expression for the fitness of AA, Aa and aa males based on female genotypic distribution

fitness.males <- function(popvect,dum,mvm=0,mvf=0,A2=1.8,d=0.5,d2=0.5,wf=0.5,surv=function(n) 0){
  
  N.m.1 <- sum(popvect*as.numeric(dum$sex=='m')*dum$p1)
  N.m.2 <- sum(popvect*as.numeric(dum$sex=='m')*dum$p1)
  N.f.1 <- sum(popvect*as.numeric(dum$sex=='f')*(1-dum$p1))
  N.f.2 <- sum(popvect*as.numeric(dum$sex=='f')*(1-dum$p1))
  # offspring per island - based on #females
  pre_off_1 <- sum(dum$p1*as.numeric(dum$sex=='f')*popvect)# before survival penalty
  pre_off_2 <- sum((1-dum$p1)*as.numeric(dum$sex=='f')*popvect)
  
  off <- c(pre_off_1*surv_off(),pre_off_2*surv_off())
  
  # female fitness (depending on offspring survival rate)
  
  # per father
  male.dist.1 <- dum$p1*as.numeric(dum$sex=='m')*popvect
  male.dist.2 <- (1-dum$p1)*as.numeric(dum$sex=='m')*popvect
  
  
  # survival of adults
}