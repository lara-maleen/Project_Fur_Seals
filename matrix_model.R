## Matrix-like model for the Fur Seals

popvect <- runif(18)
# popvect <- rep(0,18)
# popvect[c(1,10)] <- 1
dum <- expand.grid(c('a','A'),c('a','A'),'-',c('b','B'),c('b','B'),'-',c('m','f'),stringsAsFactors = FALSE)
dum <- dum[!(dum[,1] == 'A' & dum[,2] == 'a') & !(dum[,4] =='B' & dum[,5] == 'b'),]
colnames(dum) <- c('a1','a2','dash1','b1','b2','dash2','sex')
catnames <- as.character(apply(dum,1, FUN = function(x) paste(x,sep="",collapse="")))
dum2 <- dum
min_val_m <- 0.1
min_val_f <- 0.01
dum2$p1 <- min_val_m*as.numeric(dum2$sex=='m') + min_val_f*as.numeric(dum2$sex == 'f') + 
  0.5*(1-2*min_val_m)*(as.numeric(dum2$sex=='m')*(as.numeric(dum2$a1=='A') + as.numeric(dum2$a2=='A'))) + 0.5*(1-2*min_val_f)*(as.numeric(dum2$sex=='f')*(as.numeric(dum2$b1=='B') + as.numeric(dum2$b2=='B'))) # probability of going to island 1 for each indiv
dum2$Ascore <- as.numeric(dum2$a1 == 'A') + as.numeric(dum2$a2 == 'A')
dum2$Bscore <- as.numeric(dum2$b1 == 'B') + as.numeric(dum2$b2 == 'B')

names(popvect) <- catnames
surv <- 0
A.adv <- 1.5

calc_off_dist <- function(fem,dum,male.dist){
  # print(summary(male.dist))
  #which alleles does the male contribute
  # cat(sum(male.dist),"\n")
  male.dist[dum$Ascore == 2] <- A.adv*male.dist[dum$Ascore == 2]
  male.dist <- male.dist/sum(male.dist)
  
  df.male <- rbind(data.frame(a1=dum$a1,b1=dum$b1,pfat=0.25*male.dist,stringsAsFactors = FALSE),
                   data.frame(a1=dum$a2,b1=dum$b1,pfat=0.25*male.dist,stringsAsFactors = FALSE),
                   data.frame(a1=dum$a1,b1=dum$b2,pfat=0.25*male.dist,stringsAsFactors = FALSE),
                   data.frame(a1=dum$a2,b1=dum$b2,pfat=0.25*male.dist,stringsAsFactors = FALSE))
  # print(summary(df.male))
  df.offs <- as.data.frame(rbind(df.male,df.male,df.male,df.male),stringsAsFactors=FALSE)
  # print(df.offs)
  colnames(df.offs) <- c('a1','b1','pfat')
  
  df.offs$a2 <- rep(rep(c(fem$a1,fem$a2),each=nrow(df.male)),2)
  df.offs$b2 <- rep(c(fem$b1,fem$b2),each=2*nrow(df.male))

  df.offs$pmat <- 0.25
  
  offs <- rep(0,length(male.dist))
  
  df.offs$Ascore <- as.numeric(df.offs$a1 == 'A') + as.numeric(df.offs$a2 == 'A')
  df.offs$Bscore <- as.numeric(df.offs$b1 == 'B') + as.numeric(df.offs$b2 == 'B')
    # print(df.offs)
    # stop("err")
  for(i in 1:nrow(df.offs)){
    loc <- which(df.offs$Ascore[i] == dum$Ascore & df.offs$Bscore[i] == dum$Bscore)
    # print(loc)
    # cat(df.offs$pfat[i]*df.offs$pmat[i],"\n")
    # print(summary(df.offs$pfat))
    # print(summary(df.offs$pmat))
    offs[loc] <- offs[loc] + df.offs$pmat[i]*df.offs$pfat[i]/2
    
  }
  # cat(sum(offs),"\n")
  offs
}

make_mat <- function(popvect,dum){

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
  rownames(A) <- catnames
  colnames(A) <- catnames
  #A[cbind(2:length(popvect),1:(length(popvect)-1))] <- 0.99
  #A[length(popvect),length(popvect)] <- 0.99
  #A[1,length(popvect)] <- 2
  diag(A) <- surv
  N.f.1 <- sum((dum$p1*popvect)[dum$sex == 'f'])
  N.f.2 <- sum(((1-dum$p1)*popvect)[dum$sex == 'f'])
  surv.1 <- 0.9+0.1*plogis(100*(0.25-N.f.1))
  surv.2 <- 0.9+0.1*plogis(100*(0.25-N.f.2))
  
  for(i in which(dum$sex == 'f')){
    offs.dist.1 <- calc_off_dist(dum[i,],dum,male.dist.1)
    offs.dist.2 <- calc_off_dist(dum[i,],dum,male.dist.2)

    A[,i] <- A[,i] + dum$p1[i]*offs.dist.1*surv.1 + (1-dum$p1[i])*offs.dist.2*surv.2

  }
  cat(surv.1,"\t",surv.2,"\n")
  return(A)
}

# parameters
Nt <- 1000
store <- matrix(NA,nrow=length(popvect),ncol=Nt)
store[,1] <- popvect / sum(popvect)
for(t in 2:Nt){
  A <- make_mat(store[,t-1],dum2)
  store[,t] <- A %*% store[,t-1] 
  store[,t] <- store[,t]/sum(store[,t])
}


normalized <- t(store)/colSums(store)
df <- data.frame(val = as.numeric(normalized), cat = factor(rep(catnames,each=Nt),levels = catnames), time = rep(1:Nt,length(popvect)))
library(ggplot2)

# All frequencies
p0 <- ggplot(df, aes(x=time,y=val,fill=cat))+ geom_col() + theme_bw()
round(normalized[nrow(normalized),],2)
round(Re(eigen(A)$vector[,1]/sum(eigen(A)$vector[,1])),2)

# Num inds per patch
N.m.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$sex=='m']))
N.m.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$sex=='m']))
N.f.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$sex=='f']))
N.f.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$sex=='f']))

df.sex.patch <- data.frame(N = c(N.m.1,N.m.2,N.f.1,N.f.2), t = rep(1:length(N.m.1),4) ,  sex = c(rep('m',2*length(N.m.1)),rep('f',2*length(N.f.1))),patch = c(rep(c(1,2,1,2),each=length(N.m.1))))

p1 <- ggplot(df.sex.patch, aes(x=t,y=N,col = sex, lty = factor(patch))) + geom_line() + theme_bw() + ylim(c(0,1))

# tot ind per patch
N.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)))
N.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)))

df.patch <- data.frame(t = rep(1:length(N.1),2), N = c(N.1,N.2), patch = factor(rep(c(1,2),each=length(N.1))))

p2 <- ggplot(df.patch, aes(x=t, y=N, lty = patch)) + geom_line() + theme_bw() + ylim(c(0,1))


# Genotypes per patch: A and a (males and females)
N.aa.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Ascore==0]))
N.Aa.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Ascore==1]))
N.AA.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Ascore==2]))

N.aa.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Ascore==0]))
N.Aa.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Ascore==1]))
N.AA.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Ascore==2]))

df.Aa <- data.frame(t = rep(1:length(N.aa.1), 6), N = c(N.aa.1,N.Aa.1,N.AA.1,N.aa.2,N.Aa.2,N.AA.2), patch = factor(rep(c(1,2),each=3*length(N.aa.1))), alleles = factor(rep(c('aa','Aa','AA','aa','Aa','AA'),each = length(N.aa.1))))

p3 <- ggplot(df.Aa, aes(x=t, y=N, lty = patch, col = alleles)) + geom_line() + theme_bw() + ylim(c(0,1))

# Genotypes per patch: B and b (males and females)
N.bb.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Bscore==0]))
N.Bb.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Bscore==1]))
N.BB.1 <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Bscore==2]))

N.bb.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Bscore==0]))
N.Bb.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Bscore==1]))
N.BB.2 <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Bscore==2]))

df.Bb <- data.frame(t = rep(1:length(N.bb.1), 6), N = c(N.bb.1,N.Bb.1,N.BB.1,N.bb.2,N.Bb.2,N.BB.2), patch = factor(rep(c(1,2),each=3*length(N.bb.1))), alleles = factor(rep(c('bb','Bb','BB','bb','Bb','BB'),each = length(N.aa.1))))

p4 <- ggplot(df.Bb, aes(x=t, y=N, lty = patch, col = alleles)) + geom_line() + theme_bw() + ylim(c(0,1))

# Genotypes per patch: A and a (males)
N.aa.1.m <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Ascore==0 & dum2$sex == 'm']))
N.Aa.1.m <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Ascore==1 & dum2$sex == 'm']))
N.AA.1.m <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Ascore==2 & dum2$sex == 'm']))

N.aa.2.m <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Ascore==0 & dum2$sex == 'm']))
N.Aa.2.m <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Ascore==1 & dum2$sex == 'm']))
N.AA.2.m <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Ascore==2 & dum2$sex == 'm']))

df.Aa.m <- data.frame(t = rep(1:length(N.aa.1), 6), N = c(N.aa.1.m,N.Aa.1.m,N.AA.1.m,N.aa.2.m,N.Aa.2.m,N.AA.2.m), patch = factor(rep(c(1,2),each=3*length(N.aa.1))), alleles = factor(rep(c('aa','Aa','AA','aa','Aa','AA'),each = length(N.aa.1))))

p5 <- ggplot(df.Aa.m, aes(x=t, y=N, lty = patch, col = alleles)) + geom_line() + theme_bw() + ylim(c(0,1))

# Genotypes per patch: B and b (females)
N.bb.1.f <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Bscore==0 & dum2$sex == 'f']))
N.Bb.1.f <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Bscore==1 & dum2$sex == 'f']))
N.BB.1.f <- apply(normalized,1,FUN = function(x) sum((dum2$p1*x)[dum2$Bscore==2 & dum2$sex == 'f']))

N.bb.2.f <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Bscore==0 & dum2$sex == 'f']))
N.Bb.2.f <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Bscore==1 & dum2$sex == 'f']))
N.BB.2.f <- apply(normalized,1,FUN = function(x) sum(((1-dum2$p1)*x)[dum2$Bscore==2 & dum2$sex == 'f']))

df.Bb.f <- data.frame(t = rep(1:length(N.bb.1), 6), N = c(N.bb.1.f,N.Bb.1.f,N.BB.1.f,N.bb.2.f,N.Bb.2.f,N.BB.2.f), patch = factor(rep(c(1,2),each=3*length(N.bb.1))), alleles = factor(rep(c('bb','Bb','BB','bb','Bb','BB'),each = length(N.aa.1))))

p6 <- ggplot(df.Bb.f, aes(x=t, y=N, lty = patch, col = alleles)) + geom_line() + theme_bw() + ylim(c(0,1))
library(gridExtra)
grid.arrange(p0,p1,p2,p3,p4,p5,p6,nrow=3)
