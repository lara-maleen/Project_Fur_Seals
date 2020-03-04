library(ggplot2)
library(gridExtra)

plotting_rep <- function(filename){
  
  store <- read.csv(paste(filename,".csv",sep=""))
  dum2 <- read.csv(paste(filename,".dum",sep=""))
  
  normalized <- t(store)/colSums(store)
  df <- data.frame(val = as.numeric(normalized), cat = factor(rep(catnames,each=Nt),levels = catnames), time = rep(1:Nt,length(popvect)))

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
  
  pdf(paste(filename,".pdf",sep=""))
    grid.arrange(p0,p1,p2,p3,p4,p5,p6,nrow=3)
  dev.off()
  return(NA)
}

plotting_sum <- function(stats){
  
}
