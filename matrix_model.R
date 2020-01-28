## Matrix-like model for the Fur Seals

popvect <- runif(18)
# popvect <- rep(0,18)
# popvect[c(1,10)] <- 1
dum <- expand.grid(c('a','A'),c('a','A'),'-',c('b','B'),c('b','B'),'-',c('m','f'),stringsAsFactors = FALSE)
dum <- dum[!(dum[,1] == 'A' & dum[,2] == 'a') & !(dum[,4] =='B' & dum[,5] == 'b'),]
colnames(dum) <- c('a1','a2','dash1','b1','b2','dash2','sex')
catnames <- as.character(apply(dum,1, FUN = function(x) paste(x,sep="",collapse="")))
dum2 <- dum
min_val <- 0.1
dum2$p1 <- min_val + 0.5*(1-2*min_val)*(as.numeric(dum2$sex=='m')*(as.numeric(dum2$a1=='A') + as.numeric(dum2$a2=='A')) + as.numeric(dum2$sex=='f')*(as.numeric(dum2$b1=='B') + as.numeric(dum2$b2=='B'))) # probability of going to island 1 for each indiv
dum2$Ascore <- as.numeric(dum2$a1 == 'A') + as.numeric(dum2$a2 == 'A')
dum2$Bscore <- as.numeric(dum2$b1 == 'B') + as.numeric(dum2$b2 == 'B')

names(popvect) <- catnames
surv <- 0.9

calc_off_dist <- function(fem,dum,male.dist){
  # print(summary(male.dist))
  #which alleles does the male contribute
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
  for(i in which(dum$sex == 'f')){
    offs.dist.1 <- calc_off_dist(dum[i,],dum,male.dist.1)
    offs.dist.2 <- calc_off_dist(dum[i,],dum,male.dist.2)
    A[,i] <- A[,i] + dum$p1[i]*offs.dist.1 + (1-dum$p1[i])*offs.dist.2
    # print(offs.dist.2)
  }

  # print(A)
  # stop("err")
  return(A)
}

# parameters
Nt <- 100
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

ggplot(df, aes(x=time,y=val,fill=cat)) + geom_col(col="black")
round(normalized[nrow(normalized),],2)
round(Re(eigen(A)$vector[,1]/sum(eigen(A)$vector[,1])),2)
