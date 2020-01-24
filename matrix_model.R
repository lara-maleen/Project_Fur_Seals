## Matrix-like model for the Fur Seals

popvect <- rep(1,18)
dum <- expand.grid(c('a','A'),c('a','A'),'-',c('b','B'),c('b','B'),'-',c('m','f'))
dum <- dum[!(dum[,1] == 'A' & dum[,2] == 'a') & !(dum[,4] =='B' & dum[,5] == 'b'),]
colnames(dum) <- c('a1','a2','dash1','b1','b2','dash2','sex')
catnames <- as.character(apply(dum,1, FUN = function(x) paste(x,sep="",collapse="")))
dum2 <- dum
min_val <- 0.1
dum2$p1 <- min_val + 0.5*(1-2*min_val)*(as.numeric(dum2$sex=='m')*(as.numeric(dum2$a1=='A') + as.numeric(dum2$a2=='A')) + as.numeric(dum2$sex=='f')*(as.numeric(dum2$b1=='B') + as.numeric(dum2$b2=='B'))) # probability of going to island 1 for each indiv
names(popvect) <- catnames

N.males <- function(popvect,dum,island){
  # inds <- dum$sex =='m' &  
}

make_mat <- function(popvect,dum){

  
  A <- matrix(0,nrow=18,ncol=18)
  rownames(A) <- catnames
  colnames(A) <- catnames
  A[cbind(2:length(popvect),1:(length(popvect)-1))] <- 0.99
  A[length(popvect),length(popvect)] <- 0.99
  A[1,length(popvect)] <- 2
  return(A)
}

# parameters
Nt <- 200
store <- matrix(NA,nrow=length(popvect),ncol=Nt)
store[,1] <- popvect
for(t in 2:Nt){
  A <- make_mat(store[,t-1],catnames)
  store[,t] <- A %*% store[,t-1] 
}


normalized <- t(store)/colSums(store)
df <- data.frame(val = as.numeric(normalized), cat = factor(rep(catnames,each=Nt),levels = catnames), time = rep(1:Nt,length(popvect)))
library(ggplot2)

ggplot(df, aes(x=time,y=val,fill=cat)) + geom_col()
round(normalized[nrow(normalized),],2)
round(Re(eigen(A)$vector[,1]/sum(eigen(A)$vector[,1])),2)
