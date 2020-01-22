## Matrix-like model for the Fur Seals

popvect <- rep(1,18)
dum <- expand.grid(c('a','A'),c('a','A'),'-',c('b','B'),c('b','B'),'-',c('m','f'))
dum <- dum[!(dum[,1] == 'A' & dum[,2] == 'a') & !(dum[,4] =='B' & dum[,5] == 'b'),]
catnames <- as.character(apply(dum,1, FUN = function(x) paste(x,sep="",collapse="")))
A <- matrix(0,nrow=18,ncol=18)
names(popvect) <- catnames
rownames(A) <- catnames
colnames(A) <- catnames
A
popvect
A[1,1] <- 1.01
A[2,2] <- 0.99
# diag(A) <- seq(0,1,length.out = length(popvect))

A %*% (A %*% popvect)


# parameters

Nt <- 100
store <- matrix(NA,nrow=length(popvect),ncol=Nt)
store[,1] <- popvect
for(t in 2:Nt){
 store[,t] <- A %*% store[,t-1] 
}


normalized <- t(store)/colSums(store)
df <- data.frame(val = as.numeric(normalized), cat = factor(rep(catnames,each=Nt),levels = catnames), time = rep(1:Nt,length(popvect)))
library(ggplot2)

ggplot(df, aes(x=time,y=val,fill=cat)) + geom_col()
