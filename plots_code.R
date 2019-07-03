setwd("C:/Users/Lara/Documents/fur_seals/lt2_2019-06-22") 

####### Plots for model: fur-seals project. Marginal-male-hypothesis I ######

### 1) Plot over time 

library(Hmisc)

overview <- read.csv("00-param-overview.csv",stringsAsFactors = FALSE)

dummy_file <- read.csv(overview$filename[19])[,-1]

all_data <- array(NA,dim=c(nrow(dummy_file),ncol(dummy_file),nrow(overview)))

for(i in 1:nrow(overview)){
  if(file.exists(overview$filename[i])){
    all_data[,,i] <- as.matrix(read.csv(overview$filename[i])[1:nrow(dummy_file),-1])
  }
}

overview$treatment <- as.numeric(factor(paste(overview$p)))

pdf("output.pdf",  width=12, height=8)

par(mfrow=c(2,3))
color <- c("darkgrey","lightcoral","lightblue")
color_trans <- adjustcolor(color, alpha.f = 0.2) 

for(i in unique(overview$treatment)){
  inds <- which(overview$treatment ==i)
  tmp <- all_data[,,inds]
  if(all(is.na(tmp))){
    tmp <- array(0,dim=dim(tmp))
  }
  
  #Decide which patch is high density, which low density patch (mean N over time and compare values N1/N2):
  high_patches <- c()
  low_patches <- c()
  finals <- c()
  
  for(b in 1:length(inds)){
    
    finals <- cbind(mean(tmp[,2,b]),mean(tmp[,3,b])) #get mean of column 2 and 3 (N patch 1 and N patch2)
    high_patches[b] <- which.max(apply(finals,MARGIN=2,max)) #get the patch with the highest mean value of N over time, safe that in vector
    low_patches[b] <-  which.min(apply(finals,MARGIN=2,min)) 
    
    if (low_patches[b]==1){  #one equals column 2
      tmp[,,b] <- tmp[,c(1,3,2,5,4,7,6,9,8,11,10,13,12,15,14),b] #interchange the columns (1 gets 2 etc. when patch is lower) --> column 2 is always the one with the high density
    }
    
  }
  
  
  #calculate sd 
  sd_N <- c()
  sd_Nhigh <- c()
  sd_Nlow <- c()
  sd_trait_high <- c()
  sd_trait_low <- c()
  sd_trait_high_f <- c()
  sd_trait_low_f <- c()
  sd_off_high <- c()
  sd_off_low <- c()
  sd_cov_high <- c()
  sd_cov_low <- c()
  sd_fem_high <- c()
  sd_fem_low <- c()
  sd_males_high <- c()
  sd_males_low <- c()
  
  for (n in 1:nrow(dummy_file)){
    
    sd_N[n] <- sd(tmp[n,1,],na.rm=TRUE) #sd for column 1, calculated by values of the replicates
    sd_Nhigh[n] <- sd(tmp[n,2,],na.rm=TRUE) #sd for N in patch high dens
    sd_Nlow[n] <-  sd(tmp[n,3,],na.rm=TRUE) #sd for N in patch low dens
    sd_trait_high[n] <- sd(tmp[n,4,],na.rm=TRUE) #trait sd patch high dens
    sd_trait_low[n] <- sd(tmp[n,5,],na.rm=TRUE) #trait sd patch low dens
    sd_trait_high_f[n] <- sd(tmp[n,6,],na.rm=TRUE) #trait sd patch high dens
    sd_trait_low_f[n] <- sd(tmp[n,7,],na.rm=TRUE) #trait sd patch low den
    sd_off_high[n] <- sd(tmp[n,12,],na.rm=TRUE) #offspring produced sd patch high dens
    sd_off_low[n] <- sd(tmp[n,13,],na.rm=TRUE) #offspring produced sd patch low dens
    sd_cov_high[n] <- sd(tmp[n,14,],na.rm=TRUE) #covariance sd patch high dens
    sd_cov_low[n] <- sd(tmp[n,15,],na.rm=TRUE) #cov sd patch low dens
    sd_fem_high[n] <- sd(tmp[n,10,],na.rm=TRUE) #covariance sd patch high dens
    sd_fem_low[n] <- sd(tmp[n,11,],na.rm=TRUE) #cov sd patch low dens
    sd_males_high[n] <- sd(tmp[n,8,],na.rm=TRUE) #covariance sd patch high dens
    sd_males_low[n] <- sd(tmp[n,9,],na.rm=TRUE) #cov sd patch low dens
    
  }
  
  par(xpd=NA,oma=c(3,0,0,0))  
  plot(1:nrow(dummy_file),apply(tmp[,1,],1,FUN = function(x) mean(x,na.rm=TRUE)),ylim=c(0,max(tmp,na.rm=TRUE)),type="l",lwd=0.1,xlab="Time",ylab="N",main=paste(overview[inds[1],c('p')],collapse=" "))
  lines(1:nrow(dummy_file),apply(tmp[,2,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightcoral")
  lines(1:nrow(dummy_file),apply(tmp[,3,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightblue")
  errbar(1:nrow(dummy_file), apply(tmp[,1,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,1,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_N, type="l", apply(tmp[,1,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_N, errbar.col=color_trans[1], lwd=0.1, col="darkgrey", cap=0.003, add=TRUE)
  errbar(1:nrow(dummy_file), apply(tmp[,2,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,2,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_Nhigh, type="l", apply(tmp[,2,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_Nhigh, errbar.col=color_trans[2], lwd=0.1,col="lightcoral", cap=0.003, add=TRUE)
  errbar(1:nrow(dummy_file), apply(tmp[,3,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,3,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_Nlow, type="l", apply(tmp[,3,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_Nlow, errbar.col=color_trans[3], lwd=0.1,col="lightblue", cap=0.003, add=TRUE)
  legend(par("usr")[1],par("usr")[3], xjust=0, yjust=2, bty = "n", xpd = TRUE, legend=c("total", "high density patch", "low density patch"), col=c("darkgrey","lightcoral", "lightblue"), lty=1, lwd=2,cex=0.6)
  
  plot(0,0,type="n",xlim=c(0,nrow(dummy_file)),ylim=c(0,35),xlab="Time",ylab="mean male trait")
  lines(1:nrow(dummy_file),apply(tmp[,4,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightcoral", lwd=0.1)
  lines(1:nrow(dummy_file),apply(tmp[,5,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightblue" ,lwd=0.1)
  errbar(1:nrow(dummy_file), apply(tmp[,4,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,4,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_trait_high, type="l", apply(tmp[,4,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_trait_high, errbar.col=color_trans[2],lwd=0.1,col="lightcoral", cap=0.003, add=TRUE)
  errbar(1:nrow(dummy_file), apply(tmp[,5,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,5,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_trait_low, type="l", apply(tmp[,5,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_trait_low, errbar.col=color_trans[3],lwd=0.1,col="lightblue", cap=0.003, add=TRUE)
  
  plot(0,0,type="n",xlim=c(0,nrow(dummy_file)),ylim=c(-2.2,1.5),xlab="Time",ylab="mean female trait")
  lines(1:nrow(dummy_file),apply(tmp[,6,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightcoral",lwd=0.1)
  lines(1:nrow(dummy_file),apply(tmp[,7,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightblue",lwd=0.1)
  errbar(1:nrow(dummy_file), apply(tmp[,6,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,6,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_trait_high_f, type="l", apply(tmp[,6,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_trait_high_f, errbar.col=color_trans[2],lwd=0.1, col="lightcoral", cap=0.003, add=TRUE)
  errbar(1:nrow(dummy_file), apply(tmp[,7,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,7,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_trait_low_f, type="l", apply(tmp[,7,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_trait_low_f, errbar.col=color_trans[3],lwd=0.1, col="lightblue",cap=0.003, add=TRUE)
  
  # plot(0,0,type="n",xlim=c(0,nrow(dummy_file)),ylim=c(-5,5),xlab="Time",ylab="covariance male offspring ~ trait")
  # lines(1:nrow(dummy_file),apply(tmp[,14,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="blue")
  # lines(1:nrow(dummy_file),apply(tmp[,15,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="red")
  # errbar(1:nrow(dummy_file), apply(tmp[,14,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,12,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_cov1, type="l", apply(tmp[,12,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_cov1, errbar.col="lightcoral", cap=0.003, add=TRUE)
  # errbar(1:nrow(dummy_file), apply(tmp[,15,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,13,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_cov2, type="l", apply(tmp[,13,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_cov2, errbar.col="lightblue", cap=0.003, add=TRUE)
  
  plot(0,0,type="n",xlim=c(0,nrow(dummy_file)),ylim=c(-1,80),xlab="Time",ylab="offspring produced")
  lines(1:nrow(dummy_file),apply(tmp[,12,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightcoral",lwd=0.1)
  lines(1:nrow(dummy_file),apply(tmp[,13,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightblue",lwd=0.1)
  errbar(1:nrow(dummy_file), apply(tmp[,12,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,12,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_off_high, type="l", apply(tmp[,12,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_off_high, errbar.col=color_trans[2],lwd=0.1, col="lightcoral" ,cap=0.003, add=TRUE)
  errbar(1:nrow(dummy_file), apply(tmp[,13,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,13,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_off_low, type="l", apply(tmp[,13,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_off_low, errbar.col=color_trans[3],lwd=0.1,col="lightblue", cap=0.003, add=TRUE)
  
  plot(0,0,type="n",xlim=c(0,nrow(dummy_file)),ylim=c(-1,350),xlab="Time",ylab="females per patch")
  lines(1:nrow(dummy_file),apply(tmp[,10,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightcoral",lwd=0.1)
  lines(1:nrow(dummy_file),apply(tmp[,11,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightblue",lwd=0.1)
  errbar(1:nrow(dummy_file), apply(tmp[,10,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,10,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_fem_high, type="l", apply(tmp[,10,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_fem_high, errbar.col=color_trans[2],lwd=0.1, col="lightcoral" ,cap=0.003, add=TRUE)
  errbar(1:nrow(dummy_file), apply(tmp[,11,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,11,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_fem_low, type="l", apply(tmp[,11,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_fem_low, errbar.col=color_trans[3],lwd=0.1,col="lightblue", cap=0.003, add=TRUE)
  
  plot(0,0,type="n",xlim=c(0,nrow(dummy_file)),ylim=c(0,80),xlab="Time",ylab="males per patch")
  lines(1:nrow(dummy_file),apply(tmp[,8,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightcoral",lwd=0.1)
  lines(1:nrow(dummy_file),apply(tmp[,9,],1,FUN = function(x) mean(x,na.rm=TRUE)),col="lightblue",lwd=0.1)
  errbar(1:nrow(dummy_file), apply(tmp[,8,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,8,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_males_high, type="l", apply(tmp[,8,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_males_high, errbar.col=color_trans[2],lwd=0.1, col="lightcoral" ,cap=0.003, add=TRUE)
  errbar(1:nrow(dummy_file), apply(tmp[,9,],1,FUN = function(x) mean(x,na.rm=TRUE)), apply(tmp[,9,],1,FUN = function(x) mean(x,na.rm=TRUE))+sd_males_low, type="l", apply(tmp[,9,],1,FUN = function(x) mean(x,na.rm=TRUE))-sd_males_low, errbar.col=color_trans[3],lwd=0.1,col="lightblue", cap=0.003, add=TRUE)
  
}

dev.off() #safe all plots in one file


#############################################################################

#### 2) Plot over parameter values 

library(Hmisc)
library(lattice)
library(latticeExtra)
library(tidyr)
library(ggplot2)


overview$treatment <- as.numeric(factor(paste(overview$p)))
overview$finalN <- as.numeric(factor(paste(all_data[nrow(dummy_file),1,])))
overview$finalN1 <- as.numeric(factor(paste(all_data[nrow(dummy_file),2,]))) #final N for patch1
overview$finalN2 <- as.numeric(factor(paste(all_data[nrow(dummy_file),3,]))) #final N for patch2

overview$sd_N <- rep(NA, 1, nrow(overview))
overview$sd_Nhigh <- rep(NA, 1, nrow(overview)) #for high density patch
overview$sd_Nlow <- rep(NA, 1, nrow(overview)) #for low density patch
overview$delta_patches <- rep(NA, 1, nrow(overview)) #for differences in N of the patches (over time averaged)
overview$delta_trait <- rep(NA, 1, nrow(overview)) #for female trait 
overview$delta_trait_m <- rep(NA, 1, nrow(overview)) #for male trait 

for(i in 1:nrow(overview)){
  finals <- cbind(overview[overview$X==i,]$finalN1,overview[overview$X==i,]$finalN2)
  high <- which.max(apply(finals,MARGIN=2,max))
  low <-  which.min(apply(finals,MARGIN=2,min))
  overview[overview$X==i,]$sd_N <- sd(all_data[,1,i])
  overview[overview$X==i,]$sd_Nhigh <- sd(all_data[,high,i])
  overview[overview$X==i,]$sd_Nlow <- sd(all_data[,low,i])
}


for(i in 1:nrow(overview)){
  delta_p <- c()
  delta_t <- c()
  delta_t_m <- c()
  for(u in 1:nrow(dummy_file)){
    delta_p[u] <- abs((all_data[u,2,i])-(all_data[u,3,i]))
    delta_t[u] <- abs((all_data[u,10,i])-(all_data[u,11,i])) #difference in female trait value of patch 1 and 2
    delta_t_m[u] <- abs((all_data[u,8,i])-(all_data[u,9,i])) #difference in male trait value of patch 1 and 2
  }
  delta_p <- mean(delta_p)
  overview[overview$X==i,]$delta_patches <- delta_p
  delta_t <- mean(delta_t)
  overview[overview$X==i,]$delta_trait <- delta_t
  delta_t_m <- mean(delta_t_m)
  overview[overview$X==i,]$delta_trait_m <- delta_t_m
}


pdf("output-parameter-plot.pdf",width=10)

g1 <- ggplot(overview,aes(x=p,y=delta_patches,group=p)) + ylab("difference pop size patches") + xlab("p (strength of philopatry) \n assumed average density = 100") + geom_boxplot(fill = "#7ec0ee", color = "#6e7b8b") +  geom_point(alpha = 0.25) + labs(title="Delta difference in population sizes: varying p (x axis) for different total survival rates")  
g2 <- ggplot(overview,aes(x=p,y=sd_Nhigh, group=p )) + ylab("sd N patch high density") + geom_boxplot(fill = "#7ec0ee", color = "#6e7b8b") +xlab("p (strength of philopatry) \n assumed average density = 100")+labs(title="Standard deviation over time in the patch with higher N on average") 
g3 <- ggplot(overview,aes(x=p,y=sd_Nlow, group=p)) + ylab("sd N patch low density") + geom_boxplot(fill = "#7ec0ee", color = "#6e7b8b")+ xlab("p (strength of philopatry) \n assumed average density = 100")+labs(title="Standard deviation over time in the patch with lower N on average") 
g4 <- ggplot(overview,aes(x=p,y=finalN1, group=p)) + ylab("Final N1") + geom_boxplot(fill = "#7ec0ee", color = "#6e7b8b")+xlab("p (strength of philopatry) \n assumed average density = 100")+labs("Final population size of patch 1")+expand_limits(y=c(0,200))
g5 <- ggplot(overview,aes(x=p,y=finalN2, group=p)) + ylab("Final N2") + xlab("p (strength of philopatry) \n assumed average density = 100") + geom_boxplot(fill = "#7ec0ee", color = "#6e7b8b") +  facet_wrap(~ surv)+labs("Final population size of patch 1")+expand_limits(y=c(0,200))
g6 <- ggplot(overview,aes(x=p,y=finalN, group=p)) + ylab("Final pop size") + geom_boxplot(fill = "#7ec0ee", color = "#6e7b8b")+xlab("p (strength of philopatry) \n assumed average density = 100")+labs(title="Final population size: varying p (x axis) for different total survival rates")+expand_limits(y=c(0,300))
g7 <- ggplot(overview,aes(x=p,y=delta_trait, group=p)) + ylab("difference in female trait values") +  xlab("p (strength of philopatry) \n assumed average density = 100") + geom_boxplot(fill = "#7ec0ee", color = "#6e7b8b") +  labs(title="Difference in trait values: varying p (x axis) for different total survival rates")
g8 <- ggplot(overview,aes(x=p,y=delta_trait_m, group=p)) + ylab("difference in male trait values") +  xlab("p (strength of philopatry) \n assumed average density = 100") + geom_boxplot(fill = "#7ec0ee", color = "#6e7b8b") +  labs(title="Difference in trait values: varying p (x axis) for different total survival rates")

print(g1)
print(g2)
print(g3)
print(g4)
print(g5)
print(g6)
print(g7)
print(g8)

dev.off() #safe all plot in one file
