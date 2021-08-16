rm(list=ls())
maindir <- "~/Documents/projects/Project_Fur_Seals/Matrix-Model/"

combset <- 1 # to prevent having to regroup the simulations always
# 'out-maxiscan-3/'="plain",'out-randomfatherscan-3/'='random father',
# 'out-maxiscan-2/'="plain",'out-randomfatherscan-2/'='random father',
# 'out-maxiscan-1/'="plain",'out-randomfatherscan-1/'='random father',
if(combset==1){
  subdirs <- c("out-maxiscan-2/","out-randomfatherscan-1/")
}else if(combset==2){
  subdirs <- c("out-maxiscan-1/","out-randomfatherscan-2/")
}else if(combset==3){
  subdirs <- c("out-maxiscan-3/","out-randomfatherscan-3/")
}
odir <- paste(maindir,subdirs[1],sep="")

df <- read.csv(paste(maindir,subdirs[1],"/simruns.stat",sep=""))
if(!"random_father" %in% colnames(df)){
  df <- cbind(df[,1:12],FALSE,df[,13:34])
  colnames(df)[13] <- "random_father"
}
df$sbd <- subdirs[1]
if(length(subdirs) > 1){
for(fl in subdirs[-1]){
  tmp <- read.csv(paste(maindir,fl,"/simruns.stat",sep=""))
  if(!"random_father" %in% colnames(tmp)){
    tmp <- cbind(tmp[,1:12],FALSE,tmp[,13:34])
    colnames(tmp)[13] <- "random_father"
  }
  tmp$sbd <- fl
  df <- rbind(df,tmp)
}
}
# df <- df[df$Apenalty!=0.30,]
the_big_rename <- function(x){
  inform <- list(Apenalty = "Male survival penalty",a2="Male advantage",sval="Density effect on offspring",mvm="Uncertainty male patch choice",mvf="Uncertainty female patch choice")
  
  x[x%in%names(inform)] <- inform[x[x%in%names(inform)]]
  
  x
}

varvar <- lapply(as.list(df[,c('a2','sval','sslope','stype','wm','d','d2','mvm','mvf','Apenalty')]),FUN = function(x) if(any(is.character(x))){!all(duplicated(x)[-1])}else{return(length(unique(round(x,2)))>3)})
varvar2 <- unlist(varvar)
varvar3 <- names(varvar2)[varvar2]
tmp <- expand.grid(var1=1:length(varvar3),var2=1:length(varvar3))
tmp <- tmp[tmp$var1<tmp$var2,]
library(ggplot2)
library(ggnewscale)
# creating variables
df$DNf <- df$N.1.f - df$N.2.f
df$DNm <- df$N.1.m - df$N.2.m
df$DN <- df$N1 - df$N2

df$r1 <- df$N.1.f/df$N.1.m
df$r2 <- df$N.2.f/df$N.2.m

library(tidyr)
df2 <- gather(df,condition,val,DNf:r2)

pdf(paste(odir,"comb-graph.pdf",sep=""),width=6,height=4)
for(i in 1:nrow(tmp)){
  var1 <- varvar3[tmp$var1[i]]
  var2 <- varvar3[tmp$var2[i]]
  df2$xvar <- df2[,var1]
  df2$yvar <- df2[,var2]
  
  ## check spread of replicates
  p1 <- ggplot(df2,aes(x=xvar,group=paste(xvar,yvar),col=yvar,y=val)) + geom_boxplot(outlier.size = 0.2) +theme_classic() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    xlab(the_big_rename(var1)) + ylab("value") + scale_color_gradient(name=the_big_rename(var2))+facet_grid(sbd~condition,scales = "free")#scale_fill_gradient2(name=var2,low="blue",high="red",mid="grey",midpoint=0)
  print(p1)
  
  ## plot first replicate (faster now - later mean is probably better)
  labfun <- function(x){
    c(DN="Total",DNf="Females",DNm="Males",r1="island 1",r2="island 2",
      'out-maxiscan-3/'="plain",'out-randomfatherscan-3/'='random father',
      'out-maxiscan-2/'="plain",'out-randomfatherscan-2/'='random father',
      'out-maxiscan-1/'="plain",'out-randomfatherscan-1/'='random father',
      'out-maxiscan-1-lowpenalty/'='weak density effect','out-maxiscan-2-lowpenalty/'='weak density effect',
      'out-maxiscan-2-long-lowpanelty/'='10x longer',
      'out-offspringsurvives-1/'="no survival penalty")[x]
  }  
  df3 <- aggregate(df2$val,by=list(condition=df2$condition,a2=df2$a2,sval=df2$sval,sslope=df2$sslope,stype=df2$stype,wm=df2$wm,wf=df2$wf,d=df2$d,d2=df2$d2,mvm=df2$mvm,mvf=df2$mvf,Apenalty=df2$Apenalty,sbd=df2$sbd),mean)
  colnames(df3)[colnames(df3)=='x'] <- 'val'
  df3$xvar <- df3[,var1]
  df3$yvar <- df3[,var2]
  p2 <- ggplot(df3[df3$condition %in% c('DN','DNf','DNm'),],aes(x=xvar,y=yvar)) +theme_classic() + scale_x_continuous(expand=c(0,0),n.breaks = 3) + scale_y_continuous(expand=c(0,0))+
    xlab(the_big_rename(var1)) + ylab(the_big_rename(var2)) #+ scale_fill_gradient2(name="Difference",high="orange",low="dodgerblue",mid="grey",midpoint = 0)
  p2 <- p2 + geom_tile(aes(fill=val), data = ~ subset(., sbd == subdirs[1])) + scale_fill_gradient2(name="Difference 1",high="orange",low="dodgerblue",mid="grey",midpoint = 0,guide=guide_colorbar(order=1))
  p2 <- p2 + new_scale_fill()
  p2 <- p2 + geom_tile(aes(fill=val), data = ~ subset(., sbd == subdirs[2])) 
  p2 <- p2 + scale_fill_gradient2(name="Difference 2",high="orange",low="dodgerblue",mid="grey",midpoint = 0,guide=guide_colorbar(order=2))
  p2 <- p2 + facet_grid(sbd~condition,scales = "free",labeller=as_labeller(labfun))+theme(panel.spacing.x = unit(4, "mm"))
  #scale_fill_gradient2(name=var2,low="blue",high="red",mid="grey",midpoint=0)
  print(p2)  
  

  p3 <- ggplot(df3[df3$condition %in% c('r1','r2'),],aes(x=xvar,y=yvar)) +theme_classic() + scale_x_continuous(expand=c(0,0),n.breaks = 3) + scale_y_continuous(expand=c(0,0))+
    xlab(the_big_rename(var1)) + ylab(the_big_rename(var2)) 
  p3 <- p3  + geom_tile(aes(fill=val), data = ~ subset(., sbd == subdirs[1])) + scale_fill_gradient2(trans="log",labels = function(x)round(x,2),name="Ratio f/m 1",high="orange",low="dodgerblue",mid="grey",midpoint = 0,guide=guide_colorbar(order=1))
  p3 <- p3  + new_scale_fill() + geom_tile(aes(fill=val), data = ~ subset(., sbd == subdirs[2])) + scale_fill_gradient2(trans="log",labels = function(x)round(x,2),name="Ratio f/m 2",high="orange",low="dodgerblue",mid="grey",midpoint = 0,guide=guide_colorbar(order=2))
    p3 <- p3 +facet_grid(sbd~condition,scales = "free",labeller = as_labeller(labfun))+theme(panel.spacing.x = unit(4, "mm"))#scale_fill_gradient2(name=var2,low="blue",high="red",mid="grey",midpoint=0)
  print(p3)  
  
}
dev.off()
