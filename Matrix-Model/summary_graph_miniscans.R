rm(list=ls())
odir <- "~/Documents/projects/Project_Fur_Seals/Matrix-Model/out-maxiscan-2"
df <- read.csv(paste(odir,"/simruns.stat",sep=""))
# df <- df[df$Apenalty!=0.30,]
the_big_rename <- function(x){
  inform <- list(Apenalty = "Male survival penalty",a2="Male advantage",sval="Density effect on offspring",mvm="Uncertainty male patch choice",mvf="Uncertainty female patch choice")
  
  x[x%in%names(inform)] <- inform[x[x%in%names(inform)]]
  
  x
}
varvar <- lapply(as.list(df[,c('a2','sval','sslope','stype','wm','d','d2','mvm','mvf','Apenalty')]),FUN = function(x) if(any(is.character(x))){!all(duplicated(x)[-1])}else{return(var(x)>1e-6)})
varvar2 <- unlist(varvar)
varvar3 <- names(varvar2)[varvar2]
tmp <- expand.grid(var1=1:length(varvar3),var2=1:length(varvar3))
tmp <- tmp[tmp$var1<tmp$var2,]
library(ggplot2)

# creating variables
df$DNf <- df$N.1.f - df$N.2.f
df$DNm <- df$N.1.m - df$N.2.m
df$DN <- df$N1 - df$N2

df$r1 <- df$N.1.f/df$N.1.m
df$r2 <- df$N.2.f/df$N.2.m

library(tidyr)
df2 <- gather(df,condition,val,DNf:r2)

pdf(paste(odir,"/graph.pdf",sep=""),width=5,height=3)
for(i in 1:nrow(tmp)){
  var1 <- varvar3[tmp$var1[i]]
  var2 <- varvar3[tmp$var2[i]]
  df2$xvar <- df2[,var1]
  df2$yvar <- df2[,var2]
  
  ## check spread of replicates
  p1 <- ggplot(df2,aes(x=xvar,group=paste(xvar,yvar),col=yvar,y=val)) + geom_boxplot(outlier.size = 0.2) +theme_classic() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    xlab(the_big_rename(var1)) + ylab("value") + scale_color_gradient(name=the_big_rename(var2))+facet_wrap(~condition,scales = "free")#scale_fill_gradient2(name=var2,low="blue",high="red",mid="grey",midpoint=0)
  print(p1)
  
  ## plot first replicate (faster now - later mean is probably better)
  labfun <- function(x){
    c(DN="Total",DNf="Females",DNm="Males")[x]
  }  
  
  df3 <- aggregate(df2$val,by=list(condition=df2$condition,a2=df2$a2,sval=df2$sval,sslope=df2$sslope,stype=df2$stype,wm=df2$wm,wf=df2$wf,d=df2$d,d2=df2$d2,mvm=df2$mvm,mvf=df2$mvf,Apenalty=df2$Apenalty),mean)
  colnames(df3)[colnames(df3)=='x'] <- 'val'
  df3$xvar <- df3[,var1]
  df3$yvar <- df3[,var2]
  p2 <- ggplot(df3[df3$condition %in% c('DN','DNf','DNm'),],aes(x=xvar,y=yvar,fill=val)) + geom_tile() +theme_classic() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    xlab(the_big_rename(var1)) + ylab(the_big_rename(var2)) + scale_fill_gradient2(name="Difference",high="orange",low="dodgerblue",mid="grey",midpoint = 0)+facet_wrap(~condition,scales = "free",labeller=as_labeller(labfun))#scale_fill_gradient2(name=var2,low="blue",high="red",mid="grey",midpoint=0)
  print(p2)  
  
  labfun2 <- function(x){
    c(r1="island 1",r2="island 2")[x]
  }  
  p3 <- ggplot(df3[df3$condition %in% c('r1','r2'),],aes(x=xvar,y=yvar,fill=val)) + geom_tile() +theme_classic() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    xlab(the_big_rename(var1)) + ylab(the_big_rename(var2)) + scale_fill_gradient2(name="Ratio female/male",high="orange",low="dodgerblue",mid="grey",midpoint = 1)+facet_wrap(~condition,scales = "free",labeller = as_labeller(labfun2))#scale_fill_gradient2(name=var2,low="blue",high="red",mid="grey",midpoint=0)
  print(p3)  
  
}
dev.off()
