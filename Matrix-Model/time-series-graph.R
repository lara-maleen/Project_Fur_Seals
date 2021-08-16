sfun <- function(n) 0.1*plogis(7.5*(n-0.5))
sfun2 <- function(n) 0.05*plogis(7.5*(n-0.5))
source("~/Documents/projects/Project_Fur_Seals/Matrix-Model/matrix_model_kirk_version.R")
df01 <- read.csv("~/Documents/projects/Project_Fur_Seals/Matrix-Model/out-maxiscan-1/focal.csv")

dummy <- run_sim(Nt=100,surv_off=sfun,dumgen=TRUE)
out1 <- with(df01[1,],run_sim(A.adv = a2,wm=wm,wf=1-wm,d=d,d2=d2,Nt=1e3,min_val_m=mvm,min_val_f=mvf,Apenalty = Apenalty,surv_off=sfun2,test=TRUE,tol = -10))
# out2 <- run_sim(N0=runif(nrow(dummy)),min_val_m=0.05,min_val_f=0.05,Nt=100,surv_off=sfun2,test=TRUE,tol=-10,A.adv = 1.8,Apenalty = 0.3,wm = 1,wf=0,d=0.5,d2 = 0.5)

out2 <- run_sim(N0=runif(nrow(dummy)),min_val_m=0.05,min_val_f=0.05,Nt=2e3,surv_off=sfun2,test=TRUE,tol=-10,A.adv = 1.1,Apenalty = 0.45,wm = 0.5,wf=0.5,d=0.5,d2 = 0.5)

library(ggplot2)

lines_dist <- function(dat.m,dat.f,ht=0.1,ylm=c(0,1)){
  normy <- function(dat){
    tmp <- aggregate(dat$x,by=list(t=dat$t,sex=dat$sex,isle=dat$isle),sum)
    tmp$x[match(paste(dat$t,dat$sex,dat$isle),paste(tmp$t,tmp$sex,tmp$isle))]
  }
  dat.m$norm <- normy(dat.m)
  dat.f$norm <- normy(dat.f)
  dat.m$cont <- ht*dat.m$x/dat.m$norm
  dat.f$cont <- ht*dat.f$x/dat.f$norm
 
  
  # calculate heights of the segments for each genotype
  make_seg <- function(dat){
    # dat$cont[is.na(dat$cont)] <- 0
    dat$y0 <- NA
    dat$y1 <- NA
    dat$y0[dat$genotype==0] <- dat$norm[dat$genotype==0] - 0.5*ht
    dat$y1[dat$genotype==0] <- dat$y0[dat$genotype==0] + dat$cont[dat$genotype==0]
    dat$y1[dat$genotype==2] <- dat$norm[dat$genotype==0] + 0.5*ht
    dat$y0[dat$genotype==2] <- dat$y1[dat$genotype==2] - dat$cont[dat$genotype==2]
    newd0 <- dat[dat$genotype==0,]
    newd2 <- dat[dat$genotype==2,]
    newd1 <- dat[dat$genotype==1,]
    newd1$y0 <- newd0$y1[match(paste(newd1$t,newd1$sex,newd1$isle),paste(newd0$t,newd0$sex,newd0$isle))]
    newd1$y1 <- newd2$y0[match(paste(newd1$t,newd1$sex,newd1$isle),paste(newd2$t,newd2$sex,newd2$isle))]
    
    rbind(newd0,newd1,newd2)
  }
  dat.m2 <- make_seg(dat.m)
  dat.f2 <- make_seg(dat.f)
  tm2m <- rbind(dat.m2[dat.m2$isle==1,],dat.f2[dat.f2$isle==1,])#tm2[tm2$isle==1,]
  tm2f <- rbind(dat.m2[dat.m2$isle==2,],dat.f2[dat.f2$isle==2,])
  tm2m$isle <- factor(tm2m$isle,levels=1:2)
  tm2f$isle <- factor(tm2f$isle,levels=1:2)
  
  
  # p1 <- ggplot(tm2f,aes(x=t,fill=factor(genotype),y=cont)) + geom_col()  + facet_grid(sex~.) +theme_classic() + scale_fill_manual(values=c('blue','orange','grey',NA),name='genotype',labels=c('aa/bb','Aa/Bb','AA/BB','')) + ylab("N") + xlab("time") + 
  p1 <- ggplot(tm2f,aes(xmin=t-0.5,xmax=t+0.5,ymin=y0,ymax=y1,fill=factor(genotype))) + geom_rect()  + facet_grid(sex~.) +theme_classic() + scale_fill_manual(values=c('dodgerblue','orange','grey'),name='genotype',labels=c('aa/bb','Aa/Bb','AA/BB','')) + ylab("N") + xlab("time") + 
    geom_rect(data=tm2m,aes(xmin=t-0.5,xmax=t+0.5,ymin=y0,ymax=y1,fill=factor(genotype)),inherit.aes = FALSE) + 
     geom_line(data=tm2f[tm2f$genotype==0,],aes(x=t,y=y0,lty=isle),inherit.aes = FALSE,lwd=0.5) +
    geom_line(data=tm2f[tm2f$genotype==0,],aes(x=t,y=y0+ht,lty=isle),inherit.aes = FALSE,lwd=0.5)+
    geom_line(data=tm2m[tm2f$genotype==0,],aes(x=t,y=y0,lty=isle),inherit.aes = FALSE,lwd=0.5)+
        geom_line(data=tm2m[tm2f$genotype==0,],aes(x=t,y=y0+ht,lty=isle),inherit.aes = FALSE,lwd=0.5) +
    
    coord_cartesian(xlim =c(0,100), ylim = ylm) +scale_linetype(labels=c(1,2),name='island')
# pdf("testgraph.pdf")
print(p1)
# dev.off()
}

to_long <- function(dat,dum,gty='Ascore'){
  
  tmp <- apply(dat,2,FUN = function(x){
    o1 <- aggregate(x[-1]*dum$p1, by=list(genotype=dum[,gty],sex=dum$sex),sum)
    o1$isle <- 1
    o2 <- aggregate(x[-1]*(1-dum$p1), by=list(genotype=dum[,gty],sex=dum$sex),sum)
    o2$isle <- 2
    o3 <- rbind(o1,o2)
    o3$t <- x[1]
    o3
  })
  df_out <- do.call(rbind,tmp)
  df_out[,c('t','sex','genotype','isle','x')]
}

restructure_df <- function(dat,dum){
  # outputs a list of two data frames, one for the males, one for the females
  # dat <- out1$ts
  # dum <- out1$dum
  
  datm <- dat[c(1,1+which(dum$sex=='m')),]
  datf <- dat[c(1,1+which(dum$sex=='f')),]
  dumm <- dum[dum$sex=='m',]
  dumf <- dum[dum$sex=='f',]
  
  # 1. Males
  newm <- to_long(datm,dumm,'Ascore')
  # 2. Females
  newf <- to_long(datf,dumf,'Bscore')
  
  list(m=newm,f=newf)
}



#### option 1
# lines_dist(ldat$m,ldat$f,ht=0.1)

#### option 2: stacked bars
stbar <- function(ldat){
fuldat <- do.call(rbind,ldat)
fuldat <- fuldat[order(fuldat$isle),]
sumdat <- with(fuldat[fuldat$isle==2,],aggregate(x,by=list(t=t,sex=sex),sum)) 
darken <- function(cols){
  tmp <- col2rgb(cols)*0.85/255
  apply(tmp,2,FUN = function(x) rgb(x[1],x[2],x[3]))
}
p1 <- ggplot(fuldat,aes(x=t,y=x,fill=paste("i",isle,"g",genotype,sep=""))) + geom_col() +
  facet_grid(sex~.) + theme_classic() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
scale_fill_manual(values=c('dodgerblue','orange','grey',darken(c('dodgerblue','orange','grey'))),name='genotype',labels=c('aa/bb 1','Aa/Bb 1','AA/BB 1','aa/bb 2','Aa/Bb 2','AA/BB 2')) +
  geom_line(data=sumdat,aes(x=t,y=x,sex=sex),inherit.aes = FALSE,lwd=1,col="white") + xlab("time") + ylab("N")
print(p1)
}
# stbar(ldat)

#### option 3: line graph with added (genotype) frequency graphs
library(plyr)
library(grid)
gr31 <- function(ldat,sep_isle=FALSE){
  mdat <- ldat$m
  fdat <- ldat$f
  tdat <- mdat
  # tdat$x <- mdat$x + fdat$x[match(paste(tdat$t,tdat$isle),paste(fdat$t,fdat$isle))]
  magg <- aggregate(mdat$x,by=list(t=mdat$t,isle=mdat$isle,sex=mdat$sex),sum)
  fagg <- aggregate(fdat$x,by=list(t=fdat$t,isle=fdat$isle,sex=fdat$sex),sum)
  Nagg <- magg
  Nagg$x <- magg$x + fagg$x[match(paste(magg$t,magg$isle),paste(fagg$t,fagg$isle))]
  Nagg$sex <- "Total"
  # Nagg <- aggregate(tdat$x,by=list(t=tdat$t,isle=tdat$isle,sex=rep("Total",nrow(tdat))),sum)
  aggf <- with(rbind(magg,fagg,Nagg),data.frame(t=t,sex=sex,genotype=NA,isle=isle,x=x,var='N'))
  fdat$var <- 'f'
  mdat$var <- 'm'
  if(sep_isle){
    fdat$var <- paste('f, isle',fdat$isle)
    mdat$var <- paste('m, isle',mdat$isle)
  }
  fulldat <- rbind(mdat,fdat,aggf)
  if(!sep_isle){
    fulldat$var <- factor(fulldat$var,levels=c('f','N','m'))
  }else{
    fulldat$var <- factor(fulldat$var,levels=c('f, isle 1','f, isle 2','N','m, isle 1','m, isle 2'))
  }
  p1 <- ggplot(fulldat,aes(x=t,y=x,lty=factor(isle),col=factor(sex))) + facet_grid(var~.,scales="free_y") +
        scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
        theme_classic() + xlab("time") + ylab("N") + guides(linetype = guide_legend(override.aes=list(fill=NA,lwd=1))) +   theme(panel.border=element_rect(colour="black",size=1,fill=NA))
       
  p2 <- p1 + geom_line(data=fulldat[fulldat$var=='N',],aes(x=t,y=x,lty=factor(isle),col=factor(sex)),lwd=2,inherit.aes = TRUE) +
        scale_linetype_manual(name="island",values=c("solid", "22"))+
        scale_color_manual(name="sex",values=c("violetred","lightskyblue2","plum4")) +   theme(panel.border=element_rect(colour="black",size=1,fill=NA))
  p3 <- p2 + geom_col(data=fulldat[fulldat$var!='N',],aes(x=t,y=x,fill=factor(genotype)),position="fill",col=NA)+
    scale_fill_manual(values=c('dodgerblue','orange','grey'),name='genotype',labels=c('aa/bb','Aa/Bb','AA/BB','')) +   theme(panel.border=element_rect(colour="black",size=1,fill=NA))
  g <- ggplotGrob(p3)
  
  g$heights[c(9+2*as.numeric(sep_isle))] <- unit(3,'null')
  grid.newpage()
  grid.draw(g)
}


ldat <- restructure_df(out1$ts,out1$dum)
ldat2 <- restructure_df(out2$ts,out2$dum)
pdf("time-series-3.pdf",width=5,height=5)
stbar(ldat)
gr31(ldat,TRUE)
gr31(ldat)
lines_dist(ldat$m,ldat$f,ht=0.1,ylm=c(-0.1,1))
# 
stbar(ldat2)
gr31(ldat2,TRUE)
gr31(ldat2)
lines_dist(ldat2$m,ldat2$f,ht=0.05,ylm=c(0,0.4))

dev.off()

