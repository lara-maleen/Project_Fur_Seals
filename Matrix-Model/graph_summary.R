setwd("~/Documents/projects/Project_Fur_Seals/Matrix-Model/out4-capping-symmadv/")


dat_raw <- read.csv("summary.stat")
ref_raw <- read.csv("simruns.csv")

library(ggplot2)


pdf("summary.pdf")
for(i in unique(ref_raw$maxfreq)){
# i <- unique(ref_raw$maxfreq)[1]
plot(0,0,xlim=c(-2,2),ylim=c(-2,2),type="n",axes=FALSE,xlab="",ylab="")
text(0,0,paste("maxfreq =",i))
dat <- dat_raw[ref_raw$maxfreq == i,]
ref <- ref_raw[ref_raw$maxfreq == i,]

plot_stab <- aggregate(dat$stability,by=list(surv=ref$surv,A.adv=ref$A.adv,dens_reg=ref$dens_reg), FUN = function(x) mean(x==0))

p1 <- ggplot(plot_stab,aes(x=surv,y=x)) + geom_point() + facet_grid(dens_reg ~ A.adv,labeller=label_both) + ggtitle("stability") + theme_bw()
print(p1)

byall <- function(x,fun,name=''){
  tmp <- aggregate(x,by=list(surv=ref$surv,A.adv=ref$A.adv,dens_reg=ref$dens_reg), fun)
  tmp$name <- name
  tmp
}

plot_var <- rbind(byall(dat$N.1.m, var,'N.1.m'), byall(dat$N.2.m, var,'N.2.m'), byall(dat$N.1.f, var,'N.1.f'), byall(dat$N.2.f, var,'N.2.f')) 

p2 <- ggplot(plot_var,aes(x=surv,y=x,col=name)) + geom_point() + facet_grid(dens_reg ~ A.adv,labeller=label_both) + ggtitle("variance") + theme_bw()
print(p2)
# plot_var
library(tidyr)
dat_all <- cbind(ref[,-1],dat[,-1])
plot_all <- gather(dat_all,condition, val, N.1.m:N.2.f)
plot_all$sex <- c('f','m')[1 + as.numeric(plot_all$condition %in% c("N.1.m","N.2.m"))] 
plot_all$loc <- factor(1 + as.numeric(plot_all$condition %in% c("N.2.f","N.2.m")))
p3 <- ggplot(plot_all,aes(x=surv,y=val,group=paste(condition,surv),col=sex,lty=loc)) + geom_boxplot() + facet_grid(dens_reg ~ A.adv,labeller=label_both) +theme_bw() +ggtitle("Number of individuals")
print(p3)

plot_allAA <- gather(dat_all,condition, val, N.1.AA.m:N.2.aa.m)
p4 <- ggplot(plot_allAA,aes(x=surv,y=val,group=paste(condition,surv),col=condition)) + geom_boxplot() + facet_grid(dens_reg ~ A.adv,labeller=label_both) +theme_bw() + ggtitle("Male numbers AA")
print(p4)

plot_allBB <- gather(dat_all,condition, val, N.1.BB.f:N.2.bb.f)
p5 <- ggplot(plot_allBB,aes(x=surv,y=val,group=paste(condition,surv),col=condition)) + geom_boxplot() + facet_grid(dens_reg ~ A.adv,labeller=label_both) +theme_bw() + ggtitle("Female numbers BB")
print(p5)


# pick some condition
aadv <- 1.5
surv <- 0.8
dat_all2 <- dat_all[dat_all$A.adv == aadv & dat_all$surv == surv,]

plot_all2 <- gather(dat_all2,condition, val, N.1.m:N.2.f)
plot_all2$sex <- c('f','m')[1 + as.numeric(plot_all2$condition %in% c("N.1.m","N.2.m"))] 
plot_all2$loc <- factor(1 + as.numeric(plot_all2$condition %in% c("N.2.f","N.2.m")))
p6 <- ggplot(plot_all2,aes(x=factor(dens_reg),y=val,group=paste(condition,dens_reg),col=sex,lty=loc)) + geom_boxplot() + facet_grid(min_val_m ~ min_val_f,labeller=label_both) +theme_bw() +ggtitle(paste("Number of individuals - surv = ",surv,", A.adv = ",aadv,sep=""))
print(p6)

plot_all3 <- dat_all[dat_all$A.adv == aadv & dat_all$surv == surv,]
plot_all3$relAA <- with(plot_all3,(N.1.AA.m + N.2.AA.m)/(N.1.m+N.2.m))
plot_all3$relBB <- with(plot_all3,(N.1.BB.f + N.2.BB.f)/(N.1.f+N.2.f))
p7 <- ggplot(plot_all3,aes(x=relAA,y=relBB)) + geom_smooth(data=subset(plot_all3, relAA<(1-1e-10)),method = "lm") + geom_point() + facet_grid(min_val_m ~ min_val_f,labeller=label_both) +theme_bw() 
print(p7)
}
dev.off()