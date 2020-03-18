setwd("~/Documents/projects/Project_Fur_Seals/Matrix-Model/out3")


dat <- read.csv("summary.stat")
ref <- read.csv("simruns.csv")

library(ggplot2)

head(dat)
head(ref)

pdf("summary.pdf")
plot_stab <- aggregate(dat$stability,by=list(surv=ref$surv,A.adv=ref$A.adv,dens_reg=ref$dens_reg), FUN = function(x) mean(x==0))

p1 <- ggplot(plot_stab,aes(x=surv,y=x)) + geom_point() + facet_grid(dens_reg ~ A.adv,labeller=label_both) + ggtitle("stability") + theme_bw()
print(p1)

byall <- function(x,fun){
  aggregate(x,by=list(surv=ref$surv,A.adv=ref$A.adv,dens_reg=ref$dens_reg), fun)
}

plot_var <- rbind(byall(dat$N.1.m, var), byall(dat$N.2.m, var), byall(dat$N.1.f, var), byall(dat$N.2.f, var)) 

p2 <- ggplot(plot_var,aes(x=surv,y=x)) + geom_point() + facet_grid(dens_reg ~ A.adv,labeller=label_both) + ggtitle("variance") + theme_bw()
print(p2)

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

dev.off()