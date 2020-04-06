setwd("~/Documents/projects/Project_Fur_Seals/Kirkpatrick/")

t2s <- seq(0,1,0.01)
plot(t2s,plogis(10*(t2s - 0.5)),type="l")

p2 <- function(t2,a2,sfun){
  -sfun(t2)*(1+t2*(a2-1))/((sfun(t2)-1-sfun(t2)*t2)*(a2-1))
}

p2_f <- function(t2,a2,sfun){
  (1-sfun(t2)*t2)/(2-(sfun(t2)*a2*t2)/(1-t2+a2*t2) -sfun(t2)*t2) 
}

s_lin <- function(t2){
  0.5*t2
}

s_log <- function(t2){
  0.5*plogis(10*(t2-0.5))
}

df_lin <- expand.grid(t2=seq(0,1,0.01),a2 = c(1.1,1.5,2,5))
df_lin$fun <- "linear"
df_plog <- df_lin
df_plog$fun <- 'inverse logit'
df_lin_f <- df_lin
df_plog_f <- df_plog
df_lin$sex = "male"
df_plog$sex = "male"
df_lin_f$sex = "female"
df_plog_f$sex = "female"

df_lin$p2 <- p2(df_lin$t2,df_lin$a2,s_lin)
df_plog$p2 <- p2(df_plog$t2,df_plog$a2,s_log)
df_lin_f$p2 <- p2_f(df_lin_f$t2,df_lin_f$a2,s_lin)
df_plog_f$p2 <- p2_f(df_plog_f$t2,df_plog_f$a2,s_log)

df <- rbind(df_lin,df_plog,df_lin_f,df_plog_f)
# df$p2[df$p2 > 1] <- 1
library(ggplot2)
#df$as2 <- factor(as.character(df$a2))

g1 <- ggplot(df[df$sex=='male',],aes(x=t2,y=p2,col=factor(a2))) + geom_line() + theme_bw() + facet_grid(.~fun) + ylim(c(0,1)) +  scale_colour_discrete("a2")
g2 <- ggplot(df,aes(x=t2,y=p2,col=factor(a2),lty=sex)) + geom_line() + theme_bw() + facet_grid(.~fun) + ylim(c(0,1)) +  scale_colour_discrete("a2")

pdf("t2p2-males.pdf",height=3,width=5)
print(g1)
print(g2)
dev.off()
