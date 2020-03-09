surv <- function(N.f.1,dens_reg){
  1-dens_reg+dens_reg*plogis(5*(0.25-N.f.1))
}

xs <- seq(0,2,0.025)
alp <- 0.6
setwd("~/Documents/projects/Project_Fur_Seals/Images")
library(tikzDevice)
options(tikzLatexPackages = c(getOption("tikzLatexPackages"), "\\usepackage{amsmath}\n \\usepackage{xcolor}")) 
tikz("illustration.tex",standAlone = TRUE,width=4,height=4)
plot(xs,surv(xs,alp),xlim=c(0,2),ylab="$p_{\\text{surv},i}$",xlab="$N_{f,i}$",ylim=c(0,1),type="l",axes=FALSE)
axis(1)
tks <- seq(0,1,0.2)
lbl <- tks
ind <- which.min(abs(tks-1+alp))
tks[ind] <- 1-alp
lbl[ind] <- "$\\color{red} 1-\\alpha$"
axis(2,labels = lbl,at = tks,las=1)
abline(h=1-alp,lty=2)
dev.off()
system("pdflatex illustration.tex")
