rm(list=ls())
setwd("~/Documents/projects/Project_Fur_Seals/Matrix-Model/")
source("matrix_model.R")
source("graphing.R")

#library(profvis)
# profvis({
run_sim(filename = "test",Nt=1e4)
# })
plotting_rep(filename = "test","graph title")
