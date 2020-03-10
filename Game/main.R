rm(list=ls())

library(network)
setwd("~/Documents/projects/Project_Fur_Seals/Game/")
source("fnc.R")
# dimensions of the world
xsize <- 160
ysize <- 160

# following objects can be stored and reloaded, using save and load (especially create_world is time consuming)
world <- create_world(xsize,ysize,maxsurf=3.4e3)

# init_pop also takes arguments x and y to directly specify the x and y coordinates of the pups
init_pop <- create_init_pop(world)

# to store output to pdf
# Nt: number of time points
# filename: output filename
run_game_pdf(init_pop,xsize,ysize,world,Nt=100,filename="test.pdf")

# to show output live on screen, use pause (>=0) to determine the time between two consecutive graphs
run_game_live(init_pop,xsize,ysize,world,Nt=100,pause=0.1)
