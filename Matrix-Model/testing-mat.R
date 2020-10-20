# test matrix model functions
source("~/Documents/projects/Project_Fur_Seals/Matrix-Model/matrix_model_kirk_version.R")

# generating a dummy object:
dum <- run_sim(min_val_m = 0,min_val_f = 0,dumgen=TRUE)

# testing male_dist:
# simple test
d1 <- male.dist(dum,c(0,0,1,rep(0,15)),maxfreq=0.7,normalize = FALSE)

lapply(d1,sum)
sum(d1[[1]])
sum(d1[[2]])

# different pattern

# maxfreq

# testing offs_dist

# testing make_mat

# testing run_sim