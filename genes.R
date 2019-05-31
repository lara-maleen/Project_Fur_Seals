##### GENETICS
FunctionGvalues <- function(nbLoci=10,nbAlleles=10,min_phen=0,max_phen=1) {
  
  width <- max_phen - min_phen
  width_per_locus <- width/nbLoci
  allele_vals <- seq(0,width_per_locus,length.out = nbAlleles) + min_phen/nbLoci
  g_one_locus <- outer(allele_vals,allele_vals,FUN = function(x,y){ (x+y)/2})
  
  gvalues <- array(rep(g_one_locus,nbLoci),dim = c(nbAlleles,nbAlleles,nbLoci))
  ## Initialising a matrix that will contain the genotypic effects on trait
  dimnames(gvalues)<- list(paste("A",1: nbAlleles,sep=""),
                           paste("A",1: nbAlleles,sep=""),
                           paste("L",1:nbLoci,sep="")) 
  return(gvalues)
}

gen_phen_map <- FunctionGvalues(10,10,min_phen=-3,max_phen=1)
saveRDS(gen_phen_map, file="genes.rds")
