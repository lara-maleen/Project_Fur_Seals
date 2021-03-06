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

FunctionGvaluesHetAd <- function(nbLoci=10,nbAlleles=10,min_phen=0,max_phen=1) {
  
  # width <- max_phen - min_phen
  # width_per_locus <- width/nbLoci
  # allele_vals <- seq(0,width_per_locus,length.out = nbAlleles) + min_phen/nbLoci
  #g_one_locus <- outer(allele_vals,allele_vals,FUN = function(x,y){ (x+y)/2})
  g_one_locus <- matrix(max_phen/nbLoci,ncol=nbAlleles,nrow=nbAlleles)
  diag(g_one_locus) <-min_phen/nbLoci
  gvalues <- array(rep(g_one_locus,nbLoci),dim = c(nbAlleles,nbAlleles,nbLoci))
  ## Initialising a matrix that will contain the genotypic effects on trait
  dimnames(gvalues)<- list(paste("A",1: nbAlleles,sep=""),
                           paste("A",1: nbAlleles,sep=""),
                           paste("L",1:nbLoci,sep="")) 
  return(gvalues)
}

gen_phen_map <- FunctionGvalues(10,10,min_phen=1,max_phen=50)
saveRDS(gen_phen_map, file="genes.rds")

gen_phen_map <- FunctionGvalues(10,10,min_phen=-2,max_phen=2)
saveRDS(gen_phen_map, file="genes2.rds")

gen_phen_map <- FunctionGvaluesHetAd(10,10,min_phen=1,max_phen=50)
saveRDS(gen_phen_map, file="genes_het_ad.rds")

gen_phen_map <- FunctionGvaluesHetAd(10,10,min_phen=-2,max_phen=2)
saveRDS(gen_phen_map, file="genes2_het_ad.rds")
