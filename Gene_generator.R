filename = "growth-rate.gen"
N_genes = 10
N_alleles = 10
##### GENETICS
FunctionGvalues <- function(nbLoci=10,nbAlleles=10,dominance=0.5,SDeffects=1,SDalleles=1) {
  
  ## Initialising a matrix that will contain the genotypic effects on trait
  gvalues <- array(data=NA,dim=c(nbAlleles,nbAlleles,nbLoci),
                   dimnames=list(paste("A",1: nbAlleles,sep=""),
                                 paste("A",1: nbAlleles,sep=""),
                                 paste("L",1:nbLoci,sep=""))) 
  
  for (L in 1:nbLoci) {
    ## Setting the effects for the homozygotes [all loci]
    ## alter the locus importance in a realistic way (many small-effect loci, few major loci)
    effect <- abs(rnorm(n=1,mean=0,sd=SDeffects))
    diag(gvalues[,,L]) <- 2 * rnorm(n=dim(gvalues)[1],mean=0, sd=effect * SDalleles)
    ## Setting the effects for the heterozygotes
    ## loop for off-diagonal = heterozygotes (additive and dominance effects)
    for (A in 1:(nbAlleles - 1)) {
      for (D in (A + 1):nbAlleles) {
        d <- dominance * runif(n=1,min=-0.5,max=0.5)
        ## mean of additive effects + dominance, over diagonal
        gvalues[A,D,L] <- (0.5 - d) * gvalues[A,A,L] + (0.5 + d) * gvalues[D,D,L] 
        ## the same below diagonal    
        gvalues[D,A,L] <- (0.5 - d) * gvalues[A,A,L] + (0.5 + d) * gvalues[D,D,L] 
      }
    }
  }
  return(gvalues)
}
gen_phen_map <- FunctionGvalues(nbLoci=N_genes, nbAlleles=N_alleles, dominance=1,SDeffects=0.2, SDalleles=0.2)
gen_phen_map <- gen_phen_map / sum(gen_phen_map) * (N_genes)/(N_alleles^2) # for the factors
save(gen_phen_map,file=filename)
