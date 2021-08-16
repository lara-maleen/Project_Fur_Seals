# parameters
gamma_m = 0.05
gamma_f = 0.05
d = 0.5
d2 = 0.5
delta = 0.17
wm = 0.5
P1 = 1.8

p_surv <- function(N.m,N.f,N.m.tot,N.f.tot,wm){
  1-0.05*plogis(7.5*(wm*N.m/N.m.tot+(1-wm)*N.f/N.f.tot-0.5))
}

# Population
NmAA <- 0.25
NmAa <- 0.5
Nmaa <- 0.25
NfBB <- 0.25
NfBb <- 0.5
Nfbb <- 0.25

# Normalise
Ntot <- NmAA + NmAa + Nmaa + NfBB + NfBb + Nfbb
NmAA <- NmAA/Ntot
NmAa <- NmAa/Ntot
Nmaa <- Nmaa/Ntot
NfBB <- NfBB/Ntot
NfBb <- NfBb/Ntot
Nfbb <- Nfbb/Ntot


# Fitness calculations
# Part 1: distribution on the islands
N1mAA <- (1-gamma_m)*NmAA
N1mAa <- (gamma_m + d2*(1-2*gamma_m))*NmAa 
N1maa <- gamma_m*Nmaa
N1fBB <- (1-gamma_f)*NfBB
N1fBb <- (gamma_f+d2*(1-2*gamma_f))*NfBb
N1fbb <- gamma_f*Nfbb

N2mAA <- gamma_m*NmAA
N2mAa <- (1-gamma_m + d2*(2*gamma_m-1))*NmAa 
N2maa <- (1-gamma_m)*Nmaa
N2fBB <- gamma_f*NfBB
N2fBb <- (1-gamma_f+d2*(2*gamma_f-1))*NfBb
N2fbb <- (1-gamma_f)*Nfbb

N1f <- N1fBB+N1fBb+N1fbb
N2f <- N2fBB+N2fBb+N2fbb
N1m <- N1mAA+N1mAa+N1maa
N2m <- N2mAA+N2mAa+N2maa
# Part 2: Total pups on each island
DN1 <- N1f*p_surv(N1m,N1f,N1m+N2m,N1f+N2f,wm)
DN2 <- N2f*p_surv(N2m,N2f,N1m+N2m,N1f+N2f,wm)

# Part 3: Pups per male genotype
FAA <- P1*N1mAA*DN1/(P1*N1mAA+P1^d*N1mAa+N1maa) + P1*N2mAA*DN2/(P1*N2mAA+P1^d*N2mAa+N2maa)
FAa <- P1^d*N1mAa*DN1/(P1*N1mAA+P1^d*N1mAa+N1maa) + P1^d*N2mAa*DN2/(P1*N2mAA+P1^d*N2mAa+N2maa)
Faa <- N1maa*DN1/(P1*N1mAA+P1^d*N1mAa+N1maa) + N2maa*DN2/(P1*N2mAA+P1^d*N2mAa+N2maa)
if(abs(FAA+FAa+Faa-DN1-DN2)>1e-12){
  stop("Not all offspring attributed to a father")
}
# Part 4: Total per-time step fitness (incl. survival)
wAA <- (FAA/2+(1-delta)*NmAA)/(NmAA*(DN1+DN2+N1f+N2f+(1-delta)*NmAA+(1-delta)^d*NmAa+Nmaa))
wAa <- (FAa/2+(1-delta)^d*NmAa)/(NmAa*(DN1+DN2+N1f+N2f+(1-delta)*NmAA+(1-delta)^d*NmAa+Nmaa))
waa <- (Faa/2+Nmaa)/(Nmaa*(DN1+DN2+N1f+N2f+(1-delta)*NmAA+(1-delta)^d*NmAa+Nmaa))
(waa*Nmaa + wAa*NmAa + wAA*NmAA)/(Nmaa+NmAa+NmAA)
