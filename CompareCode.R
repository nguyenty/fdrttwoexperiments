## Compare performance of two methods of FDR control (Gosh et.al) and FDR estimate (Megan et.al)####
## Gosh simulation set up ######
## independent the same mu_i, variance = 1, 2000 genes. 

library(plyr)
source("CombineVoronoi_BH.R")
source("estimate_m0s_fn.R")
source("qvalue_functions.R")
source("SimulationPvalue_Phillips.R")
source("MeganPhillipsRes.R")



simvoronoi <- function(n,pi0,mu_a, n.rep){
  ps <- sim.pvalue(n, pi0, mu_a, n.rep) 
  intsout <- ints_out(ps,pi0)
  voronoiout <- voronoi_out(ps, pi0)
  output <- rbind(Megan= intsout, Phillips = voronoiout)
  print(paste0("n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", n.rep))
  output
}




n <- c(2000, 5000, 10000)
pi0 <- c(0.95, 0.8, 0.6 )
mu_a <- c(1, 2, 4)
#simvoronoi(n=1000, pi0=0.9, mu_a=1, n.rep =20)