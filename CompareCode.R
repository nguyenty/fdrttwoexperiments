## Compare performance of two methods of FDR control (Gosh et.al) and FDR estimate (Megan et.al)####
## Gosh simulation set up ######
## independent the same mu_i, variance = 1, 2000 genes. 


source("CombineVoronoi_BH.R")
source("estimate_m0s_fn.R")
source("qvalue_functions.R")
source("SimulationPvalue_Daisy.R")
source("MeganDaisyRes.R")
source("normal_sim_pvalues.R")


simdaisy <- function(n,pi0,mu_a, nrep){
  ps <- daisy.pvalue(n, pi0, mu_a, nrep) 
  intsout <- daisy_ints_out(ps,pi0)
  voronoiout <- daisy_voronoi_out(ps, pi0)
  output <- rbind(Megan= intsout, Daisy = voronoiout)
  print(paste0("n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", nrep))
  path <- paste0("phillips_","n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", nrep, "_out.csv" )
  write.csv(output, file = path, row.names = F)
  output
}




n <- c(2000, 5000, 10000)
pi0 <- c(0.95, 0.9, 0.8, 0.7 )
mu_a <- c(1, 2, 3, 4)

for(i in 3){
  for(j in 4){
    for(k in 4)simdaisy(n[i], pi0[j], mu_a[k], 100)
  }
}


for(i in 1:3){
  for(j in 1:4){
    for(k in 1:4)simdaisy(n[i], pi0[j], mu_a[k], 100)
  }
}
#simvoronoi(n=5000, pi0=0.9, mu_a=4, nrep =4)


## Megan Simulation normal simulation data######


# nis <- 10; mvs <- c(3000,2000,2000,3000); mnds <- 2; nreps <- 5
simmegan <- function(nis, mvs, mnds, nreps){
  ps <- normal_sim_pvalues(ni=nis, d0=3.637578, s20=0.04040928, 
                           mv=mvs, mnd=mnds, nrep=nreps)
  intsout <- megan_ints_out(ps,mvs)
  voronoiout <- megan_voronoi_out(ps, mvs)
  output <- rbind(Megan= intsout, daisy = voronoiout)
  print(paste0("Megan_", "ni_",nis,"m11_", mvs[4],  "mnds_", mnds, "nrep_", nreps))
  path <- paste0("Megan_","ni_",nis, "m11_", mvs[4], "mnds_", mnds, "nrep_", nreps, "_out.csv" )
  write.csv(output, file = path, row.names = F)
  output
}


nis <- c(4, 10, 20)
mvs <- list(c(9000, 250, 250, 500), 
            c(7000, 1000, 1000, 1000), 
            c(5000, 1500, 1500, 2000), 
            c(3000,2000,2000,3000))
mnds <- c(1, 2)
nreps <- 100

for(i in 2:3){
  for(j in 1:4){
    for (k in 1:2) simmegan(nis[i], mvs[[j]], mnds[k], nreps)
  }
}