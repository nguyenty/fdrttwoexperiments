## Compare performance of two methods of FDR control (Gosh et.al) and FDR estimate (Megan et.al)####
## Gosh simulation set up ######
## independent the same mu_i, variance = 1, 2000 genes. 


source("CombineVoronoi_BH.R")
source("estimate_m0s_fn.R")
source("qvalue_functions.R")
source("SimulationPvalue_Daisy.R")
source("SimulationPvalue_Phillips.R")
source("MeganDaisyRes.R")
source("normal_sim_pvalues.R")


simdaisy <- function(n,pi0,mu_a, nrep){
  ps <- daisy.pvalue(n, pi0, mu_a, nrep) 
  intsout <- daisy_ints_out(ps,pi0)
  voronoiout <- daisy_voronoi_out(ps, pi0)
  output <- rbind(Megan= intsout, Daisy = voronoiout)
  row.names(output) <- c(paste0("n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", nrep, "_Megan"), 
                         paste0("n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", nrep, "_Daisy"))
  print(paste0("n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", nrep))
  path <- paste0("phillips_","n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", nrep, "_out.csv" )
  #write.csv(output, file = path, row.names = F)
  output
}




n <- c(2000)
pi0 <- c(0.9,  0.8 )
mu_a <- c(1, 2, 4)
nrep <- 50
df <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 

for(i in 1){
  for(j in 1:2){
    for(k in 1:3) {
      dd <- simdaisy(n[i], pi0[j], mu_a[k], nrep) # n<- 500; i <- 1; pi0 <- 0.95; j <- 3; mu_a <- 4; k <- 4
    df <- rbind(dd, df)}
  }
}

write.csv(df, file ="DaisySim.csv", row.names = T)
#simvoronoi(n=5000, pi0=0.9, mu_a=4, nrep =4)




#simvoronoi(n=5000, pi0=0.9, mu_a=4, nrep =4)

## Megan Simulation normal simulation data######


# nis <- 10; mvs <- c(3000,2000,2000,3000); mnds <- 2; nreps <- 5
simmegan <- function(nis, mvs, mnds, nreps){
  ps <- normal_sim_pvalues(ni=nis, d0=3.637578, s20=0.04040928, 
                           mv=mvs, mnd=mnds, nrep=nreps)
  intsout <- megan_ints_out(ps,mvs)
  voronoiout <- megan_voronoi_out(ps, mvs)
  output <- rbind(Megan= intsout, daisy = voronoiout)
  row.names(output) <- c(paste0("Megan_", "ni_",nis,"m11_", mvs[4],  "mnds_", mnds, "nrep_", nreps, "_Megan"), 
                         paste0("Megan_", "ni_",nis,"m11_", mvs[4],  "mnds_", mnds, "nrep_", nreps, "_Daisy"))
  print(paste0("Megan_", "ni_",nis,"m11_", mvs[4],  "mnds_", mnds, "nrep_", nreps))
  path <- paste0("Megan_","ni_",nis, "m11_", mvs[4], "mnds_", mnds, "nrep_", nreps, "_out.csv" )
  #write.csv(output, file = path, row.names = F)
  output
}


nis <- c(4, 10)
mvs <- list(c(9000, 250, 250, 500)/5, 
            c(5000, 1500, 1500, 2000)/5)
mnds <- c(1, 2)
nreps <- 50

for(i in 1:2){ # i <- 2
  for(j in 1:2){
    for (k in 1:2) {
      dd <- simmegan(nis[i], mvs[[j]], mnds[k], nreps)
    df <- rbind(dd, df)}
  }
}

df

write.csv(df, file ="MeganSim.csv", row.names = T)
daisy_sim <- read.csv("DaisySim.csv")
?read.csv
row.names(daisy_sim) <- daisy_sim[,1]
megandaisy_sim <- rbind(df, daisy_sim[,-1])
write.csv(megandaisy_sim, "MeganDaisySim.csv", row.names = T)




## simmegan2


simmegan2 <- function(mvs,mu_a,nrep){
  ps <- megan.pvalue2(mvs, mu_a, nrep) 
  intsout <- megan_ints_out(ps,mvs)
  voronoiout <- megan_voronoi_out(ps, mvs)
  output <- rbind(Megan= intsout, daisy = voronoiout)
  row.names(output) <- c(paste0("Megan2_", "m", sum(mvs), "m11_", mvs[4], "mua", mu_a, "nrep_", nrep, "_Megan"), 
                         paste0("Megan2_", "m", sum(mvs),"m11_", mvs[4],"mua", mu_a,  "nrep_", nrep, "_Daisy"))
  print(paste0("Megan_","m11_", mvs[4],"mua_", mu_a, "nrep_", nrep))
  path <- paste0("Megan_", "m11_", mvs[4],"mua_", mu_a, "nrep_", nrep, "_out.csv" )
  #write.csv(output, file = path, row.names = F)
  output
}

mvs <- list(c(1000, 300, 300, 400),
            c(1600, 100, 100, 200)
            )
mu_a <- c(1,2,4)
nrep <- 50

df <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 
for(i in 1:length(mvs)){ # i <- 2
  for(j in 1:length(mu_a)){
    dd <- simmegan2(mvs[[i]],mu_a[j], nrep )
    df <- rbind(dd, df)}
}
}

df

write.csv(df, file ="Megan2Sim.csv", row.names = T)
