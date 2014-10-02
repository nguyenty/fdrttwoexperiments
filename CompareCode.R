## Compare performance of two methods of FDR control (Gosh et.al) and FDR estimate (Megan et.al)####
## Gosh simulation set up ######
## independent the same mu_i, variance = 1, 2000 genes. 
library(MASS)
library(plyr)
library(qvalue)
source("CombineVoronoi_BH.R")
source("estimate_m0s_fn.R")
source("qvalue_functions.R")
source("SimulationPvalue_Phillips.R")


## Simulation in Phillips and Gosh paper ####
##calculate p-values for each pair of simulated experiments
n <- 5000
pi0 <- 0.9
mu_a <- 4
n.rep <- 50
ps <- sim.pvalue(n, pi0, mu_a, n.rep)

##calculate proposed q-values (q_I) for each pair of experiments. (This takes a while.
##Between 15 and 20 minutes on my laptop).
pm1 <- proc.time()
qints <- apply(ps, 2, function(x) qvalue_intersection(pv1=x[1:n], pv2=x[(n+1):(2*n)]))
proc.time() -pm1

##Calculate V (number of false discoveries), R (number of discoveries), and
##S (number of true discoveries) and observed FDR for each data set for each method



#####################
##proposed q-value####
#####################

Vints <- apply(qints[de+ 1:(n-de),], 2, function(x) sum(x <= 0.05))
Rints <- apply(qints, 2, function(x) sum(x <= 0.05))
Sints <- Rints - Vints
FDRints <- apply(cbind(Vints, Rints), 1, function(x) x[1]/max(x[2],1))

##mean number of true discoveries (mean S)
meanSints <- mean(Sints)
seSints <- sd(Sints)/sqrt(n.rep)
meanSints
##mean observed FDR
meanFDRints <- mean(FDRints)
seFDRints <- sd(FDRints)/sqrt(n.rep)
meanFDRints


## Calculate voronoi area for each genes in each simulation data set
pm1 <- proc.time()
voronoi <- apply(ps, 2, function(x) CombineVoronoi(x.values=x[1:n], y.values=x[(n+1):(2*n)]))
proc.time() -pm1

VSRvoronoi <- laply(1:n.rep, function(i) {
  d <- as.data.frame(cbind(voronoi[[i]][[1]], voronoi[[i]][[2]]))
  dd <- (arrange(d, V1)) # rearrange the gene as the original order 1, 2, 3, ,,,
  
  bh <- BH(dd[,2])
  R.t <- bh$k
  S.t <- sum(bh$index <=de)
  V.t <- R.t - S.t
  return(c(V.t = V.t, S.t = S.t, R.t = R.t))
  })


FDRvoronoi <- apply(VSRvoronoi[,c(1,3)], 1, function(x) x[1]/max(x[2],1))

##mean number of true discoveries (mean S)
meanFDRvoronoi <- mean(FDRvoronoi)
seFDRvoronoi <- sd(FDRvoronoi)/sqrt(n.rep)
meanFDRvoronoi
##power of the procedure #####
powerFDRvoronoi <- VSRvoronoi[,2]
mean(powerFDRvoronoi)
meanFDRints
meanSints
