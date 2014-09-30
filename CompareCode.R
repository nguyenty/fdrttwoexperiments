## Compare performance of two methods of FDR control (Gosh et.al) and FDR estimate (Megan et.al)####
## Gosh simulation set up ######
## independent the same mu_i, variance = 1, 2000 genes. 
library(MASS)
set.seed(1)
n = 5000
pi0 = .9
ee = round(n*pi0)
de = n - ee
mu_a =c(2,3,4) 
n.rep <- 100
sim.pvalue <- matrix(0,  nrow = 2*n, ncol = n.rep)
for(i in 1:n.rep){
  t.value <- mvrnorm(n = n, mu = c(0,0), Sigma = diag(1,2))
  t.value[1:de, ] <-mvrnorm(n = de, mu = c(2,2), Sigma = diag(1, 2)) 
  p.value  <- apply(t.value, c(1, 2), function(x) {2*(1 - pnorm(abs(x)))})
  sim.pvalue[,i] <- c(p.value[,1], p.value[,2])  
}

#hist(sim.pvalue[1:n,1], nclass = 30)

# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")
library(qvalue)
source("qvalue_functions.R")
source("estimate_m0s_fn.R")
#source("simulation_example.R")


##calculate p-values for each pair of simulated experiments
ps <- sim.pvalue

##calculate maximum q-values (q_max) for each pair of experiment
ms <- n

qmaxs <- apply(ps, 2, function(x) qvalue_maximum(pv1=x[1:ms], pv2=x[(ms+1):(2*ms)]))


##calculate proposed q-values (q_I) for each pair of experiments. (This takes a while.
##Between 15 and 20 minutes on my laptop).
qints <- apply(ps, 2, function(x) qvalue_intersection(pv1=x[1:ms], pv2=x[(ms+1):(2*ms)]))


##Calculate V (number of false discoveries), R (number of discoveries), and
##S (number of true discoveries) and observed FDR for each data set for each method

m11s <- de

#####################
##max q-value
#####################

Vmaxs <- apply(qmaxs[1:(ms-m11s),], 2, function(x) sum(x <= 0.05))
Rmaxs <- apply(qmaxs, 2, function(x) sum(x <= 0.05))
Smaxs <- Rmaxs - Vmaxs
FDRmaxs <- apply(cbind(Vmaxs, Rmaxs), 1, function(x) x[1]/max(x[2],1))

##mean number of true discoveries (mean S)
meanSmaxs <- mean(Smaxs)
seSmaxs <- sd(Smaxs)/sqrt(nreps)

##mean observed FDR
meanFDRmaxs <- mean(FDRmaxs) 
seFDRmaxs <- sd(FDRmaxs)/sqrt(nreps)

#####################
##proposed q-value
#####################

Vints <- apply(qints[1:(ms-m11s),], 2, function(x) sum(x <= 0.05))
Rints <- apply(qints, 2, function(x) sum(x <= 0.05))
Sints <- Rints - Vints
FDRints <- apply(cbind(Vints, Rints), 1, function(x) x[1]/max(x[2],1))

##mean number of true discoveries (mean S)
meanSints <- mean(Sints)
seSints <- sd(Sints)/sqrt(nreps)

##mean observed FDR
meanFDRints <- mean(FDRints)
seFDRints <- sd(FDRints)/sqrt(nrep)













########