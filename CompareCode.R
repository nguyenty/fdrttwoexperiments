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
##proposed q-value####
#####################

Vints <- apply(qints[1:(ms-m11s),], 2, function(x) sum(x <= 0.05))
Rints <- apply(qints, 2, function(x) sum(x <= 0.05))
Sints <- Rints - Vints
FDRints <- apply(cbind(Vints, Rints), 1, function(x) x[1]/max(x[2],1))

##mean number of true discoveries (mean S)
meanSints <- mean(Sints)
seSints <- sd(Sints)/sqrt(n.rep)

##mean observed FDR
meanFDRints <- mean(FDRints)
seFDRints <- sd(FDRints)/sqrt(n.rep)


CombineVoronoi = function(x.values,y.values,ranking=1){
  
  my.data = cbind(x.values,y.values)
  #get voronoi tessellation and extract cell areas for each p-vector
  areas = deldir(x.values,y.values,rw=c(0,1,0,1),digits=20,eps=1e-13) #Get Voronoi tessellation
  tess.areas = areas$summary$dir.area #extract cell areas 
  tess.areas[tess.areas<0]= 0 #In rare cases, deldir gives negative cell areas.  If this occurs, look more closely at inputs 
  
  #get distances according to ranking scheme
  if(ranking==1){distance = apply(my.data,1,function(x){sqrt(x[1]^2+x[2]^2)})} #1: Euclidean
  if(ranking==2){distance = apply(my.data,1,max)} #2: maximum
  if(ranking==3){distance = apply(my.data,1,sum)} #3: sum
  if(ranking==4){distance = apply(my.data,1,function(x){prod(x)*(1+(x[1]/.001)^2)*(1+(x[2]/.001)^2)})} #4: de lichtenberg
  
  mrank = sort(distance,index.return=T)$ix #rank p-vectors according to the distances
  csum = cumsum(tess.areas[mrank]) #find cumulative sums
  
  return(list(index=mrank,cumareas=csum)) 
}

#-------------------------------------------------------------------------------------#
# BH: Function that performs the B-H procedure on a single set of values.
# INPUTS: p.values (vector) = values between 0 and 1 on which to perform B-H procedure
#         alpha (value) = nominal level of FDR control
# OUTPUTS: k (value) = number of hypothesess to reject
#          index (Vector) = indices of rejected hypothesis in relation to original 
#                           vector of p-values
#-------------------------------------------------------------------------------------#
BH = function(p.values, alpha=.05)
{
  n = length(p.values)
  temp = sort(p.values,index.return=T)
  p.vals = temp$x
  k = which(p.vals/(1:n) <= (alpha/n))
  index = 0
  if(length(k)==0){k=0}
  if(length(k)>0){k = max(k);index = temp$ix[1:max(k)] }
  return(list(k=k,index=index))
}
## Calculate voronoi area for each genes in each simulation data set
voronoi <- apply(ps, 2, function(x) CombineVoronoi(x.values=x[1:ms], y.values=x[(ms+1):(2*ms)]))
str(voronoi[[100]])
library(plyr)
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

##power of the procedure #####
powerFDRvoronoi <- VSRvoronoi[,2]/de
mean(powerFDRvoronoi)
