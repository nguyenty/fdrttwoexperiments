library(plyr)
source("normal_sim_pvalues.R")
source("micro_sim_pvalues.R")
source("qvalue_functions.R")


##Different combinations of the following values were used to create
##the different simulation settings.
  ##ni = 4, 10, 20
  ##mv = c(3000,2000,2000,3000), c(5000, 1500, 1500, 2000), ## m00, m01, m10, m11
  ##     c(7000,1000,1000,1000), c(9000,250,250,500)
  ##mnd = 1, 2

############################################################################
############################################################################
##Example of normal data simulation and analysis for one simulation setting
##with nrep=100.
############################################################################
############################################################################


  nis <- 10
  mvs <- c(4500,1,1,500); ms <- sum(mvs)
  mnds <- 1
  nreps <- 50

  ##calculate p-values for each pair of simulated experiments
  ps <- normal_sim_pvalues(ni=nis, d0=3.637578, s20=0.04040928, 
        mv=mvs, mnd=mnds, nrep=nreps)

  ##calculate maximum q-values (q_max) for each pair of experiment
    qmaxs <- apply(ps, 2, function(x) qvalue_maximum(pv1=x[1:ms], pv2=x[(ms+1):(2*ms)]))

    
  ##calculate proposed q-values (q_I) for each pair of experiments. (This takes a while.
  ##Between 15 and 20 minutes on my laptop).
    qints <- apply(ps, 2, function(x) qvalue_intersection(pv1=x[1:ms], pv2=x[(ms+1):(2*ms)]))


  ##Calculate V (number of false discoveries), R (number of discoveries), and
  ##S (number of true discoveries) and observed FDR for each data set for each method

    m11s <- mvs[4]

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
    seFDRints <- sd(FDRints)/sqrt(nreps)

## Phillip and Gosh paper #####

library(deldir)
library(MASS)

#-------------------------------------------------------------------------------------#
#     PART 1. Functions needed for Voronoi P-value Combination
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# CombineVoronoi: Function to perform P-value combination for 2-d vectors of p-values
# INPUTS: x.values (vector) = first set of p-values 
#         y.values (vector) = second set of p-values
#         ranking (value of 1,2,3,4) = indication of which ordering scheme to use.
#                 1 -> Euclidean
#                 2 -> Maximum
#                 3 -> Summation
#                 4 -> De Lichtenberg (not reccommended for data analysis)
# OUTPUTS: index (vector) = indices of ranked p-values.  i.e. the first value of 'index' 
#           gives the index of the first ranked gene in the original vectors of p-values
#          cum.areas (vector) = ordered cumulative areas calculated using specified 
#           rankind scheme 
#-------------------------------------------------------------------------------------#

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



voronoi <- apply(ps, 2, function(x) CombineVoronoi(x.values=x[1:ms], y.values=x[(ms+1):(2*ms)]))

VSRvoronoi <- laply(1:nreps, function(i) {
  d <- as.data.frame(cbind(voronoi[[i]][[1]], voronoi[[i]][[2]]))
  dd <- (arrange(d, V1)) # rearrange the gene as the original order 1, 2, 3, ,,,
#   head(dd) 
#   d[which(d[,2] == dd[1,2]),]
#   d[1,]
# dim(d)
  bh <- BH(dd[,2])
  R.t <- bh$k
  S.t <- sum(bh$index >sum(mvs)- mvs[4])
  V.t <- R.t - S.t
  return(c(V.t = V.t, S.t = S.t, R.t = R.t))
})

VSRvoronoi
FDRvoronoi <- apply(VSRvoronoi[,c(1,3)], 1, function(x) x[1]/max(x[2],1))

##mean number of true discoveries (mean S)
meanFDRvoronoi <- mean(FDRvoronoi, na.rm = T)

seFDRvoronoi <- sd(FDRvoronoi, na.rm = T)/sqrt(nreps - sum(is.na(FDRvoronoi)))

##power of the procedure #####
powerFDRvoronoi <- VSRvoronoi[,2]
mean(powerFDRvoronoi)



############################################################################
############################################################################
##Example of microarray data simulation and analysis for one simulation setting
##with nrep=100.
############################################################################
############################################################################

  y1 <- read.table("micro_GSE12417_samp1.txt")
  y2 <- read.table("micro_GSE12417_samp2.txt")

  nis <- 10
  mvs <- c(3000,2000,2000,3000); ms <- sum(mvs)
  mnds <- 1
  nreps <- 100

  ##calculate p-values for each pair of simulated experiments
  ps <- micro_sim_pvalues(y1=y1, y2=y2, ni=nis, mv=mvs, mnd=mnds, nrep=nreps)

  ##calculate maximum q-values (q_max) for each pair of experiment
    qmaxs <- apply(ps, 2, function(x) qvalue_maximum(pv1=x[1:ms], pv2=x[(ms+1):(2*ms)]))
    
  ##calculate proposed q-values (q_I) for each pair of experiments. (This takes a while.
  ##About 25 minutes on my laptop).
    qints <- apply(ps, 2, function(x) qvalue_intersection(pv1=x[1:ms], pv2=x[(ms+1):(2*ms)]))

  ##Calculate V (number of false discoveries), R (number of discoveries), and
  ##S (number of true discoveries) and observed FDR for each data set for each method

    m11s <- mvs[4]

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
    seFDRints <- sd(FDRints)/sqrt(nreps)





