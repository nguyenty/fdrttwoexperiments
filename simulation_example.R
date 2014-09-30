
source("normal_sim_pvalues.R")
source("micro_sim_pvalues.R")
source("qvalue_functions.R")


##Different combinations of the following values were used to create
##the different simulation settings.
  ##ni = 4, 10, 20
  ##mv = c(3000,2000,2000,3000), c(5000, 1500, 1500, 2000), 
  ##     c(7000,1000,1000,1000), c(9000,250,250,500)
  ##mnd = 1, 2

############################################################################
############################################################################
##Example of normal data simulation and analysis for one simulation setting
##with nrep=100.
############################################################################
############################################################################


  nis <- 10
  mvs <- c(3000,2000,2000,3000); ms <- sum(mvs)
  mnds <- 1
  nreps <- 100

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
    seFDRints <- sd(FDRints)/sqrt(nrep)







  
 




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





