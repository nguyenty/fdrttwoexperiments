library(limma)
library(pscl)


###################################################################
###################################################################
###################################################################
###################################################################

##The function "normal_sim_pvalues" simulates data from independent 
##normal distributions and calculates limma p-values.  The output
##of this function is a matrix of p-values with 2*m rows and nrep 
##columns. The columns correspond to p-values from a pair of
##experiments.  The first m values correspond to the p-values from
##the first experiment and the second m p-values correspond to
##p-values from the second experiment.  The inputs for this function
##are as follows:

  ##ni: the sample size in each treatment group
  ##d0: the "prior" degrees of freedom.  The value of d0 = 3.637578
    ##was used in the simulations.  This value was estimated from
    ##data described in Hannenhalli et al. (2006)
  ##s20: the "prior" variance.  This value was also estimated from
    ##the data described in Hannenhalli et al. (2006).  The inputs
    ##d0 and s20 define the inverse gamma distribution the gene
    ##variances are drawn from.
  ##mv: a vector with the following values: (m00, m01, m10, m11),
    ##where m00 is the number of genes that are EE in both
    ##experiments, m01 is the number of genes that are EE in Exp1
    ##but DE in Exp2, m10 is the number of genes that are DE in
    ##Exp1 but EE in Exp2, and m11 is the number of genes that are
    ##DE in both experiments.  NOTE: m = sum(mv).
  ##mnd: the (absolute) mean of the distribution of relative effect
    ##sizes
  ##nrep: the number of pairs of data sets to simulate.

normal_sim_pvalues <- function(ni, d0, s20, mv, mnd, nrep){
  pmat <- NULL

  m00 <- mv[1]
  m01 <- mv[2]
  m10 <- mv[3]
  m11 <- mv[4]

  m <- sum(mv)
  print(paste("m =",m))

  ##treatment indicators for each data set
  trt1 <- as.factor(c(sapply(c(1,2), rep, times=ni)))
  trt2 <- as.factor(c(sapply(c(1,2), rep, times=ni)))

  ##design matrix for each data set
  design1 <- model.matrix(~trt1+0)
  design2 <- model.matrix(~trt2+0)

  ##contrasts to be tested for each data set (two-sided t-test)
  colnames(design1)=c("t1","t2")
  colnames(design2)=c("t1","t2")
  contr.mat1 <- makeContrasts(t2-t1,levels=design1)
  contr.mat2 <- makeContrasts(t2-t1,levels=design2)

  ##Simulate data set and calculate p-values nrep times
  for(rep in 1:nrep){ 

   #########################################
   ##  Simulate data from first experiment## 
   #########################################

    ###########################################################
    ##Simulate data from genes that are EE in both experiments (both treatments
    ##have a mean of 0)

      ##Simulate standard deviations.
      s00.1 <- sqrt(rigamma(n=m00, alpha=d0/2, beta=d0*s20/2))
 
      ##Simulate data
      dat00.1 <- t(sapply(s00.1, rnorm, n=2*ni, mean=0))

    ###########################################################
    ##Simulate data from genes that are EE in Exp1 and DE in Exp2 (NOTE: this is
    ##data from Exp1, so all genes are EE and both treatment means are 0)

      ##Simulate standard deviations.
      s01.1 <- sqrt(rigamma(n=m01, alpha=d0/2, beta=d0*s20/2))
 
      ##Simulate data.
      dat011.1 <- t(sapply(s01.1, rnorm, n=ni, mean=0))
      dat012.1 <- t(sapply(s01.1, rnorm, n=ni, mean=0))
      dat01.1 <- cbind(dat011.1, dat012.1)

    ###########################################################
    ##Simulate data from genes that are DE in Exp1 and EE in Exp2 (NOTE: this is
    ##data from Exp1, so genes are DE)

      ##Simulate standard deviations.
      s10.1 <- sqrt(rigamma(n=m10, alpha=d0/2, beta=d0*s20/2))

      ##Simulate effect sizes.
      mn102.1 <- rnorm(n=m10, mean=mnd)*sample(c(rep(1,m10/2), rep(-1,m10/2)))*s10.1

      ##Simulate data.
      dat101.1 <- t(sapply(s10.1, rnorm, n=ni, mean=0))
      dat102n.1 <- t(sapply(s10.1, rnorm, n=ni, mean=0))
      dat102.1 <- sweep(dat102n.1, 1, mn102.1, "+")
      dat10.1 <- cbind(dat101.1, dat102.1)

    if(m11 == 0) {dat11.1 <- NULL}
    if(m11 != 0){

      ##Simulate data from genes that are DE in Exp1 and EE in Exp2 (NOTE: this is
      ##data from Exp1, so genes are DE)

      ##Simulate standard deviations
      s11.1 <- sqrt(rigamma(n=m11, alpha=d0/2, beta=d0*s20/2))

      ##Simulate effect sizes.
      mn112.1 <- rnorm(n=m11, mean=mnd)*sample(c(rep(1,m11/2), rep(-1,m11/2)))*s11.1

      ##Simulate data
      dat111.1 <- t(sapply(s11.1, rnorm, n=ni, mean=0))
      dat112n.1 <- t(sapply(s11.1, rnorm, n=ni, mean=0))
      dat112.1 <- sweep(dat112n.1, 1, mn112.1, "+")
      dat11.1 <- cbind(dat111.1, dat112.1)
    }


    ##Combine all data from Experiment 1 into one matrix
    dat1 <- rbind(dat00.1, dat01.1, dat10.1, dat11.1)


   #########################################
   ##  Simulate data from second experiment#
   #########################################

     ###########################################################################
     ##Simulate data from genes that are EE in both experiments (both treatments
     ##have a mean of 0)

     ##Simulate standard deviations.
      s00.2 <- sqrt(rigamma(n=m00, alpha=d0/2, beta=d0*s20/2))

     ##Simulate data.
      dat00.2 <- t(sapply(s00.2, rnorm, n=2*ni, mean=0))


      ###########################################################
      ##Simulate data from genes that are EE in Exp1 and DE in Exp2 (NOTE: this is
      ##data from Exp2, so genes are DE)

      ##Simulate standard deviations.
        s01.2 <- sqrt(rigamma(n=m01, alpha=d0/2, beta=d0*s20/2))

      ##Simulate effect sizes.
        mn013.2 <- rnorm(n=m01, mean=mnd)*sample(c(rep(1,m01/2), rep(-1,m01/2)))*s01.2

      ##Simulate data.
        dat013n.2 <- t(sapply(s01.2, rnorm, n=ni, mean=0))
        dat013.2 <- sweep(dat013n.2, 1, mn013.2, "+")
        dat014.2 <- t(sapply(s01.2, rnorm, n=ni, mean=0))
        dat01.2 <- cbind(dat013.2, dat014.2)

    ###########################################################
    ##Simulate data from genes that are DE in Exp1 and EE in Exp2 (NOTE: this is
    ##data from Exp2, so all genes are EE and both treatment means are 0)

      ##Simulate standard deviations.
        s10.2 <- sqrt(rigamma(n=m10, alpha=d0/2, beta=d0*s20/2))

      ##Simulate data.
        dat103.2 <- t(sapply(s10.2, rnorm, n=ni, mean=0))
        dat104.2 <- t(sapply(s10.2, rnorm, n=ni, mean=0))
        dat10.2 <- cbind(dat103.2, dat104.2)


    if(m11 == 0) {dat11.2 <- NULL}
    if(m11 != 0){

      ##Simulate data from genes that are DE in Exp1 and EE in Exp2 (NOTE: this is
      ##data from Exp1, so genes are DE)

      ##Simulate standard deviations.
        s11.2 <- sqrt(rigamma(n=m11, alpha=d0/2, beta=d0*s20/2))

      ##Simulate effect sizes.
        mn113.2 <- rnorm(n=m11, mean=mnd)*sample(c(rep(1,m11/2), rep(-1,m11/2)))*s11.2

      ##Simulate data
        dat113n.2 <- t(sapply(s11.2, rnorm, n=ni, mean=0))
        dat113.2 <- sweep(dat113n.2, 1, mn113.2, "+")
        dat114.2 <- t(sapply(s11.2, rnorm, n=ni, mean=0))
        dat11.2 <- cbind(dat113.2, dat114.2)
    }

    ##Combine all data from Experiment 2 into one matrix
    dat2 <- rbind(dat00.2, dat01.2, dat10.2, dat11.2)

    ###########################################################
    ##Perform moderated t-test (Smyth, 2004)
    fit11 <- lmFit(dat1,design1)
    fit12 <- contrasts.fit(fit11,contr.mat1)
    fit13 <- eBayes(fit12)

    ##p-values from Exp1.
    ps1 <- fit13$p.value[,1]


    fit21 <- lmFit(dat2,design2)
    fit22 <- contrasts.fit(fit21,contr.mat2)
    fit23 <- eBayes(fit22)

    ##p-values from Exp2.
    ps2 <- fit23$p.value[,1]


    ##Add p-values to p-values matrix from previous replications.
    pmat <- cbind(pmat, c(ps1, ps2))

#     if(rep%%10 == 0){print(paste("REPLICATION", rep, "IS COMPLETE"))}

  }

  return(pmat) 

}


