library(limma)
library(pscl)

##The function "micro_sim_pvalues" simulates data from existing
##microarray data sets and calculates limma p-values.  The output
##of this function is a matrix of p-values with 2*m rows and nrep 
##columns. The columns correspond to p-values from a pair of
##experiments.  The first m values correspond to the p-values from
##the first experiment and the second m p-values correspond to
##p-values from the second experiment.  The inputs for this function
##are as follows:

  ##y1: first microarray data set.  For the previous simulations,
    ##this data set was a portion of the data described in 
    ##Metzeler et al. (2008).  This data is in the file
    ##"micro_GSE12417_samp1.txt".
  ##y2: second microarray data set.  For the previous simulations,
    ##this data set was a portion of the data described in 
    ##Metzeler et al. (2008).  This data is in the file
    ##"micro_GSE12417_samp2.txt".
  ##ni: the sample size in each treatment group
  ##mv: a vector with the following values: (m00, m01, m10, m11),
    ##where m00 is the number of genes that are EE in both
    ##experiments, m01 is the number of genes that are EE in Exp1
    ##but DE in Exp2, m10 is the number of genes that are DE in
    ##Exp1 but EE in Exp2, and m11 is the number of genes that are
    ##DE in both experiments.  NOTE: m = sum(mv).
  ##mnd: the (absolute) mean of the distribution of relative effect
    ##sizes
  ##nrep: the number of pairs of data sets to simulate.

micro_sim_pvalues <- function(y1, y2, ni, mv, mnd, nrep){
  pmat <- NULL

  m00 <- mv[1]
  m01 <- mv[2]
  m10 <- mv[3]
  m11 <- mv[4]

  m <- sum(mv)
  print(paste("m =",m))

  ##number of units in each microarray data set
  ns1 <- dim(y1)[2]
  ns2 <- dim(y2)[2]

  ##standard deviations for each gene in each data set
  sd1 <- apply(y1, 1, sd)
  sd2 <- apply(y2, 1, sd)

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

    ##Randomly select units to be in treatment 1 (idsamp1) and
    ##treatment 2 (idsamp2).  I'm not sure why I did this,
    ##but I sample from y2 (the second data set) first and
    ##then y1 after.
      idsamp <- sample(1:ns2, size=2*ni, replace=FALSE)
      idsamp1 <- idsamp[1:ni]
      idsamp2 <- idsamp[(ni+1):(2*ni)]

    ##Randomly select the m00 genes that are EE in both
    ##experiments, the m01 genes that are EE in Exp1 but
    ##DE in Exp2, the m10 genes that are DE in Exp1 but
    ##EE in Exp2, and the m11 genes that are DE in both
    ##experiments.
      g00 <- sample(1:m, size=m00, replace=FALSE)
      gns1 <- setdiff(1:m, g00)
      g01 <- sample(gns1, size=m01, replace=FALSE)
      gns2 <- setdiff(1:m, c(g00, g01))
      g10 <- sample(gns2, size=m10, replace=FALSE)
      g11 <- setdiff(1:m, c(g00, g01, g10))


    ###########################################################
    ##Data from genes that are EE in both experiments (both treatments
    ##have a mean of 0)
      dat00.1 <- y2[g00, idsamp]


    ###########################################################
    ##Data from genes that are EE in Exp1 and DE in Exp2 (This is data
    ##from Exp1, so genes are EE and have the same population mean).
      dat011 <- y2[g01, idsamp1]
      dat012 <- y2[g01, idsamp2]
      dat01.1 <- cbind(dat011, dat012)


    ###########################################################
    ##Simulate data from genes that are DE in Exp1 and EE in Exp2 (NOTE: this is
    ##data from Exp1, so genes are DE)
      ##Simulate relative effect sizes from genes with negative effect sizes.
        mns10neg <- rnorm(n=round(m10*0.5), mean=-mnd) 

      ##Simulate relative effect sizes from genes with positive effect sizes.
        mns10pos <- rnorm(n=round(m10*0.5), mean=mnd)

      ##Randomize order of effect sizes.
        mns10 <- sample(c(mns10neg, mns10pos))*sd2[g10]

      ##Simulate data.  
        dat101 <- y2[g10, idsamp1]
        dat102 <- sweep(y2[g10, idsamp2], 1, mns10, FUN="+")
        dat10.1 <- cbind(dat101, dat102)


    ###########################################################
    ##Simulate data from genes that are DE in Exp1 and EE in Exp2 (NOTE: this is
    ##data from Exp1, so genes are DE)

      ##Simulate relative effect sizes from genes with negative effect sizes.
        mns112neg <- rnorm(n=round(m11*0.5), mean=-mnd) 

      ##Simulate relative effect sizes from genes with positive effect sizes.
        mns112pos <- rnorm(n=round(m11*0.5),, mean=mnd)

      ##Randomize order of effect sizes.
        mns112 <- sample(c(mns112neg, mns112pos))*sd2[g11]

      ##Simulate data.
        dat111 <- y2[g11, idsamp1]
        dat112 <- sweep(y2[g11, idsamp2], 1, mns112, FUN="+")
        dat11.1 <- cbind(dat111, dat112)


    ##Combine all data from Experiment 1 into one matrix.
        dat1 <- rbind(dat00.1, dat01.1, dat10.1, dat11.1)



   #########################################
   ## Simulate data from second experiment## 
   #########################################

    ##Randomly select units to be in treatment 1 (id2exp1) and
    ##treatment 2 (id2exp2).  
      id2exp <- sample(1:ns1, size=2*ni, replace=FALSE)
      id2exp1 <- id2exp[1:ni]
      id2exp2 <- id2exp[(ni+1):(2*ni)]


    ###########################################################
    ##Data from genes that are EE in both experiments (both treatments
    ##have a mean of 0)
      dat00.2 <- y1[g00, id2exp]

    
    ###########################################################
    ##Simulate data from genes that are EE in Exp1 and DE in Exp2 (NOTE: this is
    ##data from Exp2, so genes are DE)

      ##Simulate relative effect sizes from genes with negative effect sizes.
        mns01neg2exp <- rnorm(n=round(m01*0.5), mean=-mnd) 

      ##Simulate relative effect sizes from genes with positive effect sizes.  
        mns01pos2exp <- rnorm(n=round(m01*0.5), mean=mnd)

      ##Randomize order of effect sizes.
        mns012exp <- sample(c(mns01neg2exp, mns01pos2exp))*sd1[g01]


      ##Simulate data.
        dat012exp1 <- y1[g01, id2exp1]
        dat012exp2 <- sweep(y1[g01, id2exp2], 1, mns012exp, FUN="+")
        dat01.2 <- cbind(dat012exp1, dat012exp2)


    ###########################################################
    ##Data from genes that are DE in Exp1 and EE in Exp2 (This is data
    ##from Exp2, so genes are EE and have the same population mean).

      ##Simulate data.
        dat102exp1 <- y1[g10, id2exp1]
        dat102exp2 <- y1[g10, id2exp2]
        dat10.2 <- cbind(dat102exp1, dat102exp2)


    ###########################################################
    ##Simulate data from genes that are DE in Exp1 and EE in Exp2 (NOTE: this is
    ##data from Exp1, so genes are DE)

      ##Simulate relative effect sizes from genes with negative effect sizes.
        mns11neg2exp <- rnorm(n=round(m11*0.5), mean=-mnd) 

      ##Simulate relative effect sizes from genes with positive effect sizes.
        mns11pos2exp <- rnorm(n=round(m11*0.5), mean=mnd)

      ##Randomize order of effect sizes.
        mns112exp <- sample(c(mns11neg2exp, mns11pos2exp))*sd1[g11]

      ##Simulate data.
        dat112exp1 <- y1[g11, id2exp1]
        dat112exp2 <- sweep(y1[g11, id2exp2], 1, mns112exp, FUN="+")
        dat11.2 <- cbind(dat112exp1, dat112exp2)


    ##Combine all data from Experiment 2 into one matrix.
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


    ##Add p-values to matrix of p-values from previous replications.
    pmat <- cbind(pmat, c(ps1, ps2))    

    if(rep%%10 == 0) {print(paste("REPLICATION", rep, "IS COMPLETE"))}
  }

  return(pmat)
}


