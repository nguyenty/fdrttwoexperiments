source("estimate_m0s_fn.R")
library(qvalue)


##find largest p-value corresponding to largest q-value <= qmx
##(used in qvalue_intersection function below)
pvco_fn <- function(qmx, qv, pv){
  difq <- abs(qmx - qv) + (qv > qmx)
  pvco <- max(pv[difq == min(difq)])
  return(pvco)
}


##calculates max(q1, q2) for each gene
qvalue_maximum <- function(pv1, pv2){
  
  co1 <- calc.cutoff(pv1)
  qv1 <- qvalue(pv1, lambda=co1)$qvalue

  co2 <- calc.cutoff(pv2)
  qv2 <- qvalue(pv2, lambda=co2)$qvalue

  qvmax <- apply(cbind(qv1,qv2), 1, max)

  return(qvmax)
}


##calculates q-value for each gene using proposed method
qvalue_intersection <- function(pv1, pv2){

    ##calculate lambda cutoff and qvalues for experiment 1
    co1 <- calc.cutoff(pv1)
    qv1 <- qvalue(pv1, lambda=co1)$qvalue

    ##calculate lambda cutoff and qvalues for experiment 2
    co2 <- calc.cutoff(pv2)
    qv2 <- qvalue(pv2, lambda=co2)$qvalue

    m <- length(pv1)
    pv12 <- cbind(pv1, pv2)
    p0 <- estimate.m0s(p1=pv1, p2=pv2)$ms/m
    names(p0) <- NULL
    p0.1 <- p0[2]  ##pi0 estimate for experiment 1
    p0.2 <- p0[3]  ##pi0 estimate for experiment 2
    p00 <- p0[4]   ##pi00 estimate

    p01 <- max(p0.1-p00, 0)  ##pi01 estimate
    p10 <- max(p0.2-p00, 0)  ##pi10 estimate
    
    if(p0.1 >= 1){ p10 <- 0 }
    if(p0.2 >= 1){ p01 <- 0 }

    #qmx <- qvalue_maximum(pv1=pv1, pv2=pv2)
    qmx <- apply(cbind(qv1,qv2), 1, max)

    ##determine rectangular p-value rejection region for each gene
    pco1 <- sapply(qmx, pvco_fn, qv=qv1, pv=pv1)
    pco2 <- sapply(qmx, pvco_fn, qv=qv2, pv=pv2)

    pco12 <- cbind(pco1, pco2)     


    ##calculate proposed q-value for each gene
    Pp1 <- sapply(pco1, function(x) mean(pv1 <= x))
    Pp2 <- sapply(pco2, function(x) mean(pv2 <= x))
    Pp1[Pp1 < pco1] <- pco1[Pp1 < pco1]
    Pp2[Pp2 < pco2] <- pco2[Pp2 < pco2]

    Pp1p2 <- apply(pco12, 1, function(x) mean(pv1 <= x[1] & pv2 <= x[2]))
    Pp1p2[Pp1p2 < (pco1*pco2)] <- (pco1*pco2)[Pp1p2 < (pco1*pco2)]

    Pp1H1c <- rep(1, m)
    Pp2H2c <- rep(1, m)

    if(p0.1 < 1){
       Pp1H1c <- (Pp1 - pco1*p0.1)/(1 - p0.1)
       Pp1H1c[Pp1H1c < pco1] <- pco1[Pp1H1c < pco1]
    }

    if(p0.2 < 1){
       Pp2H2c <- (Pp2 - pco2*p0.2)/(1 - p0.2)
       Pp2H2c[Pp2H2c < pco2] <- pco2[Pp2H2c < pco2]
    }

    qint1 <- (pco1*pco2*p00 + pco1*Pp2H2c*p01 + pco2*Pp1H1c*p10)/Pp1p2

    qint <- sapply(qmx, function(x) min(qint1[qmx >= x]))

    return(qint)


}

