# # auc
# library("AUC")
# pauc <- function(test.vector, lab){
#   roc.out <- roc(1-test.vector, lab)
#   roc.ind <- sum(roc.out$fpr<=.1)
#   roc.min <- roc.out$cutoffs[roc.ind]
#   pauc_out <- auc(roc.out, min =roc.min)
#   return(pauc_out)
# }

ints_out <- function(ps, pi0){
  n <- dim(ps)[1]/2; n.rep <- dim(ps)[2]; ee <- round(n*pi0); de <- n - ee
  qints <- apply(ps, 2, function(x) qvalue_intersection(pv1=x[1:n], pv2=x[(n+1):(2*n)]))
  res <- apply(qints, 2, function(x){# x<- qints[,7]
    V.t <- sum(x[de+1:(n-de)] <= 0.05)
    R.t <- sum(x <= 0.05)
    S.t <- R.t - V.t
    FDR.t <- V.t/max(R.t, 1)
#     PAUC.t <- pauc(x, lab=as.factor(c(rep(1, de), rep(0, ee))))
    return(c(S.t = S.t, FDR.t = FDR.t))
  })
  meanS <- mean(res[1,]); sdmeanS <- sd(res[1,])/sqrt(n)  
  meanFDR <- mean(res[2,]); sdmeanFDR <-sd(res[2,])/sqrt(n)
#   meanPAUC <- mean(res[3,]); sdmeanPAUC <-sd(res[3,])/sqrt(n)
  
  res2 <- c(meanS = meanS, sdmeanS = sdmeanS, meanFDR = meanFDR, 
            sdmeanFDR = sdmeanFDR)
  res2
}



voronoi_out <- function(ps, pi0){ # dim(ps)
  n <- dim(ps)[1]/2; n.rep <- dim(ps)[2]; ee <- round(n*pi0); de <- n - ee
  voronoi <- apply(ps, 2, function(x) CombineVoronoi(x.values=x[1:n], y.values=x[(n+1):(2*n)]))
  res <- t(laply(1:n.rep, function(i) { # i <- 1
    d <- as.data.frame(cbind(voronoi[[i]][[1]], voronoi[[i]][[2]]))
    dd <- (arrange(d, V1)) # rearrange the gene as the original order 1, 2, 3, ,,,
    bh <- BH(dd[,2]) # which.min(dd[,2])
    R.t <- bh$k
    S.t <- ifelse(R.t >0,sum(bh$index <=de), 0)
    V.t <- R.t - S.t
    FDR.t <- V.t/max(R.t, 1)
#     PAUC.t <- pauc(dd[,2], lab=as.factor(c(rep(1, de), rep(0, ee))))
    return(c(S.t = S.t, FDR.t = FDR.t))
  }))
  meanS <- mean(res[1,]); sdmeanS <- sd(res[1,])/sqrt(n)  
  meanFDR <- mean(res[2,]); sdmeanFDR <-sd(res[2,])/sqrt(n)
#   meanPAUC <- mean(res[3,]); sdmeanPAUC <-sd(res[3,])/sqrt(n)
  
  res2 <- c(meanS = meanS, sdmeanS = sdmeanS, meanFDR = meanFDR, 
            sdmeanFDR = sdmeanFDR)
  res2
}