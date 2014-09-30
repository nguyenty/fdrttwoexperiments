

###################################################################################
###################################################################################
##calc.cutoff slightly modifies the estimate.m0 function to return the appropriate
##p-value cutoff where all p-values greater than the cutoff are assumed to come
##from null cases.

calc.cutoff = function(p, B = 20, max=1){

  m <- length(p)
  m0 <- m
  bin <- c(-0.1, (1:B)/B*max)
  bin.counts=rep(0,B)

  for(i in 1:B){
    bin.counts[i]=sum((p>bin[i])&(p<=bin[i+1]))
  }
  tail.means <- rev(cumsum(rev(bin.counts))/(1:B))
  temp <- bin.counts - tail.means
  index <- min((1:B)[temp <= 0])
  cutoff2 <- (index)/B*max
  if(cutoff2 == 1) {cutoff2 <- 1-1/B}

  return(cutoff2)

}


###################################################################################
###################################################################################

estimate.m0s <- function(p1, p2, B=20){
  m <- length(p1)

  ##find lambda cutoffs using histogram-based method
  c1 <- calc.cutoff(p1, B=B, max=1)
  c2 <- calc.cutoff(p2, B=B, max=1)

  ##estimate m0 for experiment 1
  ind1 <- (p1>=c1)
  m0.1 <- sum(ind1)/(1-c1)

  ##estimate m0 for experiment 2  
  ind2 <- (p2>=c2)
  m0.2 <- sum(ind2)/(1-c2)


  ##estimate m00
  ind12 <- ind1 & ind2
  nA <- sum(ind12)
  pA <- (1-c1)*(1-c2)
  m00 <- nA/pA


  ret <- list()
  ret$ms <- c(m, m0.1, m0.2, m00)
  names(ret$ms) <- c("m", "m0.1", "m0.2", "m00")
  ret$cutoffs <- c(c1, c2)
  return(ret)
}





