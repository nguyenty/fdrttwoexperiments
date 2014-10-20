## Simulation the pvalue from normal distrubition following Phillip paper
# this function simulate n.rep vectors of  pvalue set. Each vector has: first 
# n p-values are from experiment 1, the last n p-values are from experiment 2
# the order of genes is reserved in each half. 
#n: is the number of genes
#pi0: is the proportion of true null
#mu_a: is the mean of mu under the alternative
# n.rep: is the number of replicates of the simulation data
daisy.pvalue <-function(n, pi0, mu_a, n.rep){
  ee <- round(n*pi0)
  de <- n - ee
  sapply(1:n.rep, function(i){
  tvalue.ee <- mvrnorm(n = ee, mu = c(0,0), Sigma = diag(1,2))
  tvalue.de <- mvrnorm(n = de, mu = rep(mu_a, 2), Sigma = diag(1, 2))
  t.value = rbind(tvalue.ee, tvalue.de)
  p.value <- apply(t.value, c(1,2), function(x) 2*(1 - pnorm(abs(x))) )
  return(c(p.value[,1], p.value[,2]))
  path <- paste0("phillips_","n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", n.rep, ".csv" )
  write.csv(p.value, file = path, row.names = F)
}
)
}


daisy.pvalue2 <-function(mv,  mu_a, n.rep){
  n <- sum(mv)
  ee <- sum(mv[1:3])
  de <- n - ee
  sapply(1:n.rep, function(i){
    tvalue.m00 <- mvrnorm(n = mv[1], mu = c(0,0), Sigma = diag(1,2))
    tvalue.m10 <- mvrnorm(n = mv[2], mu = c(mu_a,0), Sigma = diag(1,2))
    tvalue.m01 <- mvrnorm(n = mv[3], mu = c(0,mu_a), Sigma = diag(1,2))
    tvalue.m11 <- mvrnorm(n = mv[4], mu = c(mu_a,mu_a), Sigma = diag(1,2))
    t.value = rbind(tvalue.m00, tvalue.m10, tvalue.m01, tvalue.m11)
    p.value <- apply(t.value, c(1,2), function(x) 2*(1 - pnorm(abs(x))) )
    return(c(p.value[,1], p.value[,2]))
    path <- paste0("phillips_","n_",n, "pi0_", pi0, "mua_", mu_a, "nrep_", n.rep, ".csv" )
    write.csv(p.value, file = path, row.names = F)
  }
  )
}