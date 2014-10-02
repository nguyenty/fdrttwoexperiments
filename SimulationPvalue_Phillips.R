## Simulation the pvalue from normal distrubition following Phillip papaer
sim.pvalue <-function(n, pi0, mu_a, n.rep){
  ee <- round(n*pi0)
  de <- n - ee
  sapply(1:n.rep, function(i){
  tvalue.ee <- mvrnorm(n = ee, mu = c(0,0), Sigma = diag(1,2))
  tvalue.de <- mvrnorm(n = de, mu = rep(mu_a, 2), Sigma = diag(1, 2))
  t.value = rbind(tvalue.de, tvalue.ee)
  p.value <- apply(t.value, c(1,2), function(x) 2*(1 - pnorm(abs(x))) )
  if (i %%10 ==0) print(paste("REPLICATION", i, "IS COMPLETE"))
  return(c(p.value[,1], p.value[,2]))
}
)
}
