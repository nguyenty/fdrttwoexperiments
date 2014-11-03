# ---------------------------------------------------------------------------------------------#
# Description of R code:
# 1. Functions needed for Voronoi P-value Combination for the disjunction hypothesis.
# 2. Code used to perform simulations as described in manuscript.
# 3. Code used to perform 'challenging situation' simulations as described in supplementary materials.
# --------------------------------------------------------------------------------------------- #

#First, load libraries
#install.packages("deldir")
library(deldir)
library(mixfdr)
library(MASS)
library(xtable)
library(AUC)
source("MeganDaisyRes.R")
source("estimate_m0s_fn.R")
source("qvalue_functions.R")
source("SimulationPvalue_Daisy.R")
source("SimulationPvalue_Phillips.R")
source("MeganDaisyRes.R")
source("normal_sim_pvalues.R")
library(qvalue)
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
  
  my.data = cbind(x.values,y.values) #x.values <- ps[1:n, 1]; y.values <- ps[n+1:n, 1]; ranking = 1
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
  csum = cumsum(tess.areas[mrank])/sum(tess.areas) #find cumulative sums sum((areas$summary$dir.area))
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

## AUC function ####

auc_out <- function(test.vector, lab){
  lab <- as.factor(lab)
  roc.out <- roc(1-test.vector, lab)
  roc.ind <- sum(roc.out$fpr<=.1)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- auc(roc.out, min =roc.min)
  return(pauc)
}


#-------------------------------------------------------------------------------------#
#     PART 3. Simulations for results presented in supplementary materials
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
# simtesthalfnull: Function that performs analysis on cumulative voronoi areas using mixFdr, and 
#          outputs quantities of interest to study properties of procedure.
# INPUTS: cum.areas (vector) = cumulative cell areas, output from CombineVoronoi
#         index (vector) = index of these areas in terms of original p-value vectors.  
#                          output from CombiniVoronoi
#         p (value) = known proportion of 'true alternative hypotheses'.
#         p1 (value) = known proportion of p-vectors with the first signal from the 
#                      alternative distribution, and the second from the null distribution.
#         p2 (value) = known proportion of p-vectors with the first signal from the 
#                      null distribution, and the second from the alternative.
#         alpha (value) = nominal level of FDR control
#         myJ (integer) = parameter to pass to mixFdr, number of assumed distributions.
#         nnull (value) = parameter to pass to mixFdr, used to distinguish null from 
#                         alternative empirical distributions
#         cal (T/F, or value) = parameter to pass to mixFdr, whether to calibrate
#         P (value) = paramater to pass to misFdr, penalization factor that influences 
#                     tendency to form one large emprical null.
#         theo.null (T/F) = parameter to pass to mixFdr, whether to fit theoretical null.
# OUTPUTS: myfdr (value) = known proportion of false discoveries ( considered false according to disjunction)
#          mypower (value) = known proportion of correctly discoveried alternative hypotheses (1-ndr) 
#          myallnull (value) = known proportion of false discoveries that came from two null signals
#          myonenull (value) = known proportion of false discoveries that came from one null, one alternative signal
#-------------------------------------------------------------------------------------#
simtesthalfnull = function(cum.areas,index,p=.1,p1=.1,p2=.1,alpha=.05,myJ=2,nnull=1,cal=F,P=NA,theo.null=F)
  {
  n = length(index)
  cum.areas.nz = cum.areas[cum.areas>0&cum.areas<0.9999683]  #Choose the areas with cumulative sum above 0 and below a cut-off (keeps transformation from going crazy)
  num.zero = sum(cum.areas==0)
  t.values = qnorm(cum.areas.nz)                         #Transform cumulative areas using standard normal quantile function
  
  #Use mixFDR to calculate estimated left-tail local fdr 
  result = mixFdr(t.values,J=myJ,nearlyNull=nnull,plots=F,P=P,calibrate=cal,theonull=theo.null)  
  k = sum(result$FDRLeft<alpha) + num.zero 
  
  # calculate error rates and power
  myfdr = 0; mypower = all.null = one.null = 0
  if(k>0)
    {
      index.sig = index[1:k]
      myfdr = sum(index.sig>(n*p))/k
      mypower = sum(index.sig<=(n*p))/(n*p)  
      all.null = sum(index.sig>(n*p+n*p1+n*p2))/k;
      one.null = sum(index.sig<(n*p+n*p1+n*p2)&(index.sig>n*p))/k
    }
  #record how many p-vectors are significant
  return(list(myfdr=myfdr,mypower=mypower,myallnull=all.null,myonenull=one.null))
}


# Preamble for all simulations with positive correllation
#rho = rep(seq(from=0,to=.8,by=.1),each=nrep)
nrep <- 100
rho <- rep(0, nrep)
n = 2000 #number of genes for each data set
p = .1 # number of true alternative signals (alt,alt)
p1 = .1 #proportion of p-vectors (alt, null)
p2 = .1 #proportion of p-vectors (null, alt)
n1 = p*n #number of p-vectors (alt, alt) m11
n11 = p1*n #number of p-vectors (alt, null) m10
n12 = p2*n #number of p-vectors (null, alt) m01
n0 = n*(1-p-p1-p2) # number of null p-vectors m00
mualt = 3 #mean of alternative signal

# --------------------------# SUPPLEMENTARY 1 #------------------------------------- #
# Simulation where mualt=3, and 10% are halfnull  (3,0), 10% are halfnull (0,3) #

#initialize matrices to store results

PMod=800
pmat <- NULL
mvs <- c(n0, n12, n11, n1)
simout <- function(mualt, n){
  fdr.test.results = ndr.test.results = matrix(NA,nrep,6)
  all.null.results = half.null.results = matrix(NA,nrep,6)
  colnames(fdr.test.results) = c("rho","M.FDR","E.FDR","S.FDR","DL.FDR","Ex.FDR")
  colnames(ndr.test.results) = c("rho","M.NDR","E.NDR","S.NDR", "DL.NDR","Ex.NDR")
  colnames(all.null.results) = c("rho","M.allFDR","E.allFDR","S.allFDR","DL.allFDR","Ex.allFDR")
  colnames(half.null.results) = c("rho","M.halfFDR","E.halfFDR","S.halfFDR", "DL.halfFDR","Ex.halfFDR")
  fdr.test.results[,1] = ndr.test.results[,1] = rho
  all.null.results[,1] = half.null.results[,1] = rho
  pauc.E <- pauc.M <- pauc.S <- pauc.DL <- pauc.Ex <- NULL
for(i in 1:nrep)
  {
    #get data
    set.seed(i) # i <- 1
    if(i==64)
    {
      set.seed(128) # something strange happens when i=64
    }
    r = rho[i]
    my.zvals = rbind(mvrnorm(n1,c(mualt,mualt),matrix(c(1,r,r,1),2,2)), # m11
                      mvrnorm(n11,c(mualt,0),matrix(c(1,r,r,1),2,2)), #m10
                      mvrnorm(n12,c(0,mualt),matrix(c(1,r,r,1),2,2)), # m01
                      mvrnorm(n0,c(0,0),matrix(c(1,r,r,1),2,2))) # m00
    my.data = 2*pnorm(-abs(my.zvals))
    x.values = my.data[,1]; y.values = my.data[,2]
    megan.data <- my.data[n+1 - 1:n, ] # m00, m01, m10, m11
    ps1 <- megan.data[, 1]
    ps2 <- megan.data[, 2]
    pmat <- cbind(pmat, c(ps1, ps2))
    #get voronoi tessellation and extract cell areas
    areas = deldir(x.values,y.values,rw=c(0,1,0,1),digits=20,eps=1e-13) 
    tess.areas = areas$summary$dir.area 
    
    #get rankings and indices
    distance.E = apply(my.data,1,function(x){sqrt(x[1]^2+x[2]^2)}); rank.E = sort(distance.E,index.return=T)$ix
    distance.M = apply(my.data,1,max) ; rank.M = sort(distance.M,index.return=T)$ix
    distance.S = apply(my.data,1,sum); rank.S = sort(distance.S,index.return=T)$ix
    distance.DL = apply(my.data,1,function(x){prod(x)*(1+(x[1]/.001)^2)*(1+(x[2]/.001)^2)}); rank.DL = sort(distance.DL,index.return=T)$ix
    
    #get cumulative sums, then myfdr and mypower for each ranking scheme
    sum.E = cumsum(tess.areas[rank.E]); E.results = simtesthalfnull(sum.E,rank.E,P=800)
    sum.M = cumsum(tess.areas[rank.M]); M.results = simtesthalfnull(sum.M,rank.M,P=800)
    sum.S = cumsum(tess.areas[rank.S]); S.results = simtesthalfnull(sum.S,rank.S,P=800)
    sum.DL = cumsum(tess.areas[rank.DL]); DL.results =simtesthalfnull(sum.DL,rank.DL,P=800)
  
    # obtain pAUC of those scores 
    lab <- c( rep(1, n1), rep(0, n-n1))
    pauc.E[i] <- auc_out(sum.E, lab)
    pauc.M[i] <- auc_out(sum.M, lab)
    pauc.S[i] <- auc_out(sum.S, lab)
    pauc.DL[i] <- auc_out(sum.DL, lab)
    
    #record results
    ndr.test.results[i,3] = E.results$mypower; fdr.test.results[i,3] = E.results$myfdr
    ndr.test.results[i,2]= M.results$mypower; fdr.test.results[i,2] = M.results$myfdr
    ndr.test.results[i,4]  = S.results$mypower; fdr.test.results[i,4]  = S.results$myfdr
    ndr.test.results[i,5] = DL.results$mypower; fdr.test.results[i,5] = DL.results$myfdr
    
    all.null.results[i,3] = E.results$myallnull; half.null.results[i,3] = E.results$myonenull
    all.null.results[i,2]= M.results$myallnull; half.null.results[i,2] = M.results$myonenull
    all.null.results[i,4]  = S.results$myallnull; half.null.results[i,4]  = S.results$myonenull
    all.null.results[i,5] = DL.results$myallnull; half.null.results[i,5] = DL.results$myonenull
    
    # get comparison results using existing method
    max.pvalues = apply(my.data,1,max)
    pauc.Ex[i] <- auc_out(max.pvalues, lab) 
    max.bh = BH(max.pvalues)$index
    k = length(max.bh)
    if(k==0)
    {
      all.null.results[i,6] = fdr.test.results[i,6] = half.null.results[i,6] = ndr.test.results[i,6] = 0
    }
    if(k>0)
    {
      fdr.test.results[i,6] = sum(max.bh>n1)/k
      ndr.test.results[i,6] = sum(max.bh<=n1)/n1
      all.null.results[i,6] = sum(max.bh>(n1+n*p1+n*p2))/k
      half.null.results[i,6] = sum((max.bh>n1)&(max.bh<(n1+n*p1+n*p2)))/k
    }
}

# get tables of means
fdr.means = ndr.means = all.null.fdr.means = half.null.fdr.means = rep(0, 6)
names(fdr.means) = names(ndr.means) = c("rho","Maximum","Euclidean","Summation","De Lichtenberg","Existing")
names(all.null.fdr.means) = names(half.null.fdr.means) = c("rho","Maximum","Euclidean","Summation","De Lichtenberg","Existing")
# 
# for(i in 0){
#   temp.low = nrep*i+1;   temp.high = nrep*(i+1)
#   ndr.means[i+1,] = apply(ndr.test.results[temp.low:temp.high,],2,mean)
#   fdr.means[i+1,] = apply(fdr.test.results[temp.low:temp.high,],2,mean)
#   all.null.fdr.means[i+1,] = apply(all.null.results[temp.low:temp.high,],2,mean)
#   half.null.fdr.means[i+1,] = apply(half.null.results[temp.low:temp.high,],2,mean)
# }



ndr.means <- apply(ndr.test.results*n11,2,mean)
ndr.se <- apply(ndr.test.results*n11,2,sd)/10
fdr.means <-  apply(fdr.test.results,2,mean)
fdr.se <-  apply(fdr.test.results,2,sd)/10
out1 <- cbind(ndr.means, ndr.se, fdr.means, fdr.se)
rownames(out1) <- c("rho", "Maximum", "Euclidean", "Summation", "De Lichtenberg", "Existing")

pauc.E.means <- mean(pauc.E)
pauc.E.se <- sd(pauc.E)/10
pauc.M.means <- mean(pauc.M)
pauc.M.se <- sd(pauc.M)/10
pauc.S.means <- mean(pauc.S)
pauc.S.se <- sd(pauc.S)/10
pauc.DL.means <- mean(pauc.DL)
pauc.DL.se <- sd(pauc.DL)/10
pauc.Ex.means <- mean(pauc.Ex)
pauc.Ex.se <- sd(pauc.Ex)/10

out2 <- cbind(pauc.mean = c(0,pauc.M.means, pauc.E.means, pauc.S.means, pauc.DL.means,pauc.Ex.means ), 
          pauc.se = c(0,pauc.M.se, pauc.E.se, pauc.S.se, pauc.DL.se,pauc.Ex.se ))
out3 <- cbind(out1, out2)
# all.null.fdr.means <- apply(all.null.results,2,mean)
# half.null.fdr.means <- apply(half.null.results,2,mean)


#print the tables
megan_out <- megan_ints_out(pmat, mvs)




out <- rbind(out3, Histogram = megan_out)
out
}

res2 <- simout(2)
res2.5 <- simout(2.5)
res2.8 <- simout(2.8)
res3 <- simout(3, 2000)
res3.3 <- simout(3.3)

res2[, c(2,4)] <- res2[, c(2,4)] /10
res2.5[, c(2,4)] <- res2.5[, c(2,4)] /10
res2.8[, c(2,4)] <- res2.8[, c(2,4)] /10
res3[, c(2,4)] <- res3[, c(2,4)] /10
res3.3[, c(2,4)] <- res3.3[, c(2,4)] /10
res2
res2.5
res2.8
res3
res3.3

