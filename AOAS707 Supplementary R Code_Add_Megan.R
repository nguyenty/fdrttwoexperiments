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
  csum = cumsum(tess.areas[mrank]) #find cumulative sums sum((areas$summary$dir.area))
  #csum = cumsum(tess.areas[mrank])/sum(tess.areas) #find cumulative sums sum((areas$summary$dir.area))
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
  roc.out <- roc(1-test.vector, lab) # plot(roc.out)
  roc.ind <- sum(roc.out$fpr<=.05)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- auc(roc.out, min =roc.min)
  return(pauc)
}
#?auc

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
  test.vector <- NULL
  test.vector[index] <- cum.areas
  n1 <- round(n*p)
  n11 = round(p1*n) #number of p-vectors (alt, null) m10
  n12 = round(p2*n)
  n0 <- n - n1 - n11 - n12
  lab <- c( rep(1, n1), rep(0, n-n1)) # sum(lab ==1)
  pauc <- auc_out(test.vector, lab)
  #Use mixFDR to calculate estimated left-tail local fdr 
  result = mixFdr(t.values,J=myJ,nearlyNull=nnull,plots=F,P=P,calibrate=cal,theonull=theo.null)  
  k = sum(result$FDRLeft<alpha) + num.zero 
  
  # calculate error rates and power
  myfdr = 0; mypower = all.null = one.null = 0
  if(k>0)
    {
      index.sig = index[1:k]
      myfdr = sum(index.sig>n1)/k
      mypower = sum(index.sig<=n1)/n1
      all.null = sum(index.sig>(n1+n11+n12))/k;
      one.null = sum(index.sig<(n1+n11+n12)&(index.sig>n1))/k
    }
  #record how many p-vectors are significant
  return(list(myfdr=myfdr,mypower=mypower, pauc = pauc, myallnull=all.null,myonenull=one.null))
}


# Preamble for all simulations with positive correllation
#rho = rep(seq(from=0,to=.8,by=.1),each=nrep)
nrep <- 10
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

simout <- function(mualt, nrep, n, p, p1, p2){
  rho <- rep(0, nrep)
  n1 = round(p*n) #number of p-vectors (alt, alt) m11
  n11 = round(p1*n) #number of p-vectors (alt, null) m10
  n12 = round(p2*n) #number of p-vectors (null, alt) m01
  n0 = n - n1 - n11 - n12 # number of null p-vectors m00
  mvs <- c(n0, n12, n11, n1)
  fdr.test.results = ndr.test.results = pauc.test.results = matrix(NA,nrep,6)
  all.null.results = half.null.results = matrix(NA,nrep,6)
  colnames(fdr.test.results) = c("rho","M.FDR","E.FDR","S.FDR","DL.FDR","Ex.FDR")
  colnames(ndr.test.results) = c("rho","M.NDR","E.NDR","S.NDR", "DL.NDR","Ex.NDR")
  colnames(pauc.test.results) = c("rho","M.NDR","E.NDR","S.NDR", "DL.NDR","Ex.NDR")
  colnames(all.null.results) = c("rho","M.allFDR","E.allFDR","S.allFDR","DL.allFDR","Ex.allFDR")
  colnames(half.null.results) = c("rho","M.halfFDR","E.halfFDR","S.halfFDR", "DL.halfFDR","Ex.halfFDR")
  fdr.test.results[,1] = ndr.test.results[,1] = rho
  all.null.results[,1] = half.null.results[,1] = rho
  pauc.test.results[,1] = rho
  
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
    #(cum.areas,index,p=.1,p1=.1,p2=.1,alpha=.05,myJ=2,nnull=1,cal=F,P=NA,theo.null=F)
    sum.E = cumsum(tess.areas[rank.E]); E.results = simtesthalfnull(sum.E,rank.E,p = p, p1= p1, p2 = p2,P=800)
    sum.M = cumsum(tess.areas[rank.M]); M.results = simtesthalfnull(sum.M,rank.M,p = p, p1= p1, p2 = p2,P=800)
    sum.S = cumsum(tess.areas[rank.S]); S.results = simtesthalfnull(sum.S,rank.S,p = p, p1= p1, p2 = p2,P=800)
    sum.DL = cumsum(tess.areas[rank.DL]); DL.results =simtesthalfnull(sum.DL,rank.DL,p = p, p1= p1, p2 = p2,P=800)
  
    # obtain pAUC of those scores 
  #record results
    ndr.test.results[i,3] = E.results$mypower; fdr.test.results[i,3] = E.results$myfdr
    ndr.test.results[i,2]= M.results$mypower; fdr.test.results[i,2] = M.results$myfdr
    ndr.test.results[i,4]  = S.results$mypower; fdr.test.results[i,4]  = S.results$myfdr
    ndr.test.results[i,5] = DL.results$mypower; fdr.test.results[i,5] = DL.results$myfdr
    
  pauc.test.results[i,3] = E.results$pauc
  pauc.test.results[i,2]= M.results$pauc
  pauc.test.results[i,4]  = S.results$pauc
  pauc.test.results[i,5] = DL.results$pauc
  
    all.null.results[i,3] = E.results$myallnull; half.null.results[i,3] = E.results$myonenull
    all.null.results[i,2]= M.results$myallnull; half.null.results[i,2] = M.results$myonenull
    all.null.results[i,4]  = S.results$myallnull; half.null.results[i,4]  = S.results$myonenull
    all.null.results[i,5] = DL.results$myallnull; half.null.results[i,5] = DL.results$myonenull
    
    # get comparison results using existing method
    max.pvalues = apply(my.data,1,max)
    pauc.test.results[i,6] <- auc_out(max.pvalues, lab) 
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
fdr.means = ndr.means = all.null.fdr.means = half.null.fdr.means = pauc.means = rep(0, 6)
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



ndr.means <- apply(ndr.test.results*n1,2,mean)
ndr.se <- apply(ndr.test.results*n1,2,sd)/sqrt(nrep)

pauc.means <- apply(pauc.test.results,2,mean)
pauc.se <- apply(pauc.test.results,2,sd)/sqrt(nrep)

fdr.means <-  apply(fdr.test.results,2,mean)
fdr.se <-  apply(fdr.test.results,2,sd)/sqrt(nrep)

out1 <- cbind(ndr.means, ndr.se, fdr.means, fdr.se, pauc.means, pauc.se)
rownames(out1) <- c("rho", "Maximum", "Euclidean", "Summation", "De Lichtenberg", "Existing")



# all.null.fdr.means <- apply(all.null.results,2,mean)
# half.null.fdr.means <- apply(half.null.results,2,mean)


#print the tables
megan_out <- megan_ints_out(pmat, mvs)




out <- rbind(out1, Histogram = megan_out)
out
}


pm1 <- proc.time()
input <- c(mualt= 3, nrep = 30,  n = 10000, p = .1, p1 = .1, p2 = .1)

res3 <- simout(input[1], input[2], input[3], input[4], input[5], input[6])

proc.time() -pm1
file <- paste0("mualt_",input[1], "nrep_", input[2], 
               "n_", input[3], "p_", input[4], 
               "p1_", input[5], "p2_", input[6], ".csv" )
write.csv(res3, file, row.names = T)
