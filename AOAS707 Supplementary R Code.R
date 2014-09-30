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
  
  my.data = cbind(x.values,y.values)
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
  csum = cumsum(tess.areas[mrank]) #find cumulative sums
  
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

#-------------------------------------------------------------------------------------#
#     PART 2. Simulations for results presented in manuscript
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
# simtest: Function that performs analysis on cumulative voronoi areas using mixFdr, and 
#          outputs quantities of interest to study properties of procedure.
# INPUTS: cum.areas (vector) = cumulative cell areas, output from CombineVoronoi
#         index (vector) = index of these areas in terms of original p-value vectors.  
#                          output from CombiniVoronoi
#         p (value) = known proportion of 'true alternative hypotheses'.
#         alpha (value) = nominal level of FDR control
#         myJ (integer) = parameter to pass to mixFdr, number of assumed distributions.
#         nnull (value) = parameter to pass to mixFdr, used to distinguish null from 
#                         alternative empirical distributions
#         cal (T/F, or value) = parameter to pass to mixFdr, whether to calibrate
#         P (value) = paramater to pass to misFdr, penalization factor that influences 
#                     tendency to form one large emprical null.
#         theo.null (T/F) = parameter to pass to mixFdr, whether to fit theoretical null.
# OUTPUTS: myfdr (value) = known proportion of false discoveries
#          mypower (value) = known proportion of correctly discoveried alternative hypotheses (1-ndr) 
#-------------------------------------------------------------------------------------#
simtest = function(cum.areas,index,p=.1,alpha=.05,myJ=2,nnull=1,cal=F,P=NA,theo.null=F)
{
  n = length(index)
  cum.areas.nz = cum.areas[cum.areas>0&cum.areas<0.9999683]  #Choose the areas with cumulative sum above 0 and below a cut-off (keeps transformation from going crazy)
  num.zero = sum(cum.areas==0)
  t.values = qnorm(cum.areas.nz)    #Transform cumulative areas using standard normal quantile function
  
  result = mixFdr(t.values,J=myJ,nearlyNull=nnull,plots=F,P=P,calibrate=cal,theonull=theo.null)#Use mixFDR to find left tail lfdr estimates
  k = sum(result$FDRLeft<alpha) + num.zero # determine number of rejections
  
  # calculate and report known values of 'fdr' and '1-ndr'
  myfdr = 0; mypower = 0
  if(k>0){index.sig = index[1:k]
          myfdr = sum(index.sig>(n*p))/k
          mypower = sum(index.sig<=(n*p))/(n*p) }                                           
  return(list(myfdr=myfdr,mypower=mypower))
}

# Preamble for all simulations
rho = rep(seq(from=0,to=.8,by=.1),each=100) #vector of values for correlation
n = 2000 # number of genes
p = .1 # proportion of true alternative signals
n1 = p*n # number of true alternative signals for each data set
n0 = n*(1-p) # number of true null signals for each data set

# -----------------------# Weak alternative simulation #--------------------------------- #

# Note - these simulations are written without using CombineVoronoi, they include all steps 
# without the function #

mualt <- 2     # set mean of alternative signals
Pweak <- 300  # set penalization factor to pass to mixFdr

# initialize matrices to store results
fdr.weak.results = ndr.weak.results = matrix(NA,length(rho),7)
colnames(fdr.weak.results) = c("alt","rho","E.FDR","M.FDR","S.FDR","DL.FDR","EX.FDR")
colnames(ndr.weak.results) = c("alt","rho","E.NDR","M.NDR","S.NDR", "DL.NDR","EX.NDR")
fdr.weak.results[,1] = ndr.weak.results[,1] = mualt
fdr.weak.results[,2] = ndr.weak.results[,2] = rho

#perform simulations
for(i in 1:900)
  {
    #get data
    set.seed(i)
    my.zvals = rbind(mvrnorm(n1,c(mualt,mualt),matrix(c(1,rho[i],rho[i],1),2,2)),mvrnorm(n0,c(0,0),matrix(c(1,rho[i],rho[i],1),2,2)))
    my.data = 2*pnorm(-abs(my.zvals))
    x.values = my.data[,1]; y.values = my.data[,2]

    #get voronoi tessellation and extract cell areas
    areas = deldir(x.values,y.values,rw=c(0,1,0,1),digits=20,eps=1e-13) 
    tess.areas = areas$summary$dir.area 

    #get rankings and indices
    distance.E = apply(my.data,1,function(x){sqrt(x[1]^2+x[2]^2)}); rank.E = sort(distance.E,index.return=T)$ix
    distance.M = apply(my.data,1,max) ; rank.M = sort(distance.M,index.return=T)$ix
    distance.S = apply(my.data,1,sum); rank.S = sort(distance.S,index.return=T)$ix
    distance.DL = apply(my.data,1,function(x){prod(x)*(1+(x[1]/.001)^2)*(1+(x[2]/.001)^2)}); rank.DL = sort(distance.DL,index.return=T)$ix

    #get cumulative sums, then myfdr and mypower for each ranking scheme
    sum.E = cumsum(tess.areas[rank.E]); E.results = simtest(sum.E,rank.E,P=Pweak)
    sum.M = cumsum(tess.areas[rank.M]); M.results = simtest(sum.M,rank.M,P=Pweak)
    sum.S = cumsum(tess.areas[rank.S]); S.results = simtest(sum.S,rank.S,P=Pweak)
    sum.DL = cumsum(tess.areas[rank.DL]); DL.results =simtest(sum.DL,rank.DL,P=Pweak)

    #record results
    ndr.weak.results[i,3] = E.results$mypower; fdr.weak.results[i,3] = E.results$myfdr
    ndr.weak.results[i,4]= M.results$mypower; fdr.weak.results[i,4] = M.results$myfdr
    ndr.weak.results[i,5]  = S.results$mypower; fdr.weak.results[i,5]  = S.results$myfdr
    ndr.weak.results[i,6] = DL.results$mypower; fdr.weak.results[i,6] = DL.results$myfdr

    #now get fdr and power from using B-H on the maximums (Comparison)
    max.pvalues = apply(my.data,1,max)
    max.bh = BH(max.pvalues)$index
    k = length(max.bh)
    if(k==0)
      {
      fdr.weak.results[i,7] = ndr.weak.results[i,7] = 0
      }
    if(k>0)
      {
      fdr.weak.results[i,7] = sum(max.bh>n1)/k
      ndr.weak.results[i,7] = sum(max.bh<=n1)/n1
      }
  }

#produce table of results
ndr.weak.means = fdr.weak.means = matrix(NA,9,7)
colnames(fdr.weak.means) = colnames(ndr.weak.means) = c("alt","rho","Euclidean","Maximum","De Lichtenberg","Summation","Existing")
for(i in 0:8)
{
  temp.low = 100*i+1;  temp.high = 100*(i+1)
  ndr.weak.means[i+1,] = apply(ndr.weak.results[temp.low:temp.high,],2,mean)
  fdr.weak.means[i+1,] = apply(fdr.weak.results[temp.low:temp.high,],2,mean)
}

#print tables of results
print.xtable(xtable(fdr.weak.means,digits=3))
print.xtable(xtable(ndr.weak.means,digits=3))


# ---------------# moderate alternative simulation #------------------------- #

mualt <- 3 # set mean of alternative signals
Pmod <- 800 # set penalization factor to pass to mixFdr

#inilialize matrices to store results
fdr.mod.results = ndr.mod.results = matrix(NA,length(rho),7)
colnames(fdr.mod.results) = c("alt","rho","E.FDR","M.FDR","S.FDR","DL.FDR","EX.FDR")
colnames(ndr.mod.results) = c("alt","rho","E.NDR","M.NDR","S.NDR", "DL.NDR","EX.NDR")
fdr.mod.results[,1] = ndr.mod.results[,1] = mualt
fdr.mod.results[,2] = ndr.mod.results[,2] = rho

#perform simulation
for(i in 1:900)
  {
    #get data
    set.seed(i)
    my.zvals = rbind(mvrnorm(n1,c(mualt,mualt),matrix(c(1,rho[i],rho[i],1),2,2)),mvrnorm(n0,c(0,0),matrix(c(1,rho[i],rho[i],1),2,2)))
    my.data = 2*pnorm(-abs(my.zvals))
    x.values = my.data[,1]; y.values = my.data[,2]
    
    #get voronoi tessellation and extract cell areas
    areas = deldir(x.values,y.values,rw=c(0,1,0,1),digits=20,eps=1e-13) 
    tess.areas = areas$summary$dir.area 
    
    #get rankings and indices
    distance.E = apply(my.data,1,function(x){sqrt(x[1]^2+x[2]^2)}); rank.E = sort(distance.E,index.return=T)$ix
    distance.M = apply(my.data,1,max) ; rank.M = sort(distance.M,index.return=T)$ix
    distance.S = apply(my.data,1,sum); rank.S = sort(distance.S,index.return=T)$ix
    distance.DL = apply(my.data,1,function(x){prod(x)*(1+(x[1]/.001)^2)*(1+(x[2]/.001)^2)}); rank.DL = sort(distance.DL,index.return=T)$ix
    
    #get cumulative sums, then myfdr and mypower for each ranking scheme
    sum.E = cumsum(tess.areas[rank.E]); E.results = simtest(sum.E,rank.E,P=Pmod)
    sum.M = cumsum(tess.areas[rank.M]); M.results = simtest(sum.M,rank.M,P=Pmod)
    sum.S = cumsum(tess.areas[rank.S]); S.results = simtest(sum.S,rank.S,P=Pmod)
    sum.DL = cumsum(tess.areas[rank.DL]); DL.results =simtest(sum.DL,rank.DL,P=Pmod)
    
    #record these results
    ndr.mod.results[i,3] = E.results$mypower; fdr.mod.results[i,3] = E.results$myfdr
    ndr.mod.results[i,4]= M.results$mypower; fdr.mod.results[i,4] = M.results$myfdr
    ndr.mod.results[i,5]  = S.results$mypower; fdr.mod.results[i,5]  = S.results$myfdr
    ndr.mod.results[i,6] = DL.results$mypower; fdr.mod.results[i,6] = DL.results$myfdr
  
    #now get fdr and power from using B-H on the maximums (Comparison)
    max.pvalues = apply(my.data,1,max)
    max.bh = BH(max.pvalues)$index
    k = length(max.bh)
    fdr.mod.results[i,7] = ndr.mod.results[i,7] = 0
    if(k>0)
      {
      fdr.mod.results[i,7] = sum(max.bh>n1)/k
      ndr.mod.results[i,7] = sum(max.bh<=n1)/n1
      }
}

# produce table of results
ndr.mod.means = fdr.mod.means = matrix(NA,9,7)
colnames(fdr.mod.means) = colnames(ndr.mod.means) = c("alt","rho","Euclidean","Maximum","De Lichtenberg","Summation","Existing")
for(i in 0:8)
{
  temp.low = 100*i+1;   temp.high = 100*(i+1)
  ndr.mod.means[i+1,] = apply(ndr.mod.results[temp.low:temp.high,],2,mean)
  fdr.mod.means[i+1,] = apply(fdr.mod.results[temp.low:temp.high,],2,mean)
}

# print tables of results
print.xtable(xtable(fdr.mod.means,digits=3))
print.xtable(xtable(ndr.mod.means,digits=3))

# ---------------------------# Strong alternative simulation #------------------------------------- #

mualt <- 4 #set mean of alternative signals
Pstrong <- 1000# Set penalization factor to pass to mixFdr

#initialize matrices to store results
fdr.strong.results = ndr.strong.results = matrix(NA,length(rho),7)
colnames(fdr.strong.results) = c("alt","rho","E.FDR","M.FDR","S.FDR","DL.FDR","EX.FDR")
colnames(ndr.strong.results) = c("alt","rho","E.NDR","M.NDR","S.NDR", "DL.NDR","EX.NDR")
fdr.strong.results[,1] = ndr.strong.results[,1] = mualt
fdr.strong.results[,2] = ndr.strong.results[,2] = rho

for(i in 1:900)
  {
    #get data
    set.seed(i)
    my.zvals = rbind(mvrnorm(n1,c(mualt,mualt),matrix(c(1,rho[i],rho[i],1),2,2)),mvrnorm(n0,c(0,0),matrix(c(1,rho[i],rho[i],1),2,2)))
    my.data = 2*pnorm(-abs(my.zvals))
    x.values = my.data[,1]; y.values = my.data[,2]
    
    #get voronoi tessellation and extract cell areas
    areas = deldir(x.values,y.values,rw=c(0,1,0,1),digits=20,eps=1e-13) 
    tess.areas = areas$summary$dir.area 
    
    #get rankings and indices
    distance.E = apply(my.data,1,function(x){sqrt(x[1]^2+x[2]^2)}); rank.E = sort(distance.E,index.return=T)$ix
    distance.M = apply(my.data,1,max) ; rank.M = sort(distance.M,index.return=T)$ix
    distance.S = apply(my.data,1,sum); rank.S = sort(distance.S,index.return=T)$ix
    distance.DL = apply(my.data,1,function(x){prod(x)*(1+(x[1]/.001)^2)*(1+(x[2]/.001)^2)}); rank.DL = sort(distance.DL,index.return=T)$ix
    
    #get cumulative sums, then myfdr and mypower for each ranking scheme
    sum.E = cumsum(tess.areas[rank.E]); E.results = simtest(sum.E,rank.E,P=Pstrong)
    sum.M = cumsum(tess.areas[rank.M]); M.results = simtest(sum.M,rank.M,P=Pstrong)
    sum.S = cumsum(tess.areas[rank.S]); S.results = simtest(sum.S,rank.S,P=Pstrong)
    sum.DL = cumsum(tess.areas[rank.DL]); DL.results =simtest(sum.DL,rank.DL,P=Pstrong)
    
    #record results
    ndr.strong.results[i,3] = E.results$mypower; fdr.strong.results[i,3] = E.results$myfdr
    ndr.strong.results[i,4]= M.results$mypower; fdr.strong.results[i,4] = M.results$myfdr
    ndr.strong.results[i,5]  = S.results$mypower; fdr.strong.results[i,5]  = S.results$myfdr
    ndr.strong.results[i,6] = DL.results$mypower; fdr.strong.results[i,6] = DL.results$myfdr
    
    #now get fdr and power from using B-H on the maximums (Comparison)
    max.pvalues = apply(my.data,1,max)
    max.bh = BH(max.pvalues)$index
    k = length(max.bh)
    fdr.strong.results[i,7] = ndr.strong.results[i,7] = 0
    if(k>0)
      {
        fdr.strong.results[i,7] = sum(max.bh>n1)/k
        ndr.strong.results[i,7] = sum(max.bh<=n1)/n1
      }
}

#get a table of results for strong
ndr.strong.means = fdr.strong.means = matrix(NA,9,7)
colnames(fdr.strong.means) = colnames(ndr.strong.means) = c("alt","rho","Euclidean","Maximum","De Lichtenberg","Summation","Existing")

for(i in 0:8)
  {
    temp.low = 100*i+1;   temp.high = 100*(i+1)
    ndr.strong.means[i+1,] = apply(ndr.strong.results[temp.low:temp.high,],2,mean)
    fdr.strong.means[i+1,] = apply(fdr.strong.results[temp.low:temp.high,],2,mean)
  }

#print results
print.xtable(xtable(fdr.strongd.means,digits=3))
print.xtable(xtable(ndr.strong.means,digits=3))

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
rho = rep(seq(from=0,to=.8,by=.1),each=100)
n = 2000 #number of genes for each data set
p = .1 # number of true alternative signals (alt,alt)
p1 = .1 #proportion of p-vectors (alt, null)
p2 = .1 #proportion of p-vectors (null, alt)
n1 = p*n #number of p-vectors (alt, alt)
n11 = p1*n #number of p-vectors (alt, null)
n12 = p2*n #number of p-vectors (null, alt)
n0 = n*(1-p-p1-p2) # number of null p-vectors
mualt = 3 #mean of alternative signal

# --------------------------# SUPPLEMENTARY 1 #------------------------------------- #
# Simulation where mualt=3, and 10% are halfnull  (3,0), 10% are halfnull (0,3) #
 
#initialize matrices to store results
fdr.test.results = ndr.test.results = matrix(NA,900,6)
all.null.results = half.null.results = matrix(NA,900,6)
colnames(fdr.test.results) = c("rho","M.FDR","E.FDR","S.FDR","DL.FDR","Ex.FDR")
colnames(ndr.test.results) = c("rho","M.NDR","E.NDR","S.NDR", "DL.NDR","Ex.NDR")
colnames(all.null.results) = c("rho","M.allFDR","E.allFDR","S.allFDR","DL.allFDR","Ex.allFDR")
colnames(half.null.results) = c("rho","M.halfFDR","E.halfFDR","S.halfFDR", "DL.halfFDR","Ex.halfFDR")
fdr.test.results[,1] = ndr.test.results[,1] = rho
all.null.results[,1] = half.null.results[,1] = rho

PMod=800

for(i in 1:900)
  {
    #get data
    set.seed(i)
    if(i==64)
    {
      set.seed(128) # something strange happens when i=64
    }
    r = rho[i]
    my.zvals = rbind(mvrnorm(n1,c(mualt,mualt),matrix(c(1,r,r,1),2,2)),
                      mvrnorm(n11,c(mualt,0),matrix(c(1,r,r,1),2,2)),
                      mvrnorm(n12,c(0,mualt),matrix(c(1,r,r,1),2,2)),
                      mvrnorm(n0,c(0,0),matrix(c(1,r,r,1),2,2)))
    my.data = 2*pnorm(-abs(my.zvals))
    x.values = my.data[,1]; y.values = my.data[,2]
    
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
fdr.means = ndr.means = all.null.fdr.means = half.null.fdr.means = matrix(NA,9,6)
colnames(fdr.means) = colnames(ndr.means) = c("rho","Maximum","Euclidean","Summation","De Lichtenberg","Existing")
colnames(all.null.fdr.means) = colnames(half.null.fdr.means) = c("rho","Maximum","Euclidean","Summation","De Lichtenberg","Existing")

for(i in 0:8){
  temp.low = 100*i+1;   temp.high = 100*(i+1)
  ndr.means[i+1,] = apply(ndr.test.results[temp.low:temp.high,],2,mean)
  fdr.means[i+1,] = apply(fdr.test.results[temp.low:temp.high,],2,mean)
  all.null.fdr.means[i+1,] = apply(all.null.results[temp.low:temp.high,],2,mean)
  half.null.fdr.means[i+1,] = apply(half.null.results[temp.low:temp.high,],2,mean)
}

#print the tables
print.xtable(xtable(fdr.means,digits=3))
print.xtable(xtable(ndr.test.means,digits=3))



## --------------------------# SUPPLEMENTARY 2 #------------------------------------- #
#  Vary the proportion of 'half-null' hypotheses.  Make them challenging.
# i.e. - mualt=(3,3), and halfnull 1 - (0,4), halfnull 2 - (4,0), letting rho=0
##----------------------------------------------------------------------------------- ##

#preliminaries
halfnull1 <- rep(seq(from=.02,to=.14,by=.02),each=100)
halfnull2 <- rep(seq(from=0.02,to=.14,by=.02),each=100)
PMod <- 800
r <- 0 

#allocate space to store results
fdr.test.results = ndr.test.results = matrix(NA,700,6)
colnames(fdr.test.results) = c("propHalfNull","M.FDR","E.FDR","S.FDR","DL.FDR","Ex.FDR")
colnames(ndr.test.results) = c("propHalfNull","M.NDR","E.NDR","S.NDR", "DL.NDR","Ex.NDR")
fdr.test.results[,1] = ndr.test.results[,1] = 2*halfnull1

set.seed(4445)
# simulation
for(i in 1:700)
  {
    if(i==227)
    {
      i= 2*i #something strange happens when i=227
    }
    n11 = halfnull1[i]*n
    n12 = halfnull2[i]*n
    n0 = n - n1 - n11 - n12
    #get data
    my.zvals = rbind(matrix(rnorm(2*n1,mean=3),n1,2),
                      cbind(rnorm(n11,mean=4),rnorm(n11,mean=0)),
                      cbind(rnorm(n11,mean=0),rnorm(n11,mean=4)),
                      matrix(rnorm(2*n0),n0,2))
    my.data = 2*pnorm(-abs(my.zvals))
    x.values = 20*my.data[,1]; y.values = 10*my.data[,2]
  
    #get voronoi tessellation and extract cell areas
    areas = deldir(x.values,y.values,rw=c(0,10,0,10),digits=20,eps=1e-13) 
    tess.areas = areas$summary$dir.area/400 
  
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
  
    #record results  
    ndr.test.results[i,3] = E.results$mypower; fdr.test.results[i,3] = E.results$myfdr
    ndr.test.results[i,2]= M.results$mypower; fdr.test.results[i,2] = M.results$myfdr
    ndr.test.results[i,4]  = S.results$mypower; fdr.test.results[i,4]  = S.results$myfdr
    ndr.test.results[i,5] = DL.results$mypower; fdr.test.results[i,5] = DL.results$myfdr
    
    #try out competitors method
    max.pvalues = apply(my.data,1,max)
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
  }
}

# get table of means
fdr.means = ndr.means = all.null.fdr.means = half.null.fdr.means = matrix(NA,7,6)
colnames(fdr.means) = colnames(ndr.means) = c("rho","Maximum","Euclidean","Summation","De Lichtenberg","Existing")
colnames(all.null.fdr.means) = colnames(half.null.fdr.means) = c("rho","Maximum","Euclidean","Summation","De Lichtenberg","Existing")

for(i in 0:6)
{
  temp.low = 100*i+1;   temp.high = 100*(i+1)
  ndr.means[i+1,] = apply(ndr.test.results[temp.low:temp.high,],2,mean)
  fdr.means[i+1,] = apply(fdr.test.results[temp.low:temp.high,],2,mean)
  all.null.fdr.means[i+1,] = apply(all.null.results[temp.low:temp.high,],2,mean)
  half.null.fdr.means[i+1,] = apply(half.null.results[temp.low:temp.high,],2,mean)
}

#print table of means
print.xtable(xtable(fdr.means,digits=3))
print.xtable(xtable(ndr.test.means,digits=3))
