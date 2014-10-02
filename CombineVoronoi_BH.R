library(deldir)
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
