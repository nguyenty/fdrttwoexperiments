# --------------------------------------------------------------------------------------- #
# The purpose of this code is to perform the proposed Voronoi combination procedure 
# on the Oliva et al. (2005) data sets, both 2-way and 3-way.
# ----------------------------------------------------------------------------------------#

#First, load all the neccessary libraries (and maybe a few extras)
library(deldir)
library(MASS)

#-------------------------------------------------------------------------------------#
# This is the actual function of interest - it takes in 2-d p-vectors and combines
# tham based on the chosen ranking scheme and their voronoi tessellation.
# It outputs the combine areas and the index of their ranks so that you can track down 
# significance...
#-------------------------------------------------------------------------------------#

CombineVoronoi <- function(x.values,y.values,ranking=1){
  
  my.data <- cbind(x.values,y.values)
  #get voronoi tessellation and extract cell areas for each p-vector
  areas <- deldir(x.values,y.values,rw=c(0,1,0,1),digits=20,eps=1e-13) #Get Voronoi tessellation
  tess.areas <- areas$summary$dir.area #extract cell areas 
  tess.areas[tess.areas<0]<- 0 #In rare cases, deldir gives negative cell areas -> BIG PROBLEM
  
  #get distances according to ranking scheme
  if(ranking==1){distance <- apply(my.data,1,function(x){sqrt(x[1]^2+x[2]^2)})} #1: Euclidean
  if(ranking==2){distance <- apply(my.data,1,max)} #2: maximum
  if(ranking==3){distance <- apply(my.data,1,sum)} #3: sum
  if(ranking==4){distance <- apply(my.data,1,function(x){prod(x)*(1+(x[1]/.001)^2)*(1+(x[2]/.001)^2)})} #4: de lichtenberg
  
  mrank <- sort(distance,index.return=T)$ix #rank p-vectors according to the distances
  csum <- cumsum(tess.areas[mrank]) #find cumulative sums
  
  return(list(index=mrank,cumareas=csum)) 
}

#-------------------------------------------------------------------------------------#
#  Load in data
#-------------------------------------------------------------------------------------#

#load in p-value matrix
p.vec <- as.matrix(read.table("C:/Users/Daisy Philtron/Documents/PSU/Research/Voronoi Combination/Data/oliva_fishers_g.txt"))

a.g <- na.omit(p.vec[,1]); b.g <- na.omit(p.vec[,2]); c.g <- na.omit(p.vec[,3]) #get complete p-value sets for each experiment

ab.index <- which(apply(p.vec,1,function(x){sum(is.na(x[-3]))==0})) #get sets of genes that have complete p-values for sets a and b
g.ab.com <- p.vec[ab.index,c(1,2)]
ac.index <- which(apply(p.vec,1,function(x){sum(is.na(x[-2]))==0})) #get sets of genes that have complete p-values for sets a and c
g.ac.com <- p.vec[ac.index,c(1,3)]
bc.index <- which(apply(p.vec,1,function(x){sum(is.na(x[-1]))==0})) #get sets of genes that have complete p-values for sets b and c
g.bc.com <- p.vec[bc.index,c(2,3)]
abc.index <- which(apply(p.vec,1,function(x){sum(is.na(x))==0})) #get sets of genes that have complete p-values for sets a, b, and c
g.abc.com <- p.vec[abc.index,]

#-------------------------------------------------------------------------------------#
# First, perform marginal analysis for each data set separately.  
#-------------------------------------------------------------------------------------#

#ELUTRIATION A
#-------------------------------------------------------------------------------------#
# % p-values under .05
n.a <- length(a.g)
sum(a.g<.05)/n.a #28.5 %

# Number and percent of p-values rejected under B-H at alpha=.05
Bh.na <- BH(a.g)$k   #number rejections
Bh.na/n.a           # % rejections: 17.3%

#Number of p-values rejected using an empirical null and a Left Tail FDR of .05
t.a <- qnorm(a.g)
empnull.a <- mixFdr(t.a, P=325,calibrate=F,J=2)
sum(empnull.a$FDRLeft<.05)  #gets... 5


#ELUTRIATION B
#-------------------------------------------------------------------------------------#
# % p-values under .05
n.b <- length(b.g)
sum(b.g<.05)/n.b # 22.8%

# Number and percent of p-values rejected under B-H at alpha=.05
Bh.nb <- BH(b.g)$k   #number rejections
Bh.nb/n.b           # % rejections: 6.5%

#Number of p-values rejected using an empirical null and a Left Tail FDR of .05
t.b <- qnorm(b.g)
empnull.b <- mixFdr(t.b, P=388.68,calibrate=F,J=2)
sum(empnull.b$FDRLeft<.05)  #gets... 21


#Cdc25
#-------------------------------------------------------------------------------------#
# % p-values under .05
n.c <- length(c.g)
sum(c.g<.05)/n.c # 66.0%

# Number and percent of p-values rejected under B-H at alpha=.05
Bh.nc <- BH(c.g)$k   #number rejections
Bh.nc/n.c           # % rejections: 60.5%

#Number of p-values rejected using an empirical null and a Left Tail FDR of .05
t.c <- qnorm(c.g)
empnull.c <- mixFdr(t.c, P=285.5,calibrate=F,J=2)
sum(empnull.c$FDRLeft<.05)  #gets... 17

#-------------------------------------------------------------------------------------#
# Second, perform analysis for p-vectors using Elutriation a and b.  
#-------------------------------------------------------------------------------------#

n.ab <- length(g.ab.com[,1])

#get combined values
Combo.E <- CombineVoronoi(g.ab.com[,1],g.ab.com[,2],ranking=1)$cumareas
Combo.M <- CombineVoronoi(g.ab.com[,1],g.ab.com[,2],ranking=2)$cumareas
Combo.S <- CombineVoronoi(g.ab.com[,1],g.ab.com[,2],ranking=3)$cumareas
Combo.DL <- CombineVoronoi(g.ab.com[,1],g.ab.com[,2],ranking=4)$cumareas

#First, apply BH to combined values
ab.bh.E <- BH(Combo.E)$k #we get 225!!
ab.bh.M <- BH(Combo.M)$k # 213
ab.bh.S <- BH(Combo.S)$k #249
ab.bh.DL <- BH(Combo.DL)$k #431

#apply existing method - BH on maximum values
ab.max <- apply(g.ab.com,1,max)
ab.ex <- BH(ab.max)$k  #we get 15!!!

t.ab.E <- qnorm(Combo.E[-n])
t.ab.M <- qnorm(Combo.M[-n])
t.ab.S <- qnorm(Combo.S[-n])
t.ab.DL <- qnorm(Combo.DL[-n])

#Apply mixFdr to combined values
empnull.ab.E <- mixFdr(t.ab.E, calibrate=F, P=275,J=2)
empnull.ab.M <- mixFdr(t.ab.M, calibrate=F, P=275,J=2)
empnull.ab.S <- mixFdr(t.ab.S, calibrate=F, P=275,J=2)
empnull.ab.DL <- mixFdr(t.ab.E, calibrate=F, P=275,J=2)

sum(empnull.ab.E$FDRLeft<.05) # 15
sum(empnull.ab.M$FDRLeft<.05) #12
sum(empnull.ab.S$FDRLeft<.05)#11
sum(empnull.ab.DL$FDRLeft<.05)#14

#-------------------------------------------------------------------------------#
# Extension to all 3 data sets
#-------------------------------------------------------------------------------------#

#first, using existing method.
max.abc <- apply(g.abc.com,1,max)
abc.ex <- BH(max.abc)$k  #k=12

#Now, with average cell areas.
p.vec <- g.abc.com

pair.1 <- p.vec[,c(1,2)] #ab
pair.2 <- p.vec[,c(2,3)] #bc
pair.3 <- p.vec[,c(1,3)] #ac

tess.temp.1 <- deldir(pair.1[,1],pair.1[,2],rw=c(0,1,0,1),digits=20,eps=1e-15) #Get Voronoi tessellation of p-vectors
tess.temp.2 <- deldir(pair.2[,1],pair.2[,2]*10^6,rw=c(0,1,0,10^6),digits=20,eps=1e-15)
tess.temp.3 <- deldir(pair.3[,1]*1000,pair.3[,2]*50000,rw=c(0,1000,0,50000),digits=20,eps=1e-15)

tess.areas.1 <- tess.temp.1$summary$dir.area
tess.areas.2 <- tess.temp.2$summary$dir.area/10^6
tess.areas.3 <- tess.temp.3$summary$dir.area/(1000*50000)


distance <- apply(p.vec,1,function(x){sum(x^2)})
if(scheme==2){distance <- apply(p.vec,1,max)}
if(scheme==3){distance <- apply(p.vec,1,sum)}
ordered.index <- sort(distance, index.return=T)$i
ordered.vectors <- p.vec[ordered.index,]


all.areas <- cbind(tess.areas.1,tess.areas.2,tess.areas.3)
average.areas <- apply(all.areas,1,mean)

ordered.av <- average.areas[ordered.index]
cum.areas <- cumsum(ordered.av)

BH.results <- BH(cum.areas)
BH.results$k

z.values <- qnorm(cum.areas[-1574])  
mix.results <- mixFdr(z.values,J=2,cal=F,P=425)
sum(mix.results$FDRLeft<.05) 