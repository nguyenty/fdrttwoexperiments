# R code for prostate cancer data
# R code for Arul's data

source("C:/Users/Daisy Philtron/Documents/PSU/Research/Voronoi Combination/R work on 2d testing/disjunction_combination.r")

ge <- read.delim("C:/Users/Daisy Philtron/Documents/PSU/Research/Voronoi Combination/Data/Cancer Data/expr_HPCA",header=T,skip=1)
copy <- read.delim("C:/Users/Daisy Philtron/Documents/PSU/Research/Voronoi Combination/Data/Cancer Data/copy_HPCA",header=T,skip=1)

# Compute T-statistics for both platforms
ge.tstat <- apply(ge[,31:47],1,mean,na.rm=T)-apply(ge[,4:30],1,mean,na.rm=T)/sqrt(apply(ge[,4:47],1,var,na.rm=T)*(1/17+1/27))
copy.tstat <- apply(copy[,31:47],1,mean,na.rm=T)-apply(copy[,4:30],1,mean,na.rm=T)/sqrt(apply(copy[,4:47],1,var,na.rm=T)*(1/17+1/27))

# Find common genes on both platforms
common = intersect(ge$UNIGENE,copy$UNIGENE)
genes = pmatch(common,ge$UNIGENE)
copy.genes = pmatch(common,copy$UNIGENE)

#get t-statistics for those common genes only
ge.tstat2 = ge.tstat[genes]
#get 2-sided p-values
ge.pvals = 2*(1-pnorm(abs(ge.tstat2)))
copy.tstat2 = copy.tstat[genes]
#get 2-sided copy pvals
copy.pvals = 2*(1-pnorm(abs(copy.tstat2)))

pval.mat = cbind(ge.pvals,copy.pvals)

# -----------------------------------------------------------------------------#
# Do analysis and exploration on data sets individually

x <- seq(from=-5,to=5,by=.01)
#get figure with histogram of p-values with emp. and actual null
postscript(file="C:/Users/Daisy Philton/Documents/PSU/Research/Spring 2013/Figures/tstatsge.eps",onefile=F,horizontal=F,width=5,height=5,paper="special")
par(mar=c(2.1,4.1,1,1))
hist(ge.tstat2,freq=F,ylim=c(0,.4),xlim=c(-5,5),breaks=25,main=" ",xlab=list(cex=1.5,"Gene expression t-statistic"))
lines(x,dnorm(x,0,sqrt(1.93)))
lines(x,dnorm(x),lty=2)
dev.off()


postscript(file="C:/Users/Daisy Philton/Documents/PSU/Research/Semester Weekly Reports/Spring 2013/Figures/tstatscopy.eps",onefile=F,horizontal=F,width=5,height=5,paper="special")
par(mar=c(2.1,4.1,1,1))
hist(copy.tstat2,freq=F,xlim=c(-5,5),ylim=c(0,.4),breaks=25,main=" ",xlab=list(cex=1.5,"Copy number t-statistic"))
lines(x,dnorm(x,0,sqrt(1.52)))
lines(x,dnorm(x),lty=2)
dev.off()
par(mfrow=c(1,1))

postscript(file="C:/Users/Daisy Philton/Documents/PSU/Research/Semester Weekly Reports/Spring 2013/Figures/gepvals.eps",onefile=F,horizontal=F,width=5,height=5,paper="special")
par(mar=c(2.1,4.1,1,1))
hist(ge.pvals,freq=F,breaks=25, main=" ",ylim=c(0,3.5))#ylim=c(0,.4),xlim=c(-5,5),breaks=25,main=" ",xlab=list(cex=1.5,"Gene expression t-statistic"))
dev.off()

postscript(file="C:/Users/Daisy Philton/Documents/PSU/Research/Spring 2013/Figures/copypvals.eps",onefile=F,horizontal=F,width=5,height=5,paper="special")
par(mar=c(2.1,4.1,1,1))
hist(copy.pvals,freq=F,breaks=25, main=" ",ylim=c(0,3.5))#ylim=c(0,.4),xlim=c(-5,5),breaks=25,main=" ",xlab=list(cex=1.5,"Gene expression t-statistic"))
dev.off()

ge.bh <- BH(ge.pvals)
ge.bh$k
copy.bh <- BH(copy.pvals)
copy.bh$k

#0------------------------------------------------------------------
# do analysis on data together

#get the Tessellation
x.values <- ge.pvals; y.values <- copy.pvals

postscript(file="C:/Users/Daisy Philton/Documents/PSU/Research/Semester Weekly Reports/Spring 2013/Figures/copyegpvectors.eps",onefile=F,horizontal=F,width=5,height=5,paper="special")
par(mar=c(4.1,4.1,1,1))
plot(x.values,y.values,main=" ",xlab=list(cex=1.5,"Gene expression p-value"),ylab=list(cex=1.5,"Copy number p-value"),cex=.5,pch=19)
dev.off()

areas <- deldir(x.values,y.values,rw=c(0,1,0,1),digits=20,eps=1e-13) #Get Voronoi tessellation
tess.areas <- areas$summary$dir.area
n <- length(x.values)
distance.E <- apply(pval.mat,1,function(x){sqrt(x[1]^2+x[2]^2)})
distance.M <- apply(pval.mat,1,max)
distance.S <- apply(pval.mat,1,sum)
distance.DL <- apply(pval.mat,1,function(x){prod(x)*(1+(x[1]/.001)^2)*(1+(x[2]/.001)^2)})

rank.E <- sort(distance.E,index.return=T)$ix
rank.M <- sort(distance.M,index.return=T)$ix
rank.S <- sort(distance.S,index.return=T)$ix
rank.DL <- sort(distance.DL,index.return=T)$ix

csum.E <- cumsum(tess.areas[rank.E])
csum.M <- cumsum(tess.areas[rank.M])
csum.S <- cumsum(tess.areas[rank.S])
csum.DL <- cumsum(tess.areas[rank.DL])

zvals.E <- qnorm(csum.E[-n])
zvals.M <- qnorm(csum.M[-n])
zvals.S <- qnorm(csum.S[-n])
zvals.DL <- qnorm(csum.DL[-n])

postscript(file="C:/Users/Daisy Philton/Documents/PSU/Research/Semester Weekly Reports/Spring 2013/Figures/transcumareas.eps",onefile=F,horizontal=F,width=5,height=5,paper="special")
par(mar=c(2.1,4.1,1,1))
hist(zvals.E,freq=F,main=" ",ylim=c(0,.4),xlab=" ")
lines(x,dnorm(x),lty=2)
lines(x,dnorm(x,-.38,1.09))
dev.off()

#First, apply BH to combined values
BH(csum.E)$k # 14
BH(csum.M)$k #12
BH(csum.S)$k #25
BH(csum.DL)$k #26

#higher alpha for using DAVID
BH(csum.E,alpha=.1)$k #82
BH(csum.M,alpha=.1)$k #68
BH(csum.S,alpha=.1)$k #95
BH(csum.DL,alpha=.1)$k #124

#much higher alpha
BH(csum.E,alpha=.2)$k #241
BH(csum.M,alpha=.2)$k #226
BH(csum.S,alpha=.2)$k #306
BH(csum.DL,alpha=.2)$k #743

sig.index <- rank.S[1:25]
sig.index1 <- rank.S[1:95]

#apply existing method - BH on maximum values
pvec.max <- apply(pval.mat,1,max)
result.ex <- BH(pvec.max)$k  #we get ZERO!!!


#Apply mixFdr to combined values and use empirical null
empnull.E <- mixFdr(zvals.E,J=2,cal=F,P=3380)
empnull.M <- mixFdr(zvals.M,J=2,cal=F, P=3766)
empnull.S <- mixFdr(zvals.S,J=2,cal=F, P=3766)
empnull.DL <- mixFdr(zvals.DL,J=2)

#see how many Left-tail FDR are less then .05... NONE
sum(empnull.E$FDRLeft<.05)
sum(empnull.M$FDRLeft<.05)
sum(empnull.S$FDRLeft<.05)
sum(empnull.DL$FDRLeft<.05)

# apply mixFdr using theoretical null
theonull.E <- mixFdr(zvals.E,J=2,theonull=T)
theonull.M <- mixFdr(zvals.M,J=2,theonull=T)
theonull.S <- mixFdr(zvals.S,J=2,theonull=T)
theonull.DL <- mixFdr(zvals.DL,J=2,theonull=T)

#now use Left-tail FDR from theoretical null...
sum(theonull.E$FDRLeft<.05)
sum(theonull.M$FDRLeft<.05)
sum(theonull.S$FDRLeft<.05)
sum(theonull.DL$FDRLeft<.05)

# to test the conjunction hypothesis using stouffer's method
p.stouf <- apply(pval.mat,1,function(x){1-pnorm(sum(qnorm(1-x)/sqrt(2)))})
BH(p.stouf)$k

# test conjunction hypothesis using Fisher's method 
p.fish <- apply(pval.mat,1,function(x){1-qchisq(-2*sum(log(x)),4)})
BH(p.fish)$k


#------------------------------------------------------------------------
#There are 25 significant genes.  I need THEIR NAMES.

sig.index <- rank.S[1:25]
copy.sig <- copy.genes[sig.index]
ge.sig <- genes[sig.index]
sig.names <- intersect(copy$UNIGENE[copy.sig],ge$UNIGENE[ge.sig])
#rbind(ge$chr[ge.sig],ge$pos[ge.sig])

#write.table(sig.names,col.names=F,row.names=F,quote=F,file="C:/Users/Daisy Philtron/Documents/PSU/Research/Semester Weekly Reports/Spring 2013/significant UNIGENE alpha05.txt")

#now get 95 significant genes using alpha=.1 and Summation
sig.index <- rank.S[1:95]
copy.sig <- copy.genes[sig.index]
ge.sig <- genes[sig.index]
sig.names <- intersect(copy$UNIGENE[copy.sig],ge$UNIGENE[ge.sig])

#write.table(sig.names,col.names=F,row.names=F,quote=F,file="C:/Users/Daisy Philtron/Documents/PSU/Research/Semester Weekly Reports/Spring 2013/significant UNIGENE alpha1.txt")

#now get 306 significant genes using alpha=.1 and Summation
sig.index <- rank.S[1:306]
copy.sig <- copy.genes[sig.index]
ge.sig <- genes[sig.index]
sig.names <- intersect(copy$UNIGENE[copy.sig],ge$UNIGENE[ge.sig])

#write.table(sig.names,col.names=F,row.names=F,quote=F,file="C:/Users/Daisy Philtron/Documents/PSU/Research/Semester Weekly Reports/Spring 2013/significant UNIGENE alpha2.txt")




# ------------------------------------------------------------------------- #
#  Everything below this point is from Dr. Ghosh's code that I'm not using...

n = dim(pca.pval.mat)[1]
jkstat = NULL
total = symmdisc(pca.pval.mat)
for (i in 1:n) {
jkstat = c(jkstat,symmdisc(pca.pval.mat[-i,]))
cat(i,"\n")

}

deljk = jkstat-total
stdjk = n*(deljk - total/n)/sqrt(total*(1-total))


n = dim(pca.pval.mat)[1]
jkstat2 = NULL
total2 = stardisc(pca.pval.mat)
for (i in 1:n) {
jkstat2 = c(jkstat2,stardisc(pca.pval.mat[-i,]))
cat(i,"\n")

}

deljk2 = jkstat2-total2
stdjk2 = n*(deljk2 - total2/n)/sqrt(total2*(1-total2))


# Do some basic bioinfo manipulations

tab2 = read.table("idconvert.txt",header=F,sep="\t")
gn1 = as.character(tab2[,1])
ug1 = as.character(tab2[,2])
tlg = common[order(tmp3$posterior[,2],decreasing=T)[1:863]]

gn.tlg = NULL
for (i in 1:863) {
   gn.tlg = c(gn.tlg,gn1[pmatch(tlg[i],ug1)])
}
