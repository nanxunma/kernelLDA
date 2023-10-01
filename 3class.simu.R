#index <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

library(ape)
library(phyloseq)
library(MASS)
library(sparseLDA)
setwd("/home/nma/Desktop/klda/0720")
source("sklda.source.0720.R")

 #ml R/4.0.2-foss-2019b


nSim=1000

do.one <- function(n=120){
p=32
labels = c(rep(0,n/3), rep(1,n/3), rep(2,n/3))
nc1 = nc2 = nc3 = n/3



smp1 = matrix(0, ncol=p, nrow=nc1/2)
smp2 = matrix(1, ncol=p, nrow=nc1/2)
sample = rbind(smp1,smp2,smp1,smp2,smp1,smp2)
col1 = rep(c(1,0,0), each=n/3)
col2 = rep(c(0,1,0), each=n/3)
col3 = rep(c(0,0,1), each=n/3)
sample[,16] = col1
#sample[,22] = 1-col1
sample[,1] = col2
#sample[,30] = 1-col2
sample[,7] = col3
#sample[,23]= 1-col3


col.name = rep(NA,p)
for(i in 1:p){
  col.name[i] = paste("t",i, sep = "")
}
colnames(sample) = col.name

tree2 <- read.tree("example.t1.tre")
physq <- phyloseq(otu_table(sample, taxa_are_rows = FALSE), tree2)
plot(tree2)

err = matrix(rnorm(n*p, sd=1),ncol=p)
err2 = matrix(rnorm(n*n, sd=0.1),ncol=n)
y1 = cbind(labels==0, labels==1,labels==2)
Y.1 = as.factor(labels)
x1 = sample + err


uw <- UniFrac(physq)
mat = as.matrix(uw)
mat.sim = nmat(mat*mat)
H1 = mat.sim + err2



fit.os <- sklda.os2(x1,Y.1,H1,Omega=0)

#fit.ray <- sklda.rayleigh2(x1,y1,diag(rep(0,n)),Omega=0)
fit.ray <- sklda.rayleigh2(x1,y1,H1,Omega=0)

#fit.os$beta/fit.ray
#fit.os$beta/fit.ray[,c(2:1)]

##U%*%x%*%beta
pd1 = pred.os(fit.os, x1, H1, Y.1)

klda.train1=match.3class(Y.1, pd1$class)

pd3 = pred.os3(fit.os, x1, H1, Y.1)

klda.train3=match.3class(Y.1, pd3$class)

#x3 = H1%*%x1
#fit0.1 = lda(Y.1~x3)
#pd0.1 = predict(fit0.1)
#lda.kx.train=match.3class(Y.1, pd0.1$class)







test.err = matrix(rnorm(n*p, sd=1),ncol=p)
test.x = sample + test.err
err4 = matrix(rnorm(n*n, sd=0.1),ncol=n)

order = sample(n)
test.y = y1[order,]
test.x = test.x[order,]
test.y1 = Y.1[order]

test.physq <- phyloseq(otu_table(sample[order,], taxa_are_rows = FALSE), tree2)
test.uw <- UniFrac(test.physq)
test.mat = as.matrix(test.uw)
test.mat.sim = nmat(test.mat*test.mat)
test.H = test.mat.sim + err4



pd1.test = pred.os(fit.os, test.x, test.H, Y.1)
klda.test1=match.3class(test.y1, pd1.test$class)

pd3.test = pred.os3(fit.os, test.x, test.H, Y.1)
klda.test3=match.3class(test.y1, pd3.test$class)

fit0.0 = lda(Y.1~x1)
pd0.0 = predict(fit0.0)
lda.train=match.3class(Y.1, pd0.0$class)
pd0.0 = predict(fit0.0, data.frame(test.x))
lda.test=match.3class(test.y1, pd0.0$class)



fit1 = sda(x1, Y.1)
pd1 = predict(fit1, x1)
sda.train=match.3class(Y.1, pd1$class)
pd1 = predict(fit1, data.frame(test.x))
sda.test=match.3class(test.y1, pd1$class)


rst = c(klda.train1, klda.test1,klda.train3, klda.test3, lda.train, lda.test, sda.train, sda.test)

return(rst)
}

result = replicate(nSim,do.one(n=120))
write.csv(result, "simu.3class.csv")






