#index <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

library(ape)
library(phyloseq)
library(MASS)
library(sparseLDA)
setwd("/home/nma/Desktop/klda/0720")
source("sklda.source.0720.R")

 #ml R/4.0.2-foss-2019b


nSim=1000

do.one <- function(n=100){
p=32
labels = c(rep(1,n/2),rep(0,n/2))

smp1 = matrix(0, ncol=p, nrow=n/4)
smp2 = matrix(1, ncol=p, nrow=n/4)
sample = rbind(smp1,smp2,smp1,smp2)
sample[,16] = rep(c(1,0), each=n/2)
sample[,17] = rep(c(0,1), each=n/2)

col.name = rep(NA,p)
for(i in 1:p){
  col.name[i] = paste("t",i, sep = "")
}
colnames(sample) = col.name

tree2 <- read.tree("example.t1.tre")
physq <- phyloseq(otu_table(sample, taxa_are_rows = FALSE), tree2)


err = matrix(rnorm(n*p, sd=1),ncol=p)
err2 = matrix(rnorm(n*n, sd=0.1),ncol=n)
y1 = cbind(labels, 1-labels)
Y.1 = labels
mu = rep(1, n)
b = rep(1,p)
x1 = sample + err
#image(x1)

uw <- UniFrac(physq)
mat = as.matrix(uw)
mat.sim = nmat(mat*mat)
H1 = mat.sim + err2
#H1 = solve(H1)
#H1 = 
H2 = scale(H1,TRUE,FALSE)


fit.os <- sklda.os2(x1,y1,H1,Omega=0)
fit.ray <- sklda.rayleigh2(x1,y1,H1,Omega=0)
#cbind(fit.os$beta, fit.ray, fit.os$beta/fit.ray)

#KLDA

pd3 = pred.os3(fit.os, x1, H1, Y.1)
klda.train=match.beta(Y.1, pd3$class)







test.err = matrix(rnorm(n*p, sd=1),ncol=p)
test.x = Y.1 + sample + test.err
test.x = sample + test.err
err4 = matrix(rnorm(n*n, sd=0.1),ncol=n)

order = sample(100)
test.y = y1[order,]
test.x = test.x[order,]

test.physq <- phyloseq(otu_table(sample[order,], taxa_are_rows = FALSE), tree2)
test.uw <- UniFrac(test.physq)
test.mat = as.matrix(test.uw)
test.mat.sim = nmat(test.mat*test.mat)
test.H = test.mat.sim + err4

pd3.test = pred.os3(fit.os, test.x, test.H, Y.1)
klda.test=match.beta(test.y[,1], pd3.test$class)

fit0.0 = lda(Y.1~x1)
pd0.0 = predict(fit0.0)
lda.train=match.beta(Y.1, pd0.0$class)
pd0.0 = predict(fit0.0, data.frame(test.x))
lda.test=match.beta(test.y[,1], pd0.0$class)



fit1 = sda(x1, y1)
pd1 = predict(fit1, x1)
sda.train=match.beta(Y.1, pd1$class)
pd1 = predict(fit1, data.frame(test.x))
sda.test=match.beta(test.y[,1], pd1$class)


rst = c(klda.train, klda.test, lda.train, lda.test, sda.train, sda.test)

return(rst)
}

result = replicate(nSim,do.one(n=100))
write.csv(result, "simu.2class.csv")






