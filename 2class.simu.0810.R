#index <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

library(ape)
cat("lib1")
library(phyloseq)
cat("lib2")
library(MASS)
cat("lib3")
library(sparseLDA)
cat("lib4")
setwd("/home/nma/Desktop/klda/0810")
source("sklda.source.0810.R")

#setwd("/Users/nanxunma/Dropbox/research/codes/phyloseq")
#source("/Users/nanxunma/Dropbox/research/simulation/sklda.source.0810.R")
 #ml R/4.0.2-foss-2019b

#install.packages("ape")


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
#H1 = mat.sim 
H2 = scale(H1,TRUE,FALSE)


fit.os <- sklda.os2(x1,y1,H1,Omega=0)
fit.ray <- sklda.rayleigh2(x1,y1,H1,Omega=0)
fit.ray.old1 <- sklda.rayleigh.old1(x1,y1,H1,Omega=0)
fit.ray.old2 <- sklda.rayleigh.old2(x1,y1,H1,Omega=0)

fit.ray.old1 <- sklda.rayleigh.new1(x1,y1,H1,Omega=0)
fit.ray.old2 <- sklda.rayleigh.new2(x1,y1,H1,Omega=0)
fit.ray.new1 <- sklda.rayleigh.new1(x1,y1,H1,Omega=0)
fit.ray.new2 <- sklda.rayleigh.new2(x1,y1,H1,Omega=0)

fit.os.inv <- sklda.os2(x1,y1,solve(H1),Omega=0)
fit.ray.old1.inv <- sklda.rayleigh.old1(x1,y1,solve(H1),Omega=0)
fit.ray.old2.inv <- sklda.rayleigh.old2(x1,y1,solve(H1),Omega=0)

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
#test.H = test.mat.sim 

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


pd.train.1 = pred.os(fit.os,x1, H1, Y.1)
pd.train.2 = pred.os2(fit.os,x1, H1, Y.1)
pd.train.3 = pred.os3(fit.os,x1, H1, Y.1)
pd.train.4 = pred.os4(fit.os,x1, H1, Y.1)

rst.train.1=match.beta(Y.1, pd.train.1$class)
rst.train.2=match.beta(Y.1, pd.train.2$class)
rst.train.3=match.beta(Y.1, pd.train.3$class)
rst.train.4=match.beta(Y.1, pd.train.4$class)


pd.train.old11 = pred.os(fit.ray.old1,x1, H1, Y.1)
pd.train.old12 = pred.os2(fit.ray.old1,x1, H1, Y.1)
pd.train.old13 = pred.os3(fit.ray.old1,x1, H1, Y.1)
pd.train.old14 = pred.os4(fit.ray.old1,x1, H1, Y.1)

rst.train.old11=match.beta(Y.1, pd.train.old11$class)
rst.train.old12=match.beta(Y.1, pd.train.old12$class)
rst.train.old13=match.beta(Y.1, pd.train.old13$class)
rst.train.old14=match.beta(Y.1, pd.train.old14$class)

pd.train.old21 = pred.os(fit.ray.old2,x1, H1, Y.1)
pd.train.old22 = pred.os2(fit.ray.old2,x1, H1, Y.1)
pd.train.old23 = pred.os3(fit.ray.old2,x1, H1, Y.1)
pd.train.old24 = pred.os4(fit.ray.old2,x1, H1, Y.1)

rst.train.old21=match.beta(Y.1, pd.train.old21$class)
rst.train.old22=match.beta(Y.1, pd.train.old22$class)
rst.train.old23=match.beta(Y.1, pd.train.old23$class)
rst.train.old24=match.beta(Y.1, pd.train.old24$class)

pd.test.1 = pred.os(fit.os,test.x, test.H, Y.1)
pd.test.2 = pred.os2(fit.os,test.x, test.H, Y.1)
pd.test.3 = pred.os3(fit.os,test.x, test.H, Y.1)
pd.test.4 = pred.os4(fit.os,test.x, test.H, Y.1)

rst.test.1=match.beta(test.y[,1], pd.test.1$class)
rst.test.2=match.beta(test.y[,1], pd.test.2$class)
rst.test.3=match.beta(test.y[,1], pd.test.3$class)
rst.test.4=match.beta(test.y[,1], pd.test.4$class)

pd.test.old11 = pred.os(fit.ray.old1,test.x, test.H, Y.1)
pd.test.old12 = pred.os2(fit.ray.old1,test.x, test.H, Y.1)
pd.test.old13 = pred.os3(fit.ray.old1,test.x, test.H, Y.1)
pd.test.old14 = pred.os4(fit.ray.old1,test.x, test.H, Y.1)

rst.test.old11=match.beta(test.y[,1], pd.test.old11$class)
rst.test.old12=match.beta(test.y[,1], pd.test.old12$class)
rst.test.old13=match.beta(test.y[,1], pd.test.old13$class)
rst.test.old14=match.beta(test.y[,1], pd.test.old14$class)


pd.test.old21 = pred.os(fit.ray.old2,test.x, test.H, Y.1)
pd.test.old22 = pred.os2(fit.ray.old2,test.x, test.H, Y.1)
pd.test.old23 = pred.os3(fit.ray.old2,test.x, test.H, Y.1)
pd.test.old24 = pred.os4(fit.ray.old2,test.x, test.H, Y.1)

rst.test.old21=match.beta(test.y[,1], pd.test.old21$class)
rst.test.old22=match.beta(test.y[,1], pd.test.old22$class)
rst.test.old23=match.beta(test.y[,1], pd.test.old23$class)
rst.test.old24=match.beta(test.y[,1], pd.test.old24$class)

cat(".")



rst = c(
rst.train.1,
rst.train.2,
rst.train.3,
rst.train.4,
rst.train.old11,
rst.train.old12,
rst.train.old13,
rst.train.old14,
rst.train.old21,
rst.train.old22,
rst.train.old23,
rst.train.old24,
rst.test.1,
rst.test.2,
rst.test.3,
rst.test.4,
rst.test.old11,
rst.test.old12,
rst.test.old13,
rst.test.old14,
rst.test.old21,
rst.test.old22,
rst.test.old23,
rst.test.old24
)

return(rst)
}

result = replicate(nSim,do.one(n=100))
write.csv(result, "simu.2class.0810.2new.pos.csv")






