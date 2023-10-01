index <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

library(glmnet)
library(KRLS)
library(xtable)
library(vegan)

setwd("/home/nma/Desktop/klda/20210618")
source("source0618.R")

autocorr.mat <- function(p = 100, rho = 0.9) {
   mat <- diag(p)
   return(rho^abs(row(mat)-col(mat)))
}

cent.bray <- function(Ky){
	n = dim(Ky)[1]
	In = diag(rep(1, n))
	K1 = rm.zero(Ky)
	K1 = (In - 1/n) %*% K1 %*% K1 %*% (In-1/n)
	return(K1)
}

rm.zero <- function(K){
  K[is.na(K)] <- 0.01
  return(K)
}

compare = function(x1,y1,z1,test.x,test.y,test.z, kernel = 1){
	# linear
	if(kernel == 1){
		H1 = z1
	test.H = test.z
}else if(kernel==2) {
	H1 = z1 %*% t(z1)
		test.H = test.z %*% t(test.z)
}else if(kernel==3){
		H1 = gausskernel(z1,sigma=1)
	test.H = gausskernel(test.z, sigma=1)
}else {
	H1 = as.matrix(vegdist(abs(z1),method="bray"))
		test.H = as.matrix(vegdist(abs(test.z),method="bray"))
		H1 = cent.bray(H1)
		H1 = scale(H1)
		test.H = cent.bray(test.H)
		test.H = scale(test.H)
}

p=30
	
	
		fit.os.l1 <- sklda.os.sparse(x1,y1,H1,)
	Y1 = y1[,1]

	
fit0.0 = lda(Y1~x1)
pd0.0 = predict(fit0.0)
lda.train=match.beta(Y1, pd0.0$class)
pd0.0 = predict(fit0.0, data.frame(test.x))
lda.test=match.beta(test.y[,1], pd0.0$class)



fit.os.old2 <- sklda.os.old2(x1,y1,H1,Omega=diag(rep(1, p)))


pd.train.old21 = pred.os(fit.os.old2,x1, H1, Y1)
#pd.train.old22 = pred.os2(fit.os.old2,x1, H1, Y1)
pd.train.old23 = pred.os3(fit.os.old2,x1, H1, Y1)
#pd.train.old24 = pred.os4(fit.os.old2,x1, H1, Y1)

rst.train.old21=match.beta(Y1, pd.train.old21$class)
#rst.train.old22=match.beta(Y1, pd.train.old22$class)
rst.train.old23=match.beta(Y1, pd.train.old23$class)
#rst.train.old24=match.beta(Y1, pd.train.old24$class)


pd.test.old21 = pred.os(fit.os.old2,test.x, test.H, test.y[,1])
#pd.test.old22 = pred.os2(fit.os.old2,test.x, test.H, test.y[,1])
pd.test.old23 = pred.os3(fit.os.old2,test.x, test.H, test.y[,1])
#pd.test.old24 = pred.os4(fit.os.old2,test.x, test.H, test.y[,1])

rst.test.old21=match.beta(test.y[,1], pd.test.old21$class)
#rst.test.old22=match.beta(test.y[,1], pd.test.old22$class)
rst.test.old23=match.beta(test.y[,1], pd.test.old23$class)
#rst.test.old24=match.beta(test.y[,1], pd.test.old24$class)



fit2 = sda(x1,y1,lambda = 1e-6,stop = -1)
pd2 = predict(fit2, data.frame(x1))
sda.train2=match.beta(Y1, pd2$class)
pd2.1 = predict(fit2, data.frame(test.x))
sda.test2=match.beta(test.y[,1], pd2.1$class)




M1 = cbind(x1,z1)
fit3 = sda(M1,y1,lambda = 1e-6,stop = -1)
pd3 = predict(fit3, data.frame(M1))
sda.train3=match.beta(Y1, pd3$class)
pd3.1 = predict(fit3, data.frame(cbind(test.x, test.z)))
sda.test3=match.beta(test.y[,1], pd3.1$class)



dta.kfda = data.frame(cbind(x1,y1[,1]))
fit.1 = kfda(dta.kfda,kernel.name = "rbfdot")
pd1 = kfda.predict(fit.1, data.frame(test.x))
pd1.train = kfda.predict(fit.1, data.frame(x1))
kda.test = match.beta(test.y[,1], pd1$class)
kda.train = match.beta(Y1, pd1.train$class)


dta.kfda2 = data.frame(cbind(x1,z1,y1[,1]))
fit.2 = kfda(dta.kfda2,kernel.name = "rbfdot")
pd2 = kfda.predict(fit.2, data.frame(cbind(test.x, test.z)))
pd2.train = kfda.predict(fit.2, data.frame(cbind(x1,z1)))
kda.test2 = match.beta(test.y[,1], pd2$class)
kda.train2 = match.beta(Y1, pd2.train$class)	
	
	


pd.train.l11 = pred.os(fit.os.l1, x1, H1, Y1)
#pd.train.old22 = pred.os2(fit.os.old2,x1, H1, Y1)
pd.train.l12 = pred.os3(fit.os.l1,x1, H1, Y1)
#pd.train.old24 = pred.os4(fit.os.old2,x1, H1, Y1)

rst.train.l11=match.beta(Y1, pd.train.l11$class)
#rst.train.old22=match.beta(Y1, pd.train.old22$class)
rst.train.l12=match.beta(Y1, pd.train.l12$class)
#rst.train.old24=match.beta(Y1, pd.train.old24$class)


pd.test.l11 = pred.os(fit.os.l1,test.x, test.H, test.y[,1])
#pd.test.old22 = pred.os2(fit.os.old2,test.x, test.H, test.y[,1])
pd.test.l12 = pred.os3(fit.os.l1,test.x, test.H, test.y[,1])
#pd.test.old24 = pred.os4(fit.os.old2,test.x, test.H, test.y[,1])

rst.test.l11=match.beta(test.y[,1], pd.test.l11$class)
#rst.test.old22=match.beta(test.y[,1], pd.test.old22$class)
rst.test.l12=match.beta(test.y[,1], pd.test.l12$class)
#rst.test.old24=match.beta(test.y[,1], pd.test.old24$class)

cvfit1 = cv.glmnet(x1, y1[,1], family = "binomial", type.measure = "class")
logit.lasso = glmnet(x1, y1[,1], family="binomial", lambda = cvfit1$lambda.min)
cvfit2 = cv.glmnet(x1, y1[,1], family = "binomial", type.measure = "class")
logit.lasso2 = glmnet(cbind(x1,z1), y1[,1], family="binomial", lambda = cvfit1$lambda.min)


pd.logit1.train = predict(logit.lasso, x1, type="class")
rst.train.logit1 = match.beta(test.y[,1], as.numeric(as.factor(pd.logit1.train)))

pd.logit1.test = predict(logit.lasso, test.x, type="class")
rst.test.logit1 = match.beta(test.y[,1], as.numeric(as.factor(pd.logit1.test)))

pd.logit2.train = predict(logit.lasso2, cbind(x1,z1), type="class")
rst.train.logit2 = match.beta(test.y[,1], as.numeric(as.factor(pd.logit2.train)))

pd.logit2.test = predict(logit.lasso2, cbind(test.x, test.z), type="class")
rst.test.logit2 = match.beta(test.y[,1], as.numeric(as.factor(pd.logit2.test)))

rst = c(
	rst.train.old21,

rst.train.old23,

lda.train,
sda.train2,
sda.train3,
kda.train,
kda.train2,

rst.train.l11,

rst.train.l12,
rst.train.logit1,
rst.train.logit2,
rst.test.old21,

rst.test.old23,

lda.test,
sda.test2,
sda.test3,
kda.test,
kda.test2,

rst.test.l11,

rst.test.l12,
rst.test.logit1,
rst.test.logit2

)
}


do.one <- function(kernel = 1){

n = 100
p = 30
r = 50
m = 4

labels = rep(c(0,1), each=n/2)
X = matrix(rnorm(n*p), ncol=p, nrow=n)
# X[,1] = X[,1] + rep(c(0,1), each=n/2)
err1 = matrix(0, n, n)
cpd1 = 0.2 + diag(rep(0.3,m))

for(i in 1:(n/m)){
	ind = (1:m)+(i-1)*m
	#print(ind)
err1[ind,ind] = cpd1
}

err2 = err1
err3 = matrix(rnorm(n*n, mean=0, sd=0.2), nrow=n,ncol=n)
err4 = matrix(rnorm(n*n, mean=0, sd=0.2), nrow=n,ncol=n)


x1 = X
zz1 = err1 + err3
# dd1 = svd(zz1)
# z1 = dd1$u%*% diag(sqrt(abs(dd1$d)))
z1 = zz1


Beta = matrix(rnorm(n*p), ncol=n, nrow=p)
dd1 = svd(err1)
z1 = x1 %*% Beta + dd1$u%*% diag(sqrt(abs(dd1$d)))



b1 = rep(c(1,0), each=p/2)
b2 = rep(c(0,1), each=n/2)

# b1 = rep(1, each=p)

b1 = rep(1, each=p)
b2 = rep(1, each=n)


linmod =  x1 %*% b1 #+ z1 %*% b2
pr = 1/(1+exp(-linmod))
y0 = rbinom(n,size = 1, prob=pr)
y1 = cbind(y0, 1-y0)

y1 = y1[1:100,]
z1 = z1[1:100,1:100]
x1 = x1[1:100,]
z1[,1:5] <- z1[,1:5] + y1[,1]
z1[,6:10] <- z1[,6:10] + y1[,2]
x1[,1:5] <- x1[,1:5] + y1[,1]
x1[,6:10] <- x1[,6:10] + y1[,2]

Xnew = matrix(rnorm(n*p), ncol=p, nrow=n)
# Xnew[,1] = labels + Xnew[,1]
Xnew = x1 + matrix(rnorm(n*p), ncol=p, nrow=n)[1:100,]
# Xnew = matrix(rnorm(n*p), ncol=p, nrow=n)[1:36,]
# X[,1] = X[,1] + rep(c(0,1), each=n/2)
order = sample(n)

Beta2 = matrix(rnorm(n*p), ncol=n, nrow=p)
test.x = Xnew[order,]
test.z = test.x %*% Beta2
test.z = test.z[1:100,1:100]
test.y = y1[order,]

test.z[,1:5] <- test.z[,1:5] + test.y[,1]
test.z[,6:10] <- test.z[,6:10] + test.y[,2]
test.x[,1:5] <- test.x[,1:5] + test.y[,1]
test.x[,6:10] <- test.x[,6:10] + test.y[,2]
# test.x = Xnew[order,]
# test.z = err2[order,]+err4[order,]
# test.z = err2[order,order]+err4[order,order]


rst = compare(x1,y1,z1,test.x,test.y,test.z, kernel = kernel)
return(rst)
}

# result = replicate(1,do.one(kernel=3))
result = replicate(100,do.one(kernel=index))
rownames(result) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2",
	"klda.test","klda2.test","lda.test","sda.test","sda2.test","kda.test","kda2.test","sklda.test","sklda2.test","logit.lasso1.test", "logit.lasso2.test")

write.csv(result, paste("simu.0618.kernelcounpound100", index, ".csv", sep=""))




