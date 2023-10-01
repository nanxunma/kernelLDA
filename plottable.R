getwd()
setwd("/Users/nanxunma/Dropbox/research/simulation/20210618/conpound")


index=3
name = paste("simu.0618.kernel.compound", index, "sparse.csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst3 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst3) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst3

index=4
name = paste("simu.0618.kernel.compound", index, "sparse.csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst4 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst4) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst4[,c("klda")]

rst4test <- rst4[,c("klda", "klda2", "sklda", "sklda2")]

index=2
name = paste("simu.0618.kernel.compound", index, "sparse.csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst2 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst2) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst2[,c("klda")]

rst2test <- rst2[,c("klda","klda2", "sklda", "sklda2")]

xtable(t(cbind(rst3, rst2test,rst4test)),digit=3)

m1 <- t(cbind(rst3, rst2test,rst4test))[,2]
########
########
########
index=3
name = paste("simu.0618.kernelconpound100", index, "sparse.csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst3 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst3) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst3

index=4
name = paste("simu.0618.kernelconpound100", index, "sparse.csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst4 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst4) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst4[,c("klda")]

rst4test <- rst4[,c("klda", "klda2", "sklda", "sklda2")]

index=2
name = paste("simu.0618.kernelconpound100", index, "sparse.csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst2 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst2) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst2[,c("klda")]

rst2test <- rst2[,c("klda","klda2", "sklda", "sklda2")]

xtable(t(cbind(rst3, rst2test,rst4test)),digit=3)

m2 <- t(cbind(rst3, rst2test,rst4test))[,2]
########
########
#########
######

index=3
name = paste("simu.0618.kernelconpound200", index, "sparse.csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst3 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst3) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst3

index=4
name = paste("simu.0618.kernelconpound200", index, "sparse.csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst4 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst4) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst4[,c("klda")]

rst4test <- rst4[,c("klda", "klda2", "sklda", "sklda2")]

index=2
name = paste("simu.0618.kernelconpound200", index, ".csv", sep="")
rst <- read.table(name, sep=",", header=TRUE)
head(rst[,1:10])

mean = apply(rst[,-1], 1, mean)
mean

rst2 = matrix(apply(rst[,-1],1,mean), nrow=2, byrow=T)
colnames(rst2) = c("klda","klda2","lda","sda","sda2","kda","kda2","sklda","sklda2","logit.lasso1", "logit.lasso2")

rst2[,c("klda")]

rst2test <- rst2[,c("klda","klda2", "sklda", "sklda2")]

xtable(t(cbind(rst3, rst2test,rst4test)),digit=3)
m3 <- t(cbind(rst3, rst2test,rst4test))[,2]

xtable(cbind(m1, m2, m3), digits=3)






\begin{table}[ht]
\centering
\begin{tabular}{rrrr}
  \hline
 & m1 & m2 & m3 \\ 
  \hline

  lda & 0.561 & 0.564 & 0.531 \\ 
  sda & 0.771 & 0.780 & 0.751 \\ 
  sda2 & 0.602 & 0.584 & 0.580 \\ 
  kda & 0.799 & 0.808 & 0.965 \\ 
  kda2 & 0.569 & 0.575 & 0.529 \\ 

  logit.lasso1 & 0.964 & 0.963 & 0.968 \\ 
  logit.lasso2 & 0.944 & 0.948 & 0.946 \\ 

klda & 0.959 & 0.562 & 0.980 \\ 
  klda2 & 0.959 & 0.811 & 0.980 \\ 
   sklda & 0.984 & 0.562 & 0.982 \\ 
  sklda2 & 0.984 & 0.830 & 0.982 \\ 
  klda.1 & 0.562 & 0.544 & 0.529 \\ 
  klda2.1 & 0.811 & 0.836 & 0.849 \\ 
  sklda.1 & 0.562 & 0.544 & 0.529 \\ 
  sklda2.1 & 0.830 & 0.836 & 0.849 \\ 
  klda.2 & 0.864 & 0.562 & 0.776 \\ 
  klda2.2 & 0.776 & 0.811 & 0.748 \\ 
  sklda.2 & 0.877 & 0.562 & 0.777 \\ 
  sklda2.2 & 0.794 & 0.830 & 0.747 \\ 
   \hline
\end{tabular}
\end{table}


