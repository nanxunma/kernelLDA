sklda.rayleigh2 <- function(x,y,H, Omega=diag(dim(x)[2]), lambda=1){
x=scale(x,TRUE,FALSE)
H = scale(H,TRUE,FALSE)
Sbk0 = t(x)%*%H%*%y%*%solve(t(y)%*%y)%*%t(y)%*%H%*%x
Swk0 = t(x)%*%x - Sbk0 + lambda*Omega
b.0class = eigen(solve(Swk0)%*%Sbk0)
#image(t((Swk0)%*%Sbk0))
#plot((Re(eigen(Swk0%*%Sbk0)$vector[,1])))
#image(t((Swk0)))
#image(t(Sbk0))
ncol = dim(y)[2]-1
return(Re(b.0class$vector[,1:ncol]))
}

sklda.os2 <- function(x, y, H, Omega=diag(dim(x)[2]), lambda=1, stop=-p, maxIte=100, Q=K-1, tol=1e-6, ...){
    class.ind <- function(cl) {
      Ik=diag(length(levels(cl)))
      x=Ik[as.numeric(cl),]
      dimnames(x) <- list(names(cl), levels(cl))
      x
    }
    orth.Q=function(Dpi,Qj,theta){
  Ik = diag(dim(Qj)[1])
  theta1 = (Ik - Qj%*%t(Qj)%*%Dpi) %*% theta
  rst = theta1/drop(sqrt(t(theta1)%*%Dpi%*%theta1))
return(rst)
    }

    rtheta=function(K,Dp){
      jj=rnorm(K);
      jj/(drop(sqrt(t(jj)%*%Dp%*%jj)))
    }


    if(is.factor(y))
      {
        classes <- levels(y)
        factorY <- y
        y <- class.ind(y)
      } else {
        classes <- colnames(y)
        factorY <- factor(colnames(y)[apply(y, 1, which.max)])
      }
    if(!is.matrix(x)) x <- as.matrix(x)
    x=scale(x,TRUE,FALSE)##This centering is essential for the trivial solution to disappear
    H = scale(H,TRUE,FALSE)
    predNames <- colnames(x)
    
    n <- dim(x)[1]
    p <- dim(x)[2]
    K <- dim(y)[2]
    ones=rep(1,K)
    if(Q>(K-1))stop("at most K-1 variates allowed")
    

    ##temp for Y
    ss = svd(H)
    U = diag(sqrt(ss$d))%*%t(ss$v)

    Dpi = 1/n*t(y)%*%y

    Yhat <- matrix(0,n,Q)
    b <- matrix(0,p,Q)
    rss <- rep(0,Q)

    theta=matrix(0,K,Q)
    Qj=matrix(ones,K,1)





    for(j in 1:Q){
      RSS <- 1e6
      RSSold <- Inf
      ite <- 0
      thetaj1=rtheta(K,Dpi)
      thetaj=orth.Q(Dpi,Qj,thetaj1)
      while (abs(RSSold-RSS)/RSS > tol & ite < maxIte){ 
        RSSold <- RSS
        ite <- ite + 1
        ## 1. Estimate beta:    
        Yc <- y%*%thetaj 
        #beta<- solvebeta(x, Yc, paras=c(lambda, abs(stop[j])),sparse=sparse) # elasticnet to estimate beta
        #beta <- solve(crossprod(x)) %*% crossprod(x, Yc)
        beta <- solve(t(x)%*%x + lambda*Omega) %*% t(x)%*%H%*%Yc
        yhatj=x%*%beta
        thetaj0 = solve(Dpi)%*%t(y)%*%H%*%x%*%beta
        thetaj=orth.Q(Dpi,Qj,thetaj0)
        RSS=sum((yhatj-Yc)^2)+lambda*sum(beta^2)
      rss[j]=RSS
    }
      Qj=cbind(Qj,thetaj)
      theta[,j]=thetaj

      
      #Theta[,j]=thetaj
      b[,j]=beta
    }

                                        # remove predictors which are not included (do not have non-zero parameter estimates)
      notZero <- apply(b, 1, function(x) any(x != 0))
      b <- as.matrix(b[notZero,])
      origP <- ncol(x)
      x <- x[, notZero, drop = FALSE]
      varNames <- colnames(x)

### remove directions with only zero elements (this can be caused by a too high weight on L1-penalty)
      if (is.vector(b)){b<-t(as.matrix(b))}
      notZeroC <- apply(b,2,function(x) any(x!=0))
      b <- as.matrix(b[,notZeroC])
      
      
      sl <- x %*% b
      sl2 <- H%*%x%*%b
      sl3 <- U%*%x%*%b
      colnames(sl) <- paste("score", 1:ncol(sl), sep = "")
      lobj<-lda(sl, factorY, ...)
      lobj1<-lda(sl, factorY)
      lobj2<-lda(sl2, factorY)
      lobj3<-lda(sl3, factorY)
    
    

    structure(
              list(call = match.call(),
                   beta = b,
                   theta = theta,
                   varNames = varNames,
                   varIndex = which(notZero),
                   origP = origP,
                   rss = rss,
                   fit = lobj2,
                   classes = classes,
                   lambda = lambda,
                   stop = stop),
              class = "klda")
  }



match.beta <- function(label, input){
  b = as.numeric(input)-1
  l1 = max(mean(b==label), mean(1-b == label))
  return(l1)
}

nmat <- function(M){
  d = dim(M); n = d[1]
  II = matrix(1,n,n)
  Mgc = (diag(n)-(1/n)*II)%*%M%*%(diag(n)-(1/n)*II)
}

pred.os <- function(obj, x, K, labels){
  x = scale(x, TRUE, FALSE)
  H = scale(K, TRUE, FALSE)
  ss = svd(H)
  U = diag(sqrt(ss$d))%*%t(ss$v)
  beta = obj$beta
  T = U %*% x %*% beta
  fit.obj = lda(T, labels)
  newX = x %*% beta
  pred = predict(fit.obj, newX)
  return(pred)
}

pred.os2 <- function(obj, x, K, labels){
  x = scale(x, TRUE, FALSE)
  H = scale(K, TRUE, FALSE)
  ss = svd(H)
  U = diag(sqrt(ss$d))%*%t(ss$v)
  beta = obj$beta
  T = U %*% x %*% beta
  fit.obj = lda(T, labels)
  newX = x %*% beta
  pred = predict(fit.obj, T)
  return(pred)
}

pred.os3 <- function(obj, x, K, labels){
  x = scale(x, TRUE, FALSE)
  H = scale(K, TRUE, FALSE)
  ss = svd(H)
  U = diag(sqrt(ss$d))%*%t(ss$v)
  beta = obj$beta
  T = H %*% x %*% beta
  fit.obj = lda(T, labels)
  newX = x %*% beta
  pred = predict(fit.obj, T)
  return(pred)
}

pred.os4 <- function(obj, x, K, labels){
  x = scale(x, TRUE, FALSE)
  H = scale(K, TRUE, FALSE)
  ss = svd(H)
  U = diag((ss$d))%*%t(ss$v)
  beta = obj$beta
  T = U %*% x %*% beta
  fit.obj = lda(T, labels)
  pred = predict(fit.obj, T)
  return(pred)
}


match.3class = function(label, input){
  L0 = label==0
  I0 = input==0
  L1 = label==1
  I1 = input==1
  L2 = label==2
  I2 = input==2
  match0 = max(match.beta(L0,I0),match.beta(L0,I1), match.beta(L0,I2))
  match1 = max(match.beta(L1,I0),match.beta(L1,I1), match.beta(L1,I2))
  match2 = max(match.beta(L2,I0),match.beta(L2,I1), match.beta(L2,I2))
  return(c(match0, match1, match2))
}

















