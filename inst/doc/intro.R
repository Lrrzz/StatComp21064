## ----eval=FALSE---------------------------------------------------------------
#  gibbsR <- function(N, thin, n, aa, bb) {
#    mat <- matrix(nrow = N, ncol = 2)
#    x <- y <- 0
#    for (i in 1:N) {
#      for (j in 1:thin) {
#        x <- rbinom(1,n,y)
#        y <- rbeta(1,x+aa,n-x+bb)
#      }
#      mat[i, ] <- c(x, y)
#    }
#    mat
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix gibbsC(int N, int thin, int n, int aa, int bb) {
#    NumericMatrix mat(N, 2);
#    double x = 0, y = 0;
#    for(int i = 0; i < N; i++) {
#      for(int j = 0; j < thin; j++) {
#        x = rbinom(1,n,y)[0];
#        y = rbeta(1,x+aa,n-x+bb)[0];
#      }
#      mat(i, 0) = x;
#      mat(i, 1) = y;
#    }
#    return(mat);
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp21064)
library(microbenchmark)
tm2 <- microbenchmark(
  vR = gibbsR(1e2, 10, 22, 2, 2),
  vC = gibbsC(1e2, 10, 22, 2, 2)
)
knitr::kable(summary(tm2)[,c(1,3,5,6)])

## ----eval=FALSE---------------------------------------------------------------
#  gradient.descent<-function(x,y,epsilon=1e-14,maxit=1000,stepsize=1e-3,alpha=0.25,b=0.8,stepmethod="backtracking",verbose=TRUE,plotLoss=TRUE){
#  
#    stopifnot(stepsize<=1e-3)
#    m<-nrow(x)
#    x<-cbind(rep(1,m),x)
#    n<-ncol(x)
#    beta<-numeric(n)
#    iter<-1
#    error<-1
#    loss<-crossprod(x%*%beta-y)
#  
#    while((error>epsilon)&&(iter<maxit)){
#      h<-x %*% beta
#      grad<-crossprod(x,(h-y))
#      beta.new<-beta-stepsize*grad
#      loss.new<-crossprod(x%*%beta.new-y)
#  
#      #update for next iter
#      if(stepmethod=="backtracking")
#        if(loss.new>(loss[iter]+alpha*stepsize*sum(grad*grad)))
#          stepsize=stepsize*b
#  
#      error<-max(abs(beta.new-beta))
#      loss<-append(loss,loss.new)
#      beta<-beta.new
#      iter<-iter+1 #thistime
#  
#      if(verbose&& iter%%10==0)
#        message(paste('Iteration:',iter))
#    }
#  
#    if(plotLoss)
#      plot(loss, type="l",bty="n")
#  
#    return(list(par=beta,RSE=sqrt(crossprod(x%*%beta.new-y)/(nrow(x)- ncol(x))),iter=iter))
#  }

## -----------------------------------------------------------------------------
n=100
x1<-rnorm(n)
x2<-rnorm(n)
y=1+0.5*x1+x2+rnorm(n,sd=0.1)
x<-cbind(x1,x2)
gradient.descent(x,y)

## ----warning=FALSE------------------------------------------------------------
library(mvtnorm)

## ----eval=FALSE---------------------------------------------------------------
#  EM<-function(x,C,iter,mu0){
#    n <- nrow(x)
#    p <- ncol(x)
#    # initialize Z, mu, Sig, pi
#    Z <- matrix(c(rep(1, C), rep(0, (n - 1) * C)), n, C)
#    Sig <- matrix(rep(c(1, 0, 0, 1), C), ncol = p, byrow = T)
#    pi <- rep(1 / C, C)
#    mu<-mu0
#    # start EM here
#    for (t in 1:iter) {
#      for (k in 1:n) {
#        for (l in 1:C) {
#          Z[k,l] <- pi[l] *mvtnorm::dmvnorm(x[k,], mu[l,], Sig[(p * (l - 1) + 1):(p * l),])
#        }
#        Z[k,] <- Z[k,] / sum(Z[k,])
#      }
#      # update Z (E-step)
#  
#      pi <- colMeans(Z)
#      # update pi (M-step-1)
#  
#      for (i in 1:C) {
#        mu[i,] <- t(Z[,i]) %*% x / sum(Z[,i])
#        # update mu (M-step-2)
#  
#        sumsig <- Z[1,i] * (x[1,] - mu[i,]) %*% t(x[1,] - mu[i,])
#        for (k in 2:n) {
#          sumsig <- sumsig + Z[k,i] * (x[k,] - mu[i,]) %*% t(x[k,] - mu[i,])
#        }
#        Sig[(p * (i - 1) + 1):(p * i),] <- sumsig / sum(Z[,i])
#        # update Sigma (M-step-3)
#      }
#    }
#    return(list(mu=mu,Sig=Sig))
#  }

## -----------------------------------------------------------------------------
library(datasets)
x<-as.matrix(datasets::faithful)
EM(x,2,20,matrix(c(2,50,4,80),2,2,byrow=T))

