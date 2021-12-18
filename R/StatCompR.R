#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions (\code{gibbsR} and \code{vaccR}) and Cpp functions (\code{gibbsC} and \code{vaccC}).
#' @examples
#' \dontrun{
#' tm1 <- microbenchmark::microbenchmark(
#'   rnR = gibbsR(100,10,22,2,2),
#'   rnC = gibbsC(100,10,22,2,2)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' }
#' @import microbenchmark datasets
#' @importFrom Rcpp evalCpp
#' @importFrom stats rbeta rbinom
#' @useDynLib StatComp21064
NULL

#' @title A Gibbs sampler using R
#' @description A Gibbs sampler using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @param n the self-selected fixed parameter 'n' in ex9.8
#' @param aa the self-selected fixed parameter 'a' in ex9.8
#' @param bb the self-selected fixed parameter 'b' in ex9.8
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' rnR <- gibbsR(100,10,22,2,2)
#' par(mfrow=c(2,1));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
gibbsR <- function(N, thin, n, aa, bb) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1,n,y)
      y <- rbeta(1,x+aa,n-x+bb)
    }
    mat[i, ] <- c(x, y)
  }
  mat
}

#' @title Gradient descent method to find the coefficients of linear regression equation
#' @description Gradient descent method to find the coefficients of linear regression equation
#' @param x a m*n data matrix 
#' @param y observations size of m*1
#' @param epsilon termination condition, difference between the coefs of two consecutive iters
#' @param maxit the maximum number of iterations, default is 1000
#' @param stepsize just stepsize, must smaller than 1e-3, default is 1e-3
#' @param alpha parameter in the backtracking method, should in (0,0.5), default is 0.25
#' @param b parameter in the backtracking method, should in (0,1), default is 0.8
#' @param stepmethod optional "backtracking" or "fixed", default is "backtracking"
#' @param verbose whether to print out iterations, default is "TRUE"
#' @param plotLoss whether to plot loss, default is "TRUE"
#' @return a list contains coefficients,RSE and number of iterations
#' @examples
#' \dontrun{
#' n=100
#' x1<-rnorm(n)
#' x2<-rnorm(n)
#' y=1+0.5*x1+x2+rnorm(n,sd=0.1)
#' x<-cbind(x1,x2)
#' gradient.descent(x,y)
#' }
#' @export
gradient.descent<-function(x,y,epsilon=1e-14,maxit=1000,stepsize=1e-3,alpha=0.25,b=0.8,stepmethod="backtracking",verbose=TRUE,plotLoss=TRUE){
  
  stopifnot(stepsize<=1e-3)
  m<-nrow(x)
  x<-cbind(rep(1,m),x)
  n<-ncol(x)
  beta<-numeric(n)
  iter<-1
  error<-1
  loss<-crossprod(x%*%beta-y)
  
  while((error>epsilon)&&(iter<maxit)){
    h<-x %*% beta
    grad<-crossprod(x,(h-y))
    beta.new<-beta-stepsize*grad
    loss.new<-crossprod(x%*%beta.new-y)
    
    #update for next iter
    if(stepmethod=="backtracking")
      if(loss.new>(loss[iter]+alpha*stepsize*sum(grad*grad)))
        stepsize=stepsize*b
    
    error<-max(abs(beta.new-beta))
    loss<-append(loss,loss.new)
    beta<-beta.new
    iter<-iter+1 #thistime
    
    if(verbose&& iter%%10==0)
      message(paste('Iteration:',iter))
  }
  
  if(plotLoss)
    plot(loss, type="l",bty="n")
  
  return(list(par=beta,RSE=sqrt(crossprod(x%*%beta.new-y)/(nrow(x)- ncol(x))),iter=iter))
}

#' @title EM algorithm estimates the joint density function of the normal mixed model
#' @description EM algorithm estimates the joint density function of the normal mixed model
#' @param x data
#' @param C a count, guess there’s a mixture of C Gaussian distributions
#' @param iter count for iteration times
#' @param mu0 initialize mu, a matrix size of C*ncol(x), you must first plot x, and estimate mu0 based on the picture, otherwise it will go wrong
#' @return a list contains mu and Sigma
#' @import mvtnorm
#' @examples
#' \dontrun{
#' x<-as.matrix(datasets::faithful)
#' plot(x)
#' EM(x,2,20,matrix(c(2,50,4,80),2,2,byrow=T))
#' }
#' @export
EM<-function(x,C,iter,mu0){
  n <- nrow(x)
  p <- ncol(x)
  # initialize Z, mu, Sig, pi
  Z <- matrix(c(rep(1, C), rep(0, (n - 1) * C)), n, C)
  Sig <- matrix(rep(c(1, 0, 0, 1), C), ncol = p, byrow = T)
  pi <- rep(1 / C, C)
  mu<-mu0
  # start EM here
  for (t in 1:iter) {
    for (k in 1:n) {
      for (l in 1:C) {
        Z[k,l] <- pi[l] *mvtnorm::dmvnorm(x[k,], mu[l,], Sig[(p * (l - 1) + 1):(p * l),])
      }
      Z[k,] <- Z[k,] / sum(Z[k,])
    }
    # update Z (E-step)
    
    pi <- colMeans(Z)
    # update pi (M-step-1) 
    
    for (i in 1:C) {
      mu[i,] <- t(Z[,i]) %*% x / sum(Z[,i])
      # update mu (M-step-2)
      
      sumsig <- Z[1,i] * (x[1,] - mu[i,]) %*% t(x[1,] - mu[i,]) 
      for (k in 2:n) {
        sumsig <- sumsig + Z[k,i] * (x[k,] - mu[i,]) %*% t(x[k,] - mu[i,])
      }
      Sig[(p * (i - 1) + 1):(p * i),] <- sumsig / sum(Z[,i])
      # update Sigma (M-step-3) 
    }
  }
  return(list(mu=mu,Sig=Sig))
}