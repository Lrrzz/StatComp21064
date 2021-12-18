## -----------------------------------------------------------------------------
set.seed(923)
rRayleigh = function(n, sigma=1) {
                      y = rexp(n, 1)
                      x = sigma*sqrt(2*y)
                      return(x)
                    }
         n = 1000
         X1 = rRayleigh(n,1)
         X3 = rRayleigh(n,3)
         X5 = rRayleigh(n,5)
         #compare mode via histogram
         hist(X1, main='sigma=1')
         hist(X3, main='sigma=3')
         hist(X5, main='sigma=5')

## -----------------------------------------------------------------------------
set.seed(923)
rRayleigh2 = function(n, sigma=1) {
                      u<-runif(n)
                      x<-sqrt(log(1/u)*2*sigma^2)
                      return(x)
                    }
         n = 1000
         X1 = rRayleigh2(n,1)
         X3 = rRayleigh2(n,3)
         X5 = rRayleigh2(n,5)
         #compare mode via histogram
         hist(X1, main='sigma=1')
         hist(X3, main='sigma=3')
         hist(X5, main='sigma=5')

## -----------------------------------------------------------------------------
set.seed(923)
n = 1000
p1 = 0.5
x1 = rnorm(n,0,1)
x2 = rnorm(n,3,1)
u = runif(n)
k = as.integer( u < p1 )
x = k*x1 + (1-k)*x2
hist( x, prob=T, ylab='', main='p1=0.5', xlab='', ylim=c(0,.3) )
lines( density(x) )

## -----------------------------------------------------------------------------
for (p1 in seq(.1,.9,.1) ){
  x1 = rnorm(n,0,1)
  x2 = rnorm(n,3,1)
  u = runif(n)
  k = as.integer( u < p1 )
  x = k*x1 + (1-k)*x2
  hist( x, prob=T, ylab='', xlab='', xlim=c(-6,6),ylim=c(0,.4),main=paste('p =',p1) )
  lines( density(x) )
}

## -----------------------------------------------------------------------------
rP_G = function(n,lambda,t0,r,beta) {
                      pp<-numeric(n)
                      ppgg<-numeric(n)
                      for(i in 1:n){
                         Tn<-rexp(100,lambda)
                         Sn<-cumsum(Tn)
                         nt0<-min(which(Sn>t0))-1
                         pp[i]<-nt0
                      }
                      for(j in 1:n){
                        ppgg[j]=sum(rgamma(pp[j],r,beta))
                      }
                      return(ppgg)
                    }


## -----------------------------------------------------------------------------
set.seed(923)
x<-rP_G(n=2000,lambda=2,t0=10,r=1,beta=2)
mean(x)
var(x)

## -----------------------------------------------------------------------------
set.seed(929)
x<-rP_G(n=2000,lambda=3,t0=10,r=2,beta=3)
mean(x)
var(x)

## -----------------------------------------------------------------------------
cdfBeta = function(x, alpha=3, beta=3, m=10000 ) {
if ( any(x < 0) ) return (0)
stopifnot( x < 1 )
t = runif( m, min=0, max=x )
h = (x-0) * (1/beta(alpha,beta)) * t^(alpha-1) * (1-t)^(beta-1)
cdf = mean( h ) 
return( min(1,cdf) )
}
set.seed(930) 
for (i in 1:9) {
         print( c(i/10 ,cdfBeta(i/10,3,3,10000), pbeta(i/10,3,3)) )
                      }

## -----------------------------------------------------------------------------
set.seed(222)
rRayleigh = function(n, sigma, antithetic=T){
                      u<-runif(n)
                      if(!antithetic) v<-runif(n)
                      else
                        v<-1-u
                      x<-(sqrt(log(1/(1-u))*2*sigma^2) +sqrt(log(1/(1-v))*2*sigma^2))/2
                      x
}

sigma <- c(1,3,5,7,9)
var_redu <- numeric(length(sigma))
for(i in 1:length(sigma)){
  MC1 <- rRayleigh(sigma=sigma[i],n=2000,antithetic = FALSE)
  MC2 <- rRayleigh(sigma=sigma[i],n=2000,antithetic = TRUE)
  var_redu[i] <- (var(MC1)-var(MC2))/var(MC1)
}

result <- rbind(sigma,var_redu)
knitr::kable(round(result,3))

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- var <- numeric(2)

g <- function(x){
  x^2*exp(-x^2/2)*(x>1)/(sqrt(2*pi))
}

# f1
u<-runif(m)
x<-1/(1-u)
g_f<-g(x)*x^2
theta.hat[1]<-mean(g_f)
var[1]<-var(g_f)

# f2
x <- rnorm(m,mean = 1.5)
g_f <- g(x)/dnorm(x,mean = 1.5)
theta.hat[2] <- mean(g_f)
var[2] <- var(g_f)

theta.hat
var

## -----------------------------------------------------------------------------
y <- seq(1,10,.01)
plot(y,g(y)*y^2,col=1,lty=1,ylim=c(-0.5,3),"l",ylab="",xlab="")
lines(y,g(y)/dnorm(y,mean = 1.5),col=2,lty=1)
legend("topright",legend=c(1:2),lty=1:2,col=1:2)

## -----------------------------------------------------------------------------
set.seed(1)
m <- 10000
theta.hat <- var <- numeric(1)
# f2
x <- rnorm(m,mean = 1.5)
g_f <- g(x)/dnorm(x,mean = 1.5)
theta.hat<- mean(g_f)

true <- integrate(g,1,upper = Inf)$value
c(theta.hat,true)

## -----------------------------------------------------------------------------
set.seed(222)
alpha = 0.05
n = 20
m = 10000

t_UCL = numeric(m)
t_LCL = numeric(m)
UCL<-numeric(m)

for(i in 1:m)
{
    x = rchisq(n, 2)
    t_LCL[i] = mean(x) - qt(1-alpha / 2, df=n-1, lower.tail = TRUE)*sd(x)/sqrt(n)
    t_UCL[i] = mean(x) - qt(alpha / 2, df=n-1, lower.tail = TRUE)*sd(x)/sqrt(n)
    UCL[i]<-(n-1)*var(x)/qchisq(alpha,df=n-1)
}

mean(t_LCL < 2 & t_UCL > 2)
mean(UCL>4)

## -----------------------------------------------------------------------------
set.seed(222)
alpha = 0.05
n = 100
m = 10000

p_chi<-p_uni<-p_exp<-numeric(m)
for(i in 1:m){
  x_chi<-rchisq(n,df=1)
  x_uni<-runif(n,0,2)
  x_exp<-rexp(n,1)
  
  p_chi[i]=abs((mean(x_chi)-1)/sd(x_chi)*sqrt(n) )>=qt(0.975,n-1,lower.tail = TRUE)
  p_uni[i]=abs((mean(x_uni)-1)/sd(x_uni)*sqrt(n) )>=qt(0.975,n-1,lower.tail = TRUE)
  p_exp[i]=abs((mean(x_exp)-1)/sd(x_exp)*sqrt(n) )>=qt(0.975,n-1,lower.tail = TRUE)
}
mean(p_chi)
mean(p_uni)
mean(p_exp)

## -----------------------------------------------------------------------------
res <- prop.test(x = c(6510, 6760), n = c(10000, 10000))
res 

## -----------------------------------------------------------------------------
Mardia.sk<-function(x){
  n <- nrow(x)
  xbar <- apply(x, 2, mean)
  Sigma_hat <- (n-1)/n*cov(x)
  ni<-solve(Sigma_hat)
  sk<-numeric(1)
  for(i in 1:n)
    for(j in 1:n)
      sk <- sk + ((x[i,]-xbar)%*%ni%*%(x[j,]-xbar))^3
  return(sk/n^2)
}

## ----warning=F,eval=FALSE-----------------------------------------------------
#  library(MASS)
#  library(knitr)
#  n<-c(10,20,30,50,100,500)
#  d <- 2
#  cv <- 6/n*qchisq(0.05, d*(d+1)*(d+2)/6, lower.tail = F)
#  p.reject<-numeric(length(n))
#  m<-10000
#  set.seed(123)
#  
#  for(i in 1:length(n)){
#    Mardia.sktests<-numeric(m)
#    for(j in 1:m){
#      x <- mvrnorm(n[i], mu = rep(0,d), Sigma = diag(1,d))
#      Mardia.sktests[j] <- as.integer(Mardia.sk(x)>= cv[i])
#    }
#    p.reject[i]<-mean(Mardia.sktests)
#  }
#  result <- rbind(n,p.reject)
#  rownames(result) <- c("n","estimate")
#  kable(result)

## -----------------------------------------------------------------------------
library(MASS)
library(knitr)
alpha<-0.1
n<-30
m<-2500
d<-2
set.seed(123)
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
cv <- 6/n*qchisq(alpha, d*(d+1)*(d+2)/6, lower.tail = F)
for (j in 1:N) {
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) {
    index<- sample(c(1, 0), replace = TRUE,size = n, prob = c(1-e, e))
    x <- index*mvrnorm(n, rep(0,d), diag(1,d))+(1-index)*mvrnorm(n, rep(0,d), diag(100,d))
    sktests[i] <- as.integer(Mardia.sk(x) >= cv)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m)
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap)
library(knitr)

## ----warning=FALSE------------------------------------------------------------
lambda_hat=eigen(cov(scor))$values
theta_hat=lambda_hat[1]/sum(lambda_hat)
theta_hat

statis<-function(x,i){
  lambda.hat=eigen(cov(x[i,]))$values
  theta.hat=lambda.hat[1]/sum(lambda.hat)
  theta.hat
}

library(boot)
set.seed(222)
B<-2000
results<-boot(data=scor,statistic=statis,R=B)
results

## -----------------------------------------------------------------------------
thetai = function(x){
  lambda<-eigen(cov(x))$values
  lambda[1]/sum(lambda)
}
n =nrow(scor)
theta.jack = numeric(n)
x = as.matrix(scor)

for (i in 1:n) {
  theta.jack[i] = thetai(x[-i,])
}
theta.jack.bar = mean(theta.jack)
theta.jack.bar
bias.jack = (n-1)*( theta.jack.bar - theta_hat )
se.jack = sqrt( (n-1)*mean( (theta.jack-theta.jack.bar)^2 ) )
result <- cbind(theta_hat,theta.jack.bar,bias.jack,se.jack)
colnames(result) <- c("original","theta.jack.bar","bias","std.error")
kable(result)

## -----------------------------------------------------------------------------
boot.ci(results, type=c('perc','bca'))

## -----------------------------------------------------------------------------
sample.skewness = function (x,i) {
    mu = mean(x[i])
    n = length(x[i])
    a = 1/n * sum((x[i] - mu)^3)
    b = (1/n * sum((x[i] - mu)^2))^1.5
    return (a/b)
}

## -----------------------------------------------------------------------------
theta=0
n<-100
m<-1e4
set.seed(222)
ci.norm<-ci.basic<-ci.perc<-matrix(0,m,2)
for(i in 1:m){
  xx<-rnorm(n,2,3)
  de<-boot(data=xx,statistic=sample.skewness,R=999)
  ci<-boot.ci(de,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}
cat("norm:","coverage=",mean(ci.norm[,1]<=theta&ci.norm[,2]>=theta),"left.miss=",mean(ci.norm[,1]>theta),"right.miss",mean(ci.norm[,2]<theta),"\n",
    "basic:","coverage=",mean(ci.basic[,1]<=theta&ci.basic[,2]>=theta),"left.miss=",mean(ci.basic[,1]>theta),"right.miss",mean(ci.basic[,2]<theta),"\n",
    "perc:","coverage=",mean(ci.perc[,1]<=theta&ci.perc[,2]>=theta),"left.miss=",mean(ci.perc[,1]>theta),"right.miss",mean(ci.perc[,2]<theta)
    )

## -----------------------------------------------------------------------------
df=5
theta=sqrt(8/df)
n<-100
m<-1e4
set.seed(222)
ci.norm<-ci.basic<-ci.perc<-matrix(0,m,2)
for(i in 1:m){
  xx<-rchisq(n, df = df)
  de<-boot(data=xx,statistic=sample.skewness,R=999)
  ci<-boot.ci(de,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}
cat("norm:","coverage=",mean(ci.norm[,1]<=theta&ci.norm[,2]>=theta),"left.miss=",mean(ci.norm[,1]>theta),"right.miss",mean(ci.norm[,2]<theta),"\n",
    "basic:","coverage=",mean(ci.basic[,1]<=theta&ci.basic[,2]>=theta),"left.miss=",mean(ci.basic[,1]>theta),"right.miss",mean(ci.basic[,2]<theta),"\n",
    "perc:","coverage=",mean(ci.perc[,1]<=theta&ci.perc[,2]>=theta),"left.miss=",mean(ci.perc[,1]>theta),"right.miss",mean(ci.perc[,2]<theta)
)

## -----------------------------------------------------------------------------
library(boot)
dCov<-function(x,y){
  x<-as.matrix(x)
  y<-as.matrix(y)
  n<-nrow(x)
  m<-nrow(y)
  if(n!=m||n<2)
    strp("Sample sizes must agree")
  if (!(all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  Akl<-function(x){
    d<-as.matrix(dist(x))
    m1<-rowMeans(d)
    m2<-colMeans(d)
    M<-mean(d)
    a<-sweep(d,1,m1)
    b<-sweep(a,2,m2)
    b+M
  }
  A<-Akl(x)
  B<-Akl(y)
  sqrt(mean(A*B))
}
ndCov2<-function(z,ix,dims){
    p<-dims[1]
    q<-dims[2]
    d<-p+q
    x<-z[,1:p]
    y<-z[ix,-(1:p)]
    return(nrow(z)*dCov(x,y)^2)
}

## -----------------------------------------------------------------------------
data("iris")
z<-as.matrix(iris[1:50,1:4])
set.seed(12345)
boot.obj<-boot(data=z,statistic=ndCov2,R=9999,sim="permutation",dims=c(2,2))
tb<-c(boot.obj$t0,boot.obj$t)
p.cor<-mean(tb>=tb[1])
hist(tb,nclass="scott",main="",freq=F,xlab="Replicates correlation",sub=list(substitute(paste(hat(p), " = ",pvalue),list(pvalue = p.cor)),col="red"))
abline(v=boot.obj$t0,col="red",lwd=2)

## -----------------------------------------------------------------------------
ndCov22<-function(z,i,dims){
    p<-dims[1]
    q<-dims[2]
    d<-p+q
    x<-z[,1:p]
    y<-z[i,-(1:p)]
    return(cor.test(x,y,method="spearman",exact=F)$estimate)
}
  set.seed(12345)
  boot.obj<-boot(data=z,statistic=ndCov22,R=9999,sim="permutation",dims=c(2,2))
  tb<-c(boot.obj$t0,boot.obj$t)
  p.cor<-mean(tb>=tb[1])
  hist(tb,nclass="scott",main="",freq=F,xlab="Replicates correlation",sub=list(substitute(paste(hat(p), " = ",pvalue),list(pvalue = p.cor)),col="red"))
  abline(v=boot.obj$t0,col="red",lwd=2)

## ----warning=F,eval=FALSE-----------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  
#  Tn<-function(z,ix,sizes,k){
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1)
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1); i2 <- sum(block2 > n1)
#    (i1 + i2) / (k * n)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  #(a)
#  set.seed(123)
#  m<-200
#  nx<-30
#  p.values<-matrix(NA,m,3)
#  colnames(p.values)<-c("NN","energy","ball")
#  for(i in 1:m){
#  x<-rnorm(nx,0,1)
#  y<-rnorm(nx,0,2)
#  z<-c(x,y)
#  N<-c(nx,nx)
#  boot.obj <- boot(data=z,statistic=Tn,R=999,sim = "permutation", sizes = N,k=3)
#  ts <- c(boot.obj$t0,boot.obj$t)
#  p.values[i,1] <- mean(ts>=ts[1])
#  p.values[i,2]<-eqdist.etest(z,sizes=N,R=999)$p.value
#  p.values[i,3]<-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  power<-numeric(3)
#  power<-colMeans(p.values<0.1)
#  power

## ----eval=FALSE---------------------------------------------------------------
#  #(b)
#  set.seed(123)
#  m<-200
#  nx<-30
#  p.values<-matrix(NA,m,3)
#  colnames(p.values)<-c("NN","energy","ball")
#  for(i in 1:m){
#  x<-rnorm(nx,0,1)
#  y<-rnorm(nx,1,2)
#  z<-c(x,y)
#  N<-c(nx,nx)
#  boot.obj <- boot(data=z,statistic=Tn,R=999,sim = "permutation", sizes = N,k=3)
#  ts <- c(boot.obj$t0,boot.obj$t)
#  p.values[i,1] <- mean(ts>=ts[1])
#  p.values[i,2]<-eqdist.etest(z,sizes=N,R=999)$p.value
#  p.values[i,3]<-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  power<-numeric(3)
#  power<-colMeans(p.values<0.1)
#  power

## ----eval=FALSE---------------------------------------------------------------
#  #(c)
#  set.seed(123)
#  m<-200
#  nx<-30
#  p.values<-matrix(NA,m,3)
#  colnames(p.values)<-c("NN","energy","ball")
#  for(i in 1:m){
#  x<-rt(nx,df=1)
#  index<-sample(c(1,0),size=nx,replace=T,prob=c(0.5,0.5))
#  y<-index*rnorm(nx,0,1)+(1-index)*rnorm(nx,1,2)
#  z<-c(x,y)
#  N<-c(nx,nx)
#  boot.obj <- boot(data=z,statistic=Tn,R=999,sim = "permutation", sizes = N,k=3)
#  ts <- c(boot.obj$t0,boot.obj$t)
#  p.values[i,1] <- mean(ts>=ts[1])
#  p.values[i,2]<-eqdist.etest(z,sizes=N,R=999)$p.value
#  p.values[i,3]<-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  power<-numeric(3)
#  power<-colMeans(p.values<0.1)
#  power

## ----eval=FALSE---------------------------------------------------------------
#  #(d)
#  set.seed(123)
#  m<-200
#  p.values<-matrix(NA,m,3)
#  colnames(p.values)<-c("NN","energy","ball")
#  for(i in 1:m){
#  x<-rnorm(5,0,1)
#  y<-rnorm(50,0,1)
#  z<-c(x,y)
#  N<-c(length(x),length(y))
#  boot.obj <- boot(data=z,statistic=Tn,R=999,sim = "permutation", sizes = N,k=3)
#  ts <- c(boot.obj$t0,boot.obj$t)
#  p.values[i,1] <- mean(ts>=ts[1])
#  p.values[i,2]<-eqdist.etest(z,sizes=N,R=999)$p.value
#  p.values[i,3]<-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  power<-numeric(3)
#  power<-colMeans(p.values<0.1)
#  power

## -----------------------------------------------------------------------------
f<-function(x,eta=0,theta=1){
  stopifnot(theta>0)
  return(1/(pi*theta*(1+((x-eta)/theta)^2)))
}

## -----------------------------------------------------------------------------
cauchy.chain<-function(f=f,m,sigma,b,x0){
x<-numeric(m)
k<-0
x[1]<-x0
u<-runif(m)
for(i in 2:m){
  xt<-x[i-1]
  y<-rnorm(1,xt,sigma)
  num<-f(y)*dnorm(xt,y,sigma)
  den<-f(xt)*dnorm(y,xt,sigma)
  if(u[i]<=num/den) x[i]<-y
  else{
    x[i]<-xt
    k<-k+1
  }
}
index<-(b+1):m
y1<-x[index]
return(y1)
}

## -----------------------------------------------------------------------------
set.seed(222)
m<-10000
sigma<-1
mu0<-0
b<-1000
x00<-rnorm(1,mu0,sigma)
y1<-cauchy.chain(f=f,m=m,sigma=sigma,b=b,x0=x00)
index<-(b+1):m
plot(index,y1,type="l",main="",ylab="x")

q<-seq(0.1,0.9,0.1)
qxt<-quantile(y1,q)
qcc<-qcauchy(q)
result<-rbind(qxt,qcc)
knitr::kable(round(result,3))

a<-ppoints(100)
QC<-qcauchy(a)
Q<-quantile(y1,a)
hist(y1,breaks="scott",main="",xlab="",freq=F)
lines(QC,f(QC))

## -----------------------------------------------------------------------------
B.B.chain<-function(a,b,n,N,burn){
  X<-matrix(0,N,2)
  y<-(0.5*n + a)/(n+a+b)
  X[1,]<-c(floor( n*y ),y)
for(i in 2:N){
  x2<-X[i-1,2]
  X[i,1]<-rbinom(1,n,y)
  x1<-X[i,1]
  X[i,2]<-rbeta(1,x1+a,n-x1+b)
}
  b<-burn+1
  x<-X[b:N,]
  return(x)
}

## -----------------------------------------------------------------------------
set.seed(123)
N<-10000
burn<-1000
a<-2
b<-2
n<-22
x<-B.B.chain(a=a,b=b,n=n,N=N,burn=burn)
plot(x,xlab="x",ylab="y",main="Target sample chain",pch=16,cex=0.8)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)

psi.means <- rowMeans(psi)
B <- n * var(psi.means)
psi.w <- apply(psi, 1, "var")
W <- mean(psi.w)
v.hat <- W*(n-1)/n + (B/n)
r.hat <- v.hat / W
return(r.hat)
}

##ex9.3
##generate the chains
set.seed(123)
n<-15000
b<-1000
k<-4
X<-matrix(0,nrow=k,ncol=n)
for(i in 1:k)
  X[i,]<-cauchy.chain(f=f,m=n,sigma=sigma,b=0,x0=8)


#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

for (i in 1:k)
      if(i==1){
        plot((b+1):n,psi[i, (b+1):n], type="l",
            xlab='Index', ylab=bquote(psi))
      }else{
        lines(psi[i, (b+1):n], col=i)
    }
    
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R") 
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
##ex9.8
##generate the chains
set.seed(123)
b<-1000
k<-4 #number of chains
N<-15000
X<-matrix(0,nrow=k,ncol=N)
a<-2
b1<-2
n<-22
for(i in 1:k)
  X[i,]<-B.B.chain(a=a,b=b1,n=n,N=N,burn=0)[,2] #此处y全序列


#compute diagnostic statistics 
psi <- t(apply(X, 1, cumsum)) 
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi)) 
print(Gelman.Rubin(psi))

#plot psi for the four chains par(mfrow=c(2,2))
for (i in 1:k)
      if(i==1){
        plot((b+1):N,psi[i, (b+1):N], type="l",
            xlab='Index', ylab=bquote(psi))
      }else{
        lines(psi[i, (b+1):N], col=i)
    }
    par(mfrow=c(1,1)) #restore default
    
#plot the sequence of R-hat statistics 
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R") 
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
##ex9.8
##generate the chains
set.seed(123)
b<-1000
k<-4 #number of chains
N<-15000
X<-matrix(0,nrow=k,ncol=N)
a<-2
b1<-2
n<-22
for(i in 1:k)
  X[i,]<-B.B.chain(a=a,b=b1,n=n,N=N,burn=0)[,1]

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

for (i in 1:k)
      if(i==1){
        plot((b+1):N,psi[i, (b+1):N], type="l",
            xlab='Index', ylab=bquote(psi))
      }else{
        lines(psi[i, (b+1):N], col=i)
      }

#plot the sequence of R-hat statistics 
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R") 
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
compute_kth<-function(k,a,d=length(a)){
  a.norm2<-sum(abs(a)^2)
  y0<-a.norm2
  if(k==0) return(y0/2*exp(lgamma((d+1)/2)+lgamma(1.5)-lgamma(d/2+1)))
  y<-rep(0,k)
  y[1]<--1/2*a.norm2^2
  if(k>1) for(i in 2:k){
    y[i]<--y[i-1]/(2*i)*a.norm2
  }
  return(y[k]/(2*k+1)/(2*k+2)*exp(lgamma((d+1)/2)+lgamma(k+1.5)-lgamma(d/2+k+1)))
}

## -----------------------------------------------------------------------------
compute_sum<-function(k,a,d=length(a)){
  a.norm2<-sum(abs(a)^2)
  y<-rep(0,k)
  kth<-rep(0,k)
  y0<-a.norm2
  y[1]<--1/2*a.norm2^2
  kth[1]<-y[1]/3/4*exp(lgamma((d+1)/2)+lgamma(2.5)-lgamma(d/2+2))
  for(i in 2:k){
    y[i]<--y[i-1]/(2*i)*a.norm2
    kth[i]<-y[i]/(2*i+1)/(2*i+2)*exp(lgamma((d+1)/2)+lgamma(i+1.5)-lgamma(d/2+i+1))
  }
  y0/2*exp(lgamma((d+1)/2)+lgamma(1.5)-lgamma(d/2+1))+sum(kth)
}

## -----------------------------------------------------------------------------
a<-c(1,2)
compute_sum(k=100,a=a)
compute_sum(k=200,a=a)
compute_sum(k=300,a=a)

## -----------------------------------------------------------------------------
k = c( 4:25,100,500,1000 )
object = function( a, df ){
  a2 = a^2
  arg = sqrt( a2*df/(df + 1 - a2) )
  Sk = pt( q=arg, df=df, lower=F)
  arg = sqrt( a2*(df-1)/(df - a2) )
  Sk_1 = pt( q=arg, df=df-1, lower=F)
  return( Sk-Sk_1 )
}
for ( i in 1:length(k) )  {
    print( c( k[i], uniroot(object, lower=1,upper=2, df=k[i])$root ) )
                                 }

## -----------------------------------------------------------------------------
#11.5
solve_equation<-function(k){
   f=function(y,kk){
     (1+y^2/(kk-1))^(-kk/2)
     }
   object=function(a,kx){
     2*exp(lgamma(kx/2)-lgamma(kx/2-0.5))/sqrt(pi*(kx-1))*integrate(f,lower=0,upper=sqrt(a^2*(kx-1)/(kx-a^2)),kk=kx)$value-2*exp(lgamma(kx/2+0.5)-lgamma(kx/2))/sqrt(pi*kx)*integrate(f,lower=0,upper=sqrt(a^2*(kx)/(kx+1-a^2)),kk=kx+1)$value
     }
   uniroot(object, lower=1,upper=2,kx=k)$root
}

k = c( 4:25,100,500,1000 )
for ( i in 1:length(k) )
    print( c( k[i], solve_equation(k[i]) ))

## -----------------------------------------------------------------------------
em<-function(lambda,max.it=10000,eps=1e-6){
  i<-1
  lambda1<-lambda
  lambda2<-10/(3.75+3*(1/lambda1+1))
  while(abs(lambda1-lambda2)>=eps){
    lambda1<-lambda2
    lambda2<-10/(3.75+3*(1/lambda1+1))
    if(i==max.it) break
    i<-i+1
  }
  return(lambda2)
}
em(lambda=0.9)
em(lambda=1.1)

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
result.ex3.lapply <- lapply(formulas, function(x) lm(formula = x, data = mtcars))
sapply(result.ex3.lapply, rsq)

result.ex3.loops <- vector("list", length(formulas))
for (i in seq_along(formulas)){
  result.ex3.loops[[i]] <- lm(formulas[[i]], data = mtcars)
}
sapply(result.ex3.loops, rsq)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) { rows <- sample(1:nrow(mtcars), rep = TRUE) 
mtcars[rows, ]
})
result.ex4.lapply <- lapply(bootstraps, lm, formula = mpg ~ disp)
sapply(result.ex4.lapply, rsq)

result.ex4.loops <- vector("list", length(bootstraps))
for (i in seq_along(bootstraps)){
  result.ex4.loops[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}
sapply(result.ex4.loops, rsq)

## -----------------------------------------------------------------------------
vapply(cars, sd, numeric(1))

## -----------------------------------------------------------------------------
vapply(iris[vapply(iris, is.numeric, logical(1))],sd, numeric(1))

## ----warning=F----------------------------------------------------------------
library(parallel)

mcsapply<-function(x,Fun){
  cl<-makeCluster(detectCores())
  results<-parSapply(cl,x,Fun)
  stopCluster(cl)
  return(results)
}

boot_lm <- function(i) {
  boot_df <- function(x) x[sample(nrow(x), rep = T), ]
  rsquared <- function(mod) summary(mod)$r.square
  rsquared(lm(mpg ~ wt + disp, data = boot_df(mtcars)))
  }

n<-1e4
system.time(sapply(1:n,boot_lm))
system.time(mcsapply(1:n,boot_lm))

