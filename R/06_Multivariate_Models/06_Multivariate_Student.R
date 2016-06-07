## by Alexander McNeil

# Libraries required
library(xts)
library(mvtnorm)
library(QRM)
library(qrmdata)

# Simulate multivariate Student
rho <- 0.7
d <- 3
P <- matrix(rho,ncol=d,nrow=d) + (1-rho)*diag(d)
nu <- 10
data <- rmvt(2000, sigma = P, df=nu)
pairs(data)

# Test of marginal normality
apply(data,2,shapiro.test)

# Tests of joint normality
MardiaTest(data)
jointnormalTest(data)

# QQplots against Student t
qplotfunc <- function(data){
  qqplot(qt(ppoints(data),df=nu),data,xlab="Theoretical",ylab="Empirical")
  qqline(data,dist=function(p){qt(p,df=nu)})
}
par(mfrow=c(1,3))
apply(data,2,qplotfunc)
par(mfrow=c(1,1))

# QQplot of Mahalanobis distances (scaled by dimension d)
# against an F(d,nu) distribution
Ddata <- mahalanobis(data,center=FALSE,cov=var(data))/d
qqplot(qf(ppoints(Ddata),df1=d,df2=nu),Ddata,xlab="Theoretical",ylab="Empirical")
qqline(Ddata,dist=function(p){qf(p,df1=d,df2=nu)})

# Estimate multivariate Student

# Dow Jones Data
data("DJ_const")
Sdata <- DJ_const['2000-01-01/',1:10]
Xdata <- diff(log(Sdata))[-1,]
Xdata.w <- apply.weekly(Xdata,FUN=colSums)
Xdata.m <- apply.monthly(Xdata,FUN=colSums)

# write function for fitting multivariate normal by MLE
fit.mvnorm <- function(data){
  n <- dim(data)[1]
  mu.hat <- colMeans(data)
  Sigma.hat <- var(data)*(n-1)/n
  ll.max <- sum(dmvnorm(data, mean = mu.hat, sigma = Sigma.hat, log=TRUE))
  list(mu=mu.hat, Sigma=Sigma.hat, ll.max=ll.max)
}

# Fit normal
mod.w.norm = fit.mvnorm(as.matrix(Xdata.w))
mod.w.norm
# Fit Student t
mod.w.t = fit.mst(as.matrix(Xdata.w),method="BFGS")
mod.w.t
names(mod.w.t)
mod.w.t$df
# Compare
mod.w.t$ll.max - mod.w.norm$ll.max

# Fit normal
mod.m.norm = fit.mvnorm(as.matrix(Xdata.m))
mod.m.norm
# Fit Student t
mod.m.t = fit.mst(as.matrix(Xdata.m),method="BFGS")
mod.m.t
# Compare
mod.m.t$ll.max - mod.m.norm$ll.max

