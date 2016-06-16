## by Alexander McNeil

# Libraries required
library(QRM)
library(xts)
library(mvtnorm)
library(qrmdata)

# The bivariate normal distribution
BiDensPlot(func = dmnorm, xpt=c(-4,4), ypts=c(-4,4),mu = c(0, 0), Sigma = equicorr(2, -0.7))

# Some trivariate simulated data
rho = 0.7
P <- matrix(rho,ncol=3,nrow=3)+(1-rho)*diag(3)
P
data <- rmvnorm(1000, sigma = P)

# ML estimation
n <- dim(data)[1]
mu.hat <- colMeans(data)
Sigma.hat <- var(data)*(n-1)/n
P.hat <- cor(data)
# Note that var(data) gives unbiased estimator which is not MLE
ll.max <- sum(dmvnorm(data, mean = mu.hat, sigma = Sigma.hat, log=TRUE))
ll.max

# Some 10-d simulated data
rho = 0.6
P <- matrix(rho,ncol=10,nrow=10)+(1-rho)*diag(10)
P
data <- rmvnorm(1000, sigma = P)

# Test of marginal normality
apply(data,2,shapiro.test)

# Tests of joint normality
MardiaTest(data)
jointnormalTest(data)

# Dow Jones Data
data("DJ_const")
Sdata <- DJ_const['2000-01-01/',1:10]
Xdata <- diff(log(Sdata))[-1,]
plot.zoo(Xdata)
Xdata.w <- apply.weekly(Xdata,FUN=colSums)
Xdata.m <- apply.monthly(Xdata,FUN=colSums)
Xdata.q <- apply.quarterly(Xdata,FUN=colSums)

# Test of marginal normality
apply(Xdata,2,shapiro.test)
apply(Xdata.w,2,shapiro.test)
apply(Xdata.m,2,shapiro.test)
apply(Xdata.q,2,shapiro.test)

# Test for joint normality
MardiaTest(as.matrix(Xdata))
jointnormalTest(as.matrix(Xdata))

MardiaTest(as.matrix(Xdata.w))
jointnormalTest(as.matrix(Xdata.w))

MardiaTest(as.matrix(Xdata.m))
jointnormalTest(as.matrix(Xdata.m))

MardiaTest(as.matrix(Xdata.q))
jointnormalTest(as.matrix(Xdata.q))

# Constructing data which are only marginally normal

library(copula)
set.seed(133)
G.cop  <- archmCopula("Gumbel", param=3,  dim=3) 
U.G  <- rCopula(1000, copula=G.cop)
data = qnorm(U.G)
pairs(data)
apply(data,2,shapiro.test)
MardiaTest(data)
jointnormalTest(data)

