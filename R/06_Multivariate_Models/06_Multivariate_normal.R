# Libraries required
require(QRM)
require(xts)
require(mvtnorm)

# The bivariate normal distribution
BiDensPlot(func = dmnorm, xpt=c(-4,4), ypts=c(-4,4),mu = c(0, 0), Sigma = equicorr(2, -0.7))

# Some trivariate simulated data
rho = 0.7
P <- matrix(rho,ncol=3,nrow=3)+(1-rho)*diag(3)
P
data <- rmvnorm(1000, sigma = P)

# Re-estimate the parameters
fit.norm(data)

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
load("DJ_const.rda")
Sdata <- DJ_const['2000-01-01/',1:10]
Xdata <- diff(log(Sdata))[-1,]
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

# Talking point: data which are only marginally normal
data = rAC("gumbel", n=1000, d=3, theta=3)
data = apply(data,2,qnorm)
pairs(data)
apply(data,2,shapiro.test)
MardiaTest(data)
jointnormalTest(data)
