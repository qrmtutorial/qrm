# Library required
# require(qrmdata)
library(xts)
library(mvtnorm)
library(qrmdata)

# We will use DJ stock to obtain multivariate data
data("DJ_const")
Sdata <- DJ_const['2000-01-01/',1:4]
Xdata <- diff(log(Sdata))[-1,]
plot.zoo(Xdata)

# First we will estimate the moments of these data
# sample mean vector, sample covariance and correlation matrices
colMeans(Xdata)
cor(Xdata)
cov(Xdata)
var(Xdata)
Sigma <- var(Xdata)
mu <- colMeans(Xdata)
P <- cor(Xdata)

# Cholesky decomposition
chol(Sigma)
A <- t(chol(Sigma))
A %*% t(A)

# Symmetric decomosition
# This is an alternative way of finding a matrix 'square root'
eigen(Sigma)
Gamma <- eigen(Sigma)$vectors
Lambda <- eigen(Sigma)$values
Lambda
Lambda <- diag(sqrt(Lambda))
Lambda
B <- Gamma %*% Lambda %*% t(Gamma)
B %*% B
Sigma



# simulation mutivariate normal data
simdata <- rmvnorm(n=1000,mean=mu,sigma=Sigma)
colMeans(simdata)
mu
var(simdata)
Sigma
cor(simdata)
P


