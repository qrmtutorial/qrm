## By Marius Hofert

## Fitting a multivariate normal and t distribution to negative log-returns of
## the five components of the Dow Jones index

library(xts) # for time series manipulation
library(QRM) # for fit.mst()
library(mvtnorm) # for sampling from a multivariate normal or t distribution
library(qrmdata) # for Dow Jones constituents data
library(qrmtools) # for returns()

## Load and extract the data we work with and plot
data(DJ_const)
str(DJ_const)
S <- DJ_const['2000-01-01/', c("AAPL", "BA", "INTC", "IBM", "NKE")]
## => We work with Apple, Boeing, Intel, IBM, Nike since 2000-01-01 here
plot.zoo(S, xlab = "Time t")

## Build and plot -log-returns
X <- -returns(S) # compute -log-returns
pairs(as.matrix(X), main = "Scatter plot matrix of risk-factor changes", gap = 0, pch = ".")
plot.zoo(X, xlab = "Time t", main = "")
## For illustration purposes, we treat this data as realizations of iid
## five-dimensional random vectors.

## Fitting a multivariate normal distribution to X and simulating from it
mu <- colMeans(X) # estimated location vector
Sigma <- cov(X) # estimated scale matrix
stopifnot(all.equal(Sigma, var(X)))
P <- cor(X) # estimated correlation matrix
stopifnot(all.equal(P, cov2cor(Sigma)))
n <- nrow(X) # sample size
set.seed(271)
X.norm <- rmvnorm(n, mean = mu, sigma = Sigma) # N(mu, Sigma) samples

## Fitting a multivariate t distribution to X
fit <- fit.mst(X, method = "BFGS") # fit a multivariate t distribution
X.t <- rmvt(n, sigma = as.matrix(fit$Sigma), df = fit$df, delta = fit$mu) # t_nu samples

## Plot (sample from fitted t (red), original sample (black), sample from fitted normal (blue))
dat <- rbind(t = X.t, original = as.matrix(X), norm = X.norm)
cols <- rep(c("maroon3", "black", "royalblue3"), each = n)
pairs(dat, gap = 0, pch = ".", col = cols)

## Pick out one pair (to better see that multivariate t fits better)
plot(dat[,1:2], col = cols) # => the multivariate normal generates too few extreme losses!
legend("bottomright", bty = "n", pch = rep(1, 3), col = c("black", "royalblue3", "maroon3"),
       legend = c("-Log-returns", expression("Simulated fitted"~N(mu,Sigma)),
                  expression("Simulated fitted"~italic(t)[nu](mu,Sigma))))
