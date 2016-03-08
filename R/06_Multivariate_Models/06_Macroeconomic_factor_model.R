## By Marius Hofert and Alexander McNeil

## Fitting a factor model

library(xts) # for time series manipulation
library(qrmdata) # for Dow Jones (constituents) data
library(qrmtools) # for plot_NA() and plot_matrix()

## Load and extract the data we work with (all available since 1990) and plot
## Index
data(DJ) # index data
DJ. <- DJ['1990-01-01/'] # all since 1990
plot.zoo(DJ., xlab="Time t", ylab="Dow Jones Index")
## Constituents
data(DJ_const) # constituents data
DJ.const <- DJ_const['1990-01-01/',] # all since 1990
plot_NA(DJ.const) # => use all but the two columns with lots of NAs
DJ.const <- DJ.const[, colSums(is.na(DJ.const)) <= 0.1 * nrow(DJ.const)] # omit columns with more than 10% NA
DJ.const <- na.fill(DJ.const, fill="extend") # fill the remaining NAs
plot.zoo(DJ.const, xlab="Time t", main="Dow Jones Constituents")

## Build and plot log-returns
## Index
X. <- diff(log(DJ.))[-1,] # compute -log-returns
plot.zoo(X., xlab="Time t", ylab=expression(X[t]), main="Risk-factor changes (log-returns) of Dow Jones index")
## Constituents
X.const <- diff(log(DJ.const))[-1,] # compute -log-returns
if(FALSE) # more time-consuming
    pairs(as.matrix(X.const), gap=0, pch=".",
          main="Scatter plot matrix of risk-factor changes (log-returns) of Dow Jones constituents")
plot.zoo(X.const, xlab="Time t", main="Risk-factor changes (log-returns) of Dow Jones constituents")

## We use monthly data here as basis (and compute monthly log-returns)
X <- apply.monthly(X.const, FUN=colSums) # (312, 28)-matrix
plot.zoo(X, type="h", xlab="Time t", main="Monthly risk-factor changes (log-returns) of Dow Jones constituents")
F <- apply.monthly(X., FUN=colSums)
plot(F, type="h", xlab="Time t", ylab=expression(X[t]), main="Monthly risk-factor changes (log-returns) of Dow Jones index")

## Fit a multivariate regression model X = a + B*F + eps
(res <- lm(X ~ F)) # more details via summary()

## Get parameter estimates
names(res)
par.ests <- coefficients(res)
a <- par.ests[1,]
B <- par.ests[2,]

## Get residuals and estimate their correlation matrix
eps <- resid(res) # (312, 28)-matrix
cor.eps <- cor(eps)

## Is Cor(eps) (roughly) diagonal?
plot_matrix(cor.eps, at=seq(-1, 1, length.out=200)) # => yes (as required)

## Are the errors uncorrelated with the factors?
cor.eps.F <- cor(eps, F) # 28 correlations (idiosyncratic risk)
summary(cor.eps.F) # => yes (as required)

## Construct the implied covariance and correlation matrix
Ups <- cov(eps) # Upsilon (covariance matrix of epsilon)
Omega <- as.matrix(cov(F)) # Omega (covariance matrix of F; systematic risk)
Sigma <- B %*% Omega %*% t(B) + diag(diag(Ups)) # Cov(X); (28, 28)-matrix
rownames(Sigma) <- colnames(Sigma)
P <- cov2cor(Sigma) # Cor(X)

## Look at discrepancies between the factor model correlation matrix and the
## sample correlation matrix
err <- P-cor(X)
plot_matrix(err, at=seq(-1, 1, length.out=200))

