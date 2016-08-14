## By Marius Hofert and Alexander McNeil

## Basic copula estimation (2d example, 4 copulas)


### 0 Setup ####################################################################

library(copula)
library(qrmtools)
library(qrmdata)
library(xts)


### 1 Working with the data ####################################################

## Load the time series
data("SP500")
data("FTSE")

## Build negative log-returns
X.SP500 <- -log_returns(SP500)
X.FTSE  <- -log_returns(FTSE)

## Merge
X <- merge(X.SP500, X.FTSE, all = FALSE)
time <- c("2003-01-01", "2012-12-31")
X <- X[paste0(time, collapese = "/"),]

## Basic plot
plot.zoo(X)

## Observe that there are zero values caused by market closures
apply(X == 0, 2, sum)
rexcl <- (X[,1] == 0) | (X[,2] == 0)
X. <- X[!rexcl,]

## Aggregating by week
dim(X.w <- apply.weekly(X., FUN = colSums))
plot.zoo(X.w, type = "h")

## Compute pseudo-observations
U <- as.matrix(pobs(X.w))
plot(U, xlab = expression(U[1]), ylab = expression(U[2]))


### 2 Fitting copulas ##########################################################

## Compare various bivariate copulas
fit.N <- fitCopula(normalCopula(),  data = U)
fit.t <- fitCopula(tCopula(),       data = U) # df of freedom are estimated, too
fit.C <- fitCopula(claytonCopula(), data = U)
fit.G <- fitCopula(gumbelCopula(),  data = U)

## Comparing the likelihoods
sort(c(N = fit.N@loglik, t = fit.t@loglik, C = fit.C@loglik, G = fit.G@loglik),
     decreasing = TRUE)

