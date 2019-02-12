## By Marius Hofert and Alexander McNeil

## Fitting symmetric and skewed t, normal inverse Gaussian (NIG) and generalized
## hyperbolic (GH) models to selected Dow Jones constituents


### Setup ######################################################################

library(xts)
library(mvtnorm)
library(QRM)
library(qrmtools)
library(qrmdata)
library(ghyp)


### 1 Data preparation #########################################################

## Dow Jones constituents we work with
data("DJ_const")
d <- 10 # for the number of components we pick out
S <- DJ_const['2000-01-01/',1:d] # first d components

## Daily log-returns
X.d <- returns(S)
plot.zoo(X.d)

## Weekly log-returns
X.w <- apply.weekly(X.d, FUN = colSums)
plot.zoo(X.w)


### 2 Is the return data jointly normal? #######################################

## Formal tests of normality

## Daily log-returns
apply(X.d, 2, function(x) shapiro.test(x)$p.value) # tests of marginal normality
jointnormalTest(as.matrix(X.d), plot = FALSE) # Kolmogorov--Smirnov test of squared Mahalanobis distances being chi_d^2
MardiaTest(as.matrix(X.d)) # Mardia's test of joint normality

## Weekly log-returns
apply(X.w, 2, function(x) shapiro.test(x)$p.value) # tests of marginal normality
jointnormalTest(as.matrix(X.w), plot = FALSE) # Kolmogorov--Smirnov test of squared Mahalanobis distances being chi_d^2
MardiaTest(as.matrix(X.w)) # Mardia's test of joint normality

## => All p-values << 0.05 => The data is certainly not normally distributed.


## Visual tests of normality

## Daily log-returns
pairs(as.matrix(X.d), gap = 0, pch = ".") # visual assessment
qq_plot(X.d[,1], FUN = qnorm, method = "empirical") # first margin only => already not normal
D2.d <- mahalanobis(X.d, center = colMeans(X.d), cov = cov(X.d)) # squared Mahalanobis distances
qq_plot(D2.d, FUN = function(p) qchisq(p, df = d)) # => departure clearly visible

## Weekly log-returns
pairs(as.matrix(X.w), gap = 0, pch = ".") # visual assessment
qq_plot(X.w[,1], FUN = qnorm, method = "empirical") # better but still not too good
D2.w <- mahalanobis(X.w, center = colMeans(X.w), cov = cov(X.w)) # squared Mahalanobis distances
qq_plot(D2.w, FUN = function(p) qchisq(p, df = d)) # => departure clearly visible


### 3 Fit various multivariate models ##########################################

max.iter <- 1e4 # maximal number of iterations for the fitting procedures

## Daily log-returns

## Symmetric
fit.t.sym.d   <- fit.tmv   (X.d, symmetric = TRUE, nit = max.iter, silent = TRUE) # t
fit.NIG.sym.d <- fit.NIGmv (X.d, symmetric = TRUE, nit = max.iter, silent = TRUE) # normal inverse Gaussian
fit.GH.sym.d  <- fit.ghypmv(X.d, symmetric = TRUE, nit = max.iter, silent = TRUE) # generalized hyperbolic
## Skewed
fit.t.skw.d   <- fit.tmv   (X.d, symmetric = FALSE, nit = max.iter, silent = TRUE) # t
fit.NIG.skw.d <- fit.NIGmv (X.d, symmetric = FALSE, nit = max.iter, silent = TRUE) # normal inverse Gaussian
fit.GH.skw.d  <- fit.ghypmv(X.d, symmetric = FALSE, nit = max.iter, silent = TRUE) # generalized hyperbolic
## Note: Warnings are due to deprecated recycling of an array of length 1

fit.GH.skw.d@gamma # skewness parameters

## Comparison of log-likelihoods
likelihoods.d <- c(fit.t.sym.d@llh, fit.NIG.sym.d@llh, fit.GH.sym.d@llh,
                   fit.t.skw.d@llh, fit.NIG.skw.d@llh, fit.GH.skw.d@llh)
which.max(likelihoods.d) # => skewed generalized hyperbolic

## Likelihood-ratio tests (for comparing the goodness-of-fit of two models)
## H0: symmetric t is no worse than skewed t
LRstat <- 2 * (fit.t.skw.d@llh - fit.t.sym.d@llh) # test statistic
1 - pchisq(LRstat, 10) # => H0 not rejected at 5% level
## H0: symmetric generalized hyperbolic is no worse than skewed generalized hyperbolic
LRstat <- 2 * (fit.GH.skw.d@llh - fit.GH.sym.d@llh) # test statistic
1 - pchisq(LRstat, 10) # => H0 not rejected at 5% level
## H0: symmetric t is no worse than symmetric generalized hyperbolic
LRstat <- 2 * (fit.GH.sym.d@llh - fit.t.sym.d@llh) # test statistic
1 - pchisq(LRstat, 1) # => H0 rejected at 5% level => symmetric generalized hyperbolic preferred


## Weekly log-returns

## Symmetric
fit.t.sym.w   <- fit.tmv   (X.w, symmetric = TRUE, nit = max.iter, silent = TRUE) # t
fit.NIG.sym.w <- fit.NIGmv (X.w, symmetric = TRUE, nit = max.iter, silent = TRUE) # normal inverse Gaussian
fit.GH.sym.w  <- fit.ghypmv(X.w, symmetric = TRUE, nit = max.iter, silent = TRUE) # generalized hyperbolic
## Skewed
fit.t.skw.w   <- fit.tmv   (X.w, symmetric = FALSE, nit = max.iter, silent = TRUE) # t
fit.NIG.skw.w <- fit.NIGmv (X.w, symmetric = FALSE, nit = max.iter, silent = TRUE) # normal inverse Gaussian
fit.GH.skw.w  <- fit.ghypmv(X.w, symmetric = FALSE, nit = max.iter, silent = TRUE) # generalized hyperbolic

## Comparison of log-likelihoods
likelihoods.w <- c(fit.t.sym.w@llh, fit.NIG.sym.w@llh, fit.GH.sym.w@llh,
                   fit.t.skw.w@llh, fit.NIG.skw.w@llh, fit.GH.skw.w@llh)
which.max(likelihoods.w) # => also skewed generalized hyperbolic

## Likelihood-ratio tests (for comparing the goodness-of-fit of two models)
## H0: symmetric t is no worse than skewed t
LRstat <- 2 * (fit.t.skw.w@llh - fit.t.sym.w@llh) # test statistic
1 - pchisq(LRstat, 10) # => H0 not rejected at 5% level
## H0: symmetric generalized hyperbolic is no worse than skewed generalized hyperbolic
LRstat <- 2 * (fit.GH.skw.w@llh - fit.GH.sym.w@llh) # test statistic
1 - pchisq(LRstat, 10) # => H0 not rejected at 5% level
## H0: symmetric t is no worse than symmetric generalized hyperbolic
LRstat <- 2 * (fit.GH.sym.w@llh - fit.t.sym.w@llh) # test statistic
1 - pchisq(LRstat, 1) # => H0 rejected at 5% level => symmetric generalized hyperbolic preferred


### 4 Visual assessment ########################################################

X <- as.matrix(X.d) # data we consider

## Visually compare two component samples against the corresponding pairwise
## implied contour lines
ind <- c(4, 5) # "CAT", "CSCO"
X. <- X[,ind]
plot(X.)
obj <- fit.GH.sym.d # get d-dimensional object
obj@mu <- obj@mu[ind]
obj@sigma <- obj@sigma[ind, ind]
obj@gamma <- obj@gamma[ind[2]] # (0 in symmetric case)
x. <- seq(min(X.[,1]), max(X.[,1]), length.out = 30) # x-locations of density evaluation points
y. <- seq(min(X.[,2]), max(X.[,2]), length.out = 30) # y-locations of density evaluation points
z. <- outer(x., y., function(xx, yy)
    dghyp(cbind(xx, yy), object = obj)) # implied bivariate density
contour(x., y., z = z., drawlabels = FALSE, axes = FALSE, col = "royalblue3", add = TRUE)

## Similarly for all d pairs
row <- 1 # row number in upper-right panel
col <- 1 # column number in upper-right panel
pairs(X, gap = 0, pch = ".",
      upper.panel = function(x, y, ...) {
          ## Figure out in which row and column we are in the plot (needed later)
          row <<- row + 1
          if(row >= col) {
              col <<- col + 1
              row <<- 1
          }
          ## Plot the points with indices (i, j)
          points(x, y, ...)
          ## Extract the (i, j)th 'margin' of the object of class 'ghyp'
          obj <- fit.GH.sym.d # get d-dimensional object
          ind <- c(row, col) # pairwise indices to consider
          obj@mu <- obj@mu[ind]
          obj@sigma <- obj@sigma[ind, ind]
          obj@gamma <- obj@gamma[col] # (0 in symmetric case)
          ## Compute contours based on
          x. <- seq(min(x), max(x), length.out = 30) # x-locations of density evaluation points
          y. <- seq(min(y), max(y), length.out = 30) # y-locations of density evaluation points
          z. <- outer(x., y., function(xx, yy)
              dghyp(cbind(xx, yy), object = obj)) # implied bivariate density
          ## Plot contours
          contour(x., y., z = z., nlevels = 5,
                  drawlabels = FALSE, axes = FALSE, col = "royalblue3", add = TRUE)
      })
## Note: For more than a handful of components, one can use a 'zenplot'
