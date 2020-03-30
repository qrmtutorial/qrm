## By Alexander J. McNeil and Marius Hofert

## Factor modelling of the yield curve using the PCA approach


## Setup
library(xts)
library(qrmtools)
library(qrmdata)

## Data preparation
data(ZCB_CAD)
ZCB10yr <- ZCB_CAD['2002-01-02/2011-12-30'] # "zero-yields" object with 10 years of data
X <- returns(ZCB10yr, method = "diff") # risk-factor changes

## PCA of daily changes in yields
PCA <- prcomp(X)
summary(PCA)
plot(PCA, xlab = "Principal component"); box() # show important components

## Extracting all information
Gamma <- PCA$rotation # principal axes (jth column is orthonormal eigenvector of cov(X) corresponding to jth largest eigenvalue) or 'loadings'
mu <- PCA$center # estimated centers
Y <- PCA$x # estimated principal components of X or 'scores'; (2012, 10)-matrix
var <- PCA$sdev^2 # explained variances per principal component
## Note: var equals the sorted eigenvalues of Cov(X) since diag(<sorted sigma^2>)
##       = Cov(Y) = Cov(Gamma^T (X - mu)) = Gamma^T Cov(X) Gamma = diag(<sorted lambda>)

## Proportion of variability explained by the first three principal axes
npr <- 3
prop <- cumsum(var)/sum(var)
prop[npr] # => ~= 97.46% of the variance is explained by the first npr principal axes/loadings

## Loadings (= principal axes) of the first three principal axes (= first three eigenvectors)
Gamma[,1:3] # (could be plotted as component samples, so each column)

## Pick out the first npr-many principal axes/loadings and plot each
## of them over time
Gamma. <- Gamma[,seq_len(npr)] # first npr-many principal axes/loadings
maturities <- seq_len(ncol(X))/4
plot(NA, xlim = range(maturities), ylim = range(Gamma.),
     xlab = "Time to maturity (years)",
     ylab = "Principal axes value (loading)")
abline(h = 0)
for(j in seq_len(npr))
    lines(maturities, Gamma.[,j], lty = j+1)
legend("topright", bty = "n", lty = 1:(npr+1),
       legend = c("Reference line", paste("PC", seq_len(npr))))

## Build time series of the first npr-many principal components of X and plot them
Y. <- xts(Y[,seq_len(npr)], time(X)) # factor series
plot.zoo(Y., type = "h", xlab = "Time", ylab = paste("Component", seq_len(npr)),
         main = "Principal components of X", panel = function(x, ...) {
             lines(x, ...)
             abline(v = as.Date("2008-09-15"), col = "maroon3") # Lehman Brothers bankruptcy
         })

## Checks
stopifnot(all.equal(cor(Y.), diag(3), check.attributes = FALSE)) # Cor ~= identity
acf(Y.) # => uncorrelated
acf(abs(Y.)) # ... but not independent
