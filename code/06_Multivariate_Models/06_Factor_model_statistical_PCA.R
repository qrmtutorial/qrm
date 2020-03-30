## By Marius Hofert and Alexander McNeil

## Fitting a PCA factor model


### Setup ######################################################################

library(xts) # for time series manipulation
library(qrmdata) # for Dow Jones (constituents) data
library(qrmtools) # for plot_NA() and plot_matrix()


### 1 Data preparation #########################################################

## Load and extract the data we work with (all available since 1990) and plot
## Index
data(DJ) # index data
DJ. <- DJ['1990-01-01/'] # all since 1990
plot.zoo(DJ., xlab = "Time t", ylab = "Dow Jones Index")
## Constituents
data(DJ_const) # constituents data
DJ.const <- DJ_const['1990-01-01/',] # all since 1990
NA_plot(DJ.const) # => use all but the two columns with lots of NAs
DJ.const <- DJ.const[, colSums(is.na(DJ.const)) <= 0.1 * nrow(DJ.const)] # omit columns with more than 10% NA
DJ.const <- na.fill(DJ.const, fill = "extend") # fill the remaining NAs
plot.zoo(DJ.const, xlab = "Time t", main = "Dow Jones Constituents")

## Build and plot log-returns
## Index
X. <- returns(DJ.) # compute log-returns
plot.zoo(X., xlab = "Time t", ylab = expression(X[t]),
         main = "Risk-factor changes (log-returns) of Dow Jones index")
## Constituents
X.const <- returns(DJ.const) # compute log-returns
if(FALSE) # more time-consuming
    pairs(as.matrix(X.const), gap = 0, pch = ".",
          main = "Scatter plot matrix of risk-factor changes (log-returns) of Dow Jones constituents")
plot.zoo(X.const, xlab = "Time t", main = "Risk-factor changes (log-returns) of Dow Jones constituents")

## We use monthly data here as basis (and compute monthly log-returns)
X <- apply.monthly(X.const, FUN = colSums) # (312, 28)-matrix
plot.zoo(X, type = "h", xlab = "Time t",
         main = "Monthly risk-factor changes (log-returns) of Dow Jones constituents")
Xindex <- apply.monthly(X., FUN = colSums)
plot(Xindex, type = "h", xlab = "Time t", ylab = expression(X[t]),
     main = "Monthly risk-factor changes (log-returns) of Dow Jones index")


### 2 Model fitting ############################################################

## Principal component analysis
PCA <- prcomp(X) # another option is princomp() (but numerically not as robust)
summary(PCA)
plot(PCA, xlab = "Principal component"); box() # show important components

## Extracting all information
Gamma <- PCA$rotation # principal axes (jth column is orthonormal eigenvector of cov(X) corresponding to jth largest eigenvalue) or 'loadings'
mu <- PCA$center # estimated centers
Y <- PCA$x # estimated principal components of X or 'scores'; (2012, 10)-matrix
var <- PCA$sdev^2 # explained variances per principal component
## Note: var equals the sorted eigenvalues of Cov(X) since diag(<sorted sigma^2>)
##       = Cov(Y) = Cov(Gamma^T (X - mu)) = Gamma^T Cov(X) Gamma = diag(<sorted lambda>)

## Working with the principal components
npr <- 3 # number of important principal components
Y1 <- xts(Y[,1:npr], time(X)) # grab out the first so-many principal components of X
stopifnot(all.equal(cor(Y1), diag(3), check.attributes = FALSE)) # check
plot.zoo(Y1, type = "h", xlab = "Time", ylab = paste("Component", 1:npr),
         main = "Principal components of X") # plot of the first so-many principal components of X

## Interpreting first principal component
plot.zoo(merge(Xindex, Y1[,1]), main = " ") # compare DJ index returns and first PC
plot(as.numeric(Xindex), -as.numeric(Y1[,1]), xlab = "Index", ylab = "-PC1")

## Part Y2 of (Y1, Y2) = Y
Y2 <- Y[,(npr+1):ncol(Y)] # ignored principal components of X

## Reconstruct X from Y1 and Y2
G1 <- Gamma[,1:npr] # compute G_1, see MFE (2015, (6.65))
G2 <- Gamma[,(npr+1):ncol(Gamma)]
eps <- G2 %*% t(Y2) # compute epsilon
eps.cor <- cor(t(eps)) # compute correlations of epsilons
matrix_plot(eps.cor, ran = c(-1, 1)) # check; not perfectly but okay
X. <- t(mu + G1 %*% t(Y1) + eps) # put everything together; see MFE (2015, (6.62) or (6.65))
err <- X. - X # difference to original X
summary(as.vector(err)) # => indeed all very small

## Note: We could fit a univariate GARCH model to each PCA factor.
##       This is called PC-GARCH model.
