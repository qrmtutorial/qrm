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
X. <- diff(log(DJ.))[-1,] # compute -log-returns
plot.zoo(X., xlab = "Time t", ylab = expression(X[t]),
         main = "Risk-factor changes (log-returns) of Dow Jones index")
## Constituents
X.const <- diff(log(DJ.const))[-1,] # compute -log-returns
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
PCA <- princomp(X) # determine the principal components
PCA. <- prcomp(X) # another option
summary(PCA)
plot(PCA) # show important components
loadings(PCA) # factor loadings (representing how much a factor explains a variable)

## Working with the principal components
mu <- PCA$center # estimated centers
Y <- PCA$scores # estimated principal components of X
nprin <- 3 # number of important principal components
Y1 <- xts(Y[,1:nprin], time(X)) # grab out the first so-many important principal components of X
plot.zoo(Y1, type = "h", xlab = "Time", ylab = paste("Component", 1:nprin),
         main = "Principal components of X")
cor(Y1)
Y2 <- Y[,(nprin+1):ncol(Y)] # ignored principal components of X

## Interpreting first principal component
plot.zoo(merge(Xindex,Y1[,1]),main=" ") # compare DJ index returns and first PC
plot(as.numeric(Xindex),-as.numeric(Y1[,1]),xlab="Index",ylab="PC1 (neg)")

## Reconstruct X from Y1 and Y2
G <- unclass(loadings(PCA)) # compute G, see McNeil, Frey, Embrechts (2015, (6.64))
G1 <- G[,1:nprin] # compute G_1, see McNeil, Frey, Embrechts (2015, (6.65))
G2 <- G[,(nprin+1):ncol(G)]
eps <- G2 %*% t(Y2) # compute epsilon
eps.cor <- cor(t(eps)) # compute correlations of epsilons
diag(eps.cor) <- NA # neglect diagonal elements
matrix_plot(eps.cor) # not perfectly diagonal

X. <- t(mu + G1 %*% t(Y1) + eps)
err <- X.-X
head(err) # check that differences are zero vectors
summary(as.vector(err))

## Note: We could fit a univariate GARCH model to each PCA factor.
##       This is called PC-GARCH model.
