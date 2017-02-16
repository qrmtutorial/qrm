## By Marius Hofert and Alexander McNeil

## Generate and transform multivariate data (with the probability and quantile
## transformations)


### 0 Setup ####################################################################

library(qrmtools) # for qPar()
library(copula) # for pairs2()
library(qrmdata) # for S&P 500 data etc.

set.seed(271) # for reproducibility


### 1 Normal copula sample #####################################################

## Sample from a bivariate normal distribution
P <- matrix(0.7, nrow = 2, ncol = 2) # correlation matrix; play with the entry!
diag(P) <- 1
A <- t(chol(P)) # Cholesky factor
n <- 1000 # sample size
Z <- matrix(rnorm(n * 2), ncol = 2) # iid N(0,1)
X <- Z %*% t(A) # X = AZ ~ N(0, P)
plot(X, xlab = expression(X[1]), ylab = expression(X[2])) # scatter plot

## (Negative) logarithmic stock (and other) returns can look quite similarly
data(SP500_const) # S&P 500 data
time <- c("2007-01-03", "2009-12-31") # time period
dat <- SP500_const[paste0(time, collapse = "/"), c("AAPL", "IBM")] # grab out stocks
X. <- as.matrix(-log_returns(dat)) # build negative logarithmic returns
plot(X., xlab = expression("Apple -log-returns"~X[1]),
     ylab = expression("IBM -log-returns"~X[2])) # scatter plot

## Marginally apply probability transforms (with the corresponding dfs)
U <- pnorm(X) # probability transformation
plot(U, xlab = expression(U[1]), ylab = expression(U[2])) # sample from the Gauss copula

## ... and the pseudo-sample of the underlying copula for the stock data
plot(pobs(X.), xlab = expression(U[1]), ylab = expression(U[2]))

## (Visually) check that the margins of our generated data are indeed U(0,1)
plot(U[,1], ylab = expression(U[1])) # scatter plot
hist(U[,1], xlab = expression(U[1]), probability = TRUE, main = "") # histogram
plot(U[,2], ylab = expression(U[2])) # scatter plot
hist(U[,2], xlab = expression(U[2]), probability = TRUE, main = "") # histogram

## Map the Gauss copula sample to one with t_3.5 margins
nu <- 3.5 # degrees of freedom
Y <- qt(U, df = nu) # quantile transformation
plot(Y, xlab = expression(Y[1]), ylab = expression(Y[2]), col = "royalblue3") # scatter plot
points(X)
legend("bottomright", bty = "n", pch = c(1,1), col = c("black", "royalblue3"),
       legend = c(expression(N(0,P)),
                  substitute("With"~italic(t)[nu.]~"margins", list(nu. = nu))))
## => Sample is stretched out by heavier tailed marginal distributions


### 2 t copula sample ##########################################################

## Build a block correlation matrix
d <- 5
rho <- c(0.3, 0.6, 0.9)
P <- matrix(rho[1], nrow = d, ncol = d)
P[1:2, 1:2] <- rho[2]
P[3:5, 3:5] <- rho[3]
diag(P) <- 1
P

## Sample from a multivariate t distribution
A <- t(chol(P)) # Cholesky factor
Z <- matrix(rnorm(n * d), ncol = d) # iid N(0,1)
sqrt.W <- sqrt(1/rgamma(n, shape = nu/2, rate = nu/2))
X <- sqrt.W * (Z %*% t(A)) # X = sqrt(W)AZ ~ t_{P,nu}
pairs2(X, labels.null.lab = "X", cex = 0.4, col = adjustcolor("black", alpha.f = 0.5))

## Marginally apply probability transforms (here: empirically estimated dfs)
U <- pobs(X) # probability transformation with the empirically estimated margins
pairs2(U, cex = 0.4, col = adjustcolor("black", alpha.f = 0.5))

## Map the t copula sample to one with N(0,1) and Par(3) margins
Y <- cbind(qPar(U[,1], theta = 3), qnorm(U[,2:5])) # quantile transformation
pairs2(Y, cex = 0.4, col = adjustcolor("black", alpha.f = 0.5),
       labels.null.lab = "Y")
## ... and the components have the same dependence structure as those of X above!