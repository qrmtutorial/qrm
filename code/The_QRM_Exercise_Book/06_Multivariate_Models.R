## By Marius Hofert

## R script for Chapter 6 of The QRM Exercise Book


### Setup ######################################################################

library(nvmix) # for rNorm(), rStudent()
library(copula)
library(xts)
library(qrmtools)
library(qrmdata)
library(ghyp) # for generalized hyperbolic distributions

doPDF <- require(crop)
options(digits = 10)


### Exercise 6.5 (Distinguishing distributions by scatterplots of random samples)

## Reproducing code

## Parameters
n <- 1000 # sample size
d <- 2 # dimension
tau <- 0.5 # Kendall's tau

## a) N(0, P)
rho <- iTau(ellipCopula("normal"), tau = tau)
P <- p2P(rho, d = d)
set.seed(271)
X.N <- rNorm(n, scale = P)

## b) t_nu(0, P) with negative P_{12}
nu <- 3.5
P. <- p2P(-rho, d = d)
set.seed(271)
X.t.m <- rStudent(n, df = nu, scale = P.)

## c) t_nu(0, P)
##    Note: same 'rho' (and thus P) as above
set.seed(271)
X.t <- rStudent(n, df = nu, scale = P)

## d) t_nu(0, I)
set.seed(271)
X.t.I <- rStudent(n, df = nu, scale = diag(d))

## e) Elliptical distribution with discrete radial part
set.seed(271)
Z <- matrix(rnorm(n * d), ncol = d)
S <- Z / sqrt(rowSums(Z^2))
A <- t(chol(P))
R <- 1 + rbinom(n, size = 1, prob = 2/3) # P(R = 1) = 1/3, P(R = 2) = 2/3
X.dis.R <- R * t(A %*% t(S))

if(FALSE) # unused
    plot(qnorm(pobs(X.dis.R))) # nice :-)

## f) Elliptical distribution with bounded radial part
R <- 1 + runif(n) # R ~ U(1,2)
X.bdd.R <- R * t(A %*% t(S))
if(FALSE) # unused
    plot(qnorm(pobs(X.bdd.R))) # nice :-)

## g) Mixture of normal and singular distribution
set.seed(271)
bool <- sample(c(FALSE, TRUE), size = n, replace = TRUE, prob = c(2/3, 1/3))
X.mix <- X.N
X.mix[bool, 2] <- X.mix[bool, 1]

if(FALSE) { # unused
    ## MO with N(0,1) margins
    alpha <- c(0.2, 0.6)
    set.seed(271)
    V <- matrix(runif(n * 3), ncol = 3)
    U <- cbind(pmax(V[,1]^(1/(1 - alpha[1])), V[,3]^(1/alpha[1])),
               pmax(V[,2]^(1/(1 - alpha[2])), V[,3]^(1/alpha[2])))
    X.MO.N01 <- qnorm(U)
    plot(X.MO.N01)
}

## h) Non-radially symmetric distribution
X.N.mirr <- X.N
ii <- X.N.mirr[,1] < 0
X.N.mirr[ii, 2] <- -X.N.mirr[ii, 2]

if(FALSE) { # unused
    ## Gumbel copula with N(0,1) margins
    family <- "Gumbel"
    th <- iTau(archmCopula(family), tau)
    cop <- archmCopula(family, param = th, dim = d)
    set.seed(271)
    X.G.N01 <- qnorm(rCopula(n, copula = cop))
    plot(X.G.N01)
}

## i) Normal mean-variance mixture (not radially symmetric)
set.seed(271)
W <- rexp(n) # too extreme: nu / rchisq(n, df = nu)
mu <- c(1, 1) # need a mu here (mu = 0 => not a normal mean-variance mixture)
X.nmvmix <- rep(mu, each = n) * W + sqrt(W) * t(A %*% t(Z)) # of the form X = m(W) + sqrt(W) * A * Z
## Note: to see non-exchangeability, use different components in mu (e.g., c(1, 10))

## Plot
r1 <- range(X.N, X.mix, X.N.mirr)
rmx1 <- max(abs(r1)) + 0.2
ran1 <- c(-rmx1, rmx1) # symmetricize
r2 <- range(X.dis.R, X.bdd.R)
rmx2 <- max(abs(r2)) + 0.2
ran2 <- c(-rmx2, rmx2) # symmetricize

## Plot
if(doPDF) pdf(file = (file <- "fig_06_recognizing_features.pdf"), width = 9, height = 9)
opar <- par(mar = c(5, 4, 4, 2) + 0.1 - c(1, 0, 3, 1)) # # set margins (bottom, left, top, right)
lay <- matrix(1:9, ncol = 3, byrow = TRUE) # layout matrix
layout(lay) # layout
## Note: Order obtained from set.seed(271); c("N", "t", "t.m", "t.I", "dis.R", "bdd.R", "mix", "N.mirr", "nmvmix")[sample(1:9)]
plot(X.t.m,    cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]))
plot(X.dis.R,  cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran2, ylim = ran2)
plot(X.nmvmix, cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]))
plot(X.t,      cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]))
plot(X.mix,    cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran1, ylim = ran1)
plot(X.N,      cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran1, ylim = ran1)
plot(X.bdd.R,  cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran2, ylim = ran2)
plot(X.t.I,    cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]))
plot(X.N.mirr, cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran1, ylim = ran1)
par(opar)
if(doPDF) dev.off.crop(file)


### Exercise 6.18 (Contrasting models for uncorrelated normal variables) #######

## c) Reproducing code

## Define the distribution functions of X_1+X_2 and Y_1+Y_2
pXsum <- function(x) pnorm(x/sqrt(2))
pYsum <- function(x) (ifelse(x >= 0, 1, 0) + pnorm(x/2))/2

## Define VaR_alpha(X_1+X_2) and VaR_alpha(Y_1+Y_2)
qXsum <- function(alpha) sqrt(2)*qnorm(alpha)
qYsum <- function(alpha)
    sapply(alpha, function(a) {
        if(a <= 1/4) 2*qnorm(2*a)
        else if(a <= 3/4) 0
        else 2*qnorm(2*a-1)
    })

## Plot of the distribution functions
x <- seq(-4, 4, length.out = 257)
F.Xsum <- pXsum(x)
F.Ysum <- pYsum(x)
ylim <- range(F.Xsum, F.Ysum)
if(doPDF) pdf(file = (file <- "fig_06_dfs_uncorrelated_normals.pdf"))
opar <- par(pty = "s")
plot(x, F.Xsum, type = "l", ylim = ylim, xlab = "x", ylab = "")
lines(x, F.Ysum, col = "royalblue3")
legend("bottomright", bty = "n", col = c("black", "royalblue3"), lty = c(1, 1),
       legend = c(expression(F[X[1]+X[2]](x)), expression(F[Y[1]+Y[2]](x))))
par(opar)
if(doPDF) dev.off.crop(file)

## Plot of the VaRs
alpha <- seq(0.0001, 0.9999, length.out = 257)
VaR.Xsum <- qXsum(alpha)
VaR.Ysum <- qYsum(alpha)
ylim <- range(VaR.Xsum, VaR.Ysum)
if(doPDF) pdf(file = (file <- "fig_06_VaRs_uncorrelated_normals.pdf"))
opar <- par(pty = "s")
plot(alpha, VaR.Xsum, type = "l", ylim = ylim, xlab = expression(alpha), ylab = "")
lines(alpha, VaR.Ysum, col = "royalblue3")
legend("bottomright", bty = "n", col = c("black", "royalblue3"), lty = c(1, 1),
       legend = c(expression(VaR[alpha](X[1]+X[2])), expression(VaR[alpha](Y[1]+Y[2]))))
par(opar)
if(doPDF) dev.off.crop(file)


## d)

## Values
alpha <- c(0.9, 0.99)
qXsum(alpha) # 1.812387605 3.289952714
qYsum(alpha) # 1.683242467 4.107497821


### Exercise 6.23 (Fitting generalized hyperbolic distributions to equity return data)

## Data preparation
data(DJ_const)
tickers <- c("AAPL", "CSCO", "DIS", "IBM", "INTC", "MCD", "MSFT", "NKE", "PG", "WMT")
S <- DJ_const['2005-01-01/2012-12-31', tickers]
d <- ncol(S)
X.d <- returns(S) # build daily log-returns of the constituents


## a)

## Univariate fitting
fit.uv <- c(apply(X.d, 2, fit.tuv,    symmetric = TRUE,  silent = TRUE), # symmetric t
            apply(X.d, 2, fit.NIGuv,  symmetric = TRUE,  silent = TRUE), # symmetric NIG
            apply(X.d, 2, fit.hypuv,  symmetric = TRUE,  silent = TRUE), # symmetric H
            apply(X.d, 2, fit.ghypuv, symmetric = TRUE,  silent = TRUE), # symmetric GH
            apply(X.d, 2, fit.VGuv,   symmetric = TRUE,  silent = TRUE), # symmetric VG
            apply(X.d, 2, fit.tuv,    symmetric = FALSE, silent = TRUE), # asymmetric t
            apply(X.d, 2, fit.NIGuv,  symmetric = FALSE, silent = TRUE), # asymmetric NIG
            apply(X.d, 2, fit.hypuv,  symmetric = FALSE, silent = TRUE), # asymmetric H
            apply(X.d, 2, fit.ghypuv, symmetric = FALSE, silent = TRUE), # asymmetric GH
            apply(X.d, 2, fit.VGuv,   symmetric = FALSE, silent = TRUE)) # asymmetric VG

## Summary of results
res.uv <- matrix(sapply(fit.uv, function(x) x@aic), nrow = length(tickers), ncol = 10)
rownames(res.uv) <- tickers
model.names <- c("t.sym",  "NIG.sym",  "H.sym",  "GH.sym",  "VG.sym",
                 "t.asym", "NIG.asym", "H.asym", "GH.asym", "VG.asym")
colnames(res.uv) <- model.names
model.names[apply(res.uv, 1, which.min)] # for each stock, which is the model with smallest AIC?


## b)

## Multivariate fitting
fit.mv <- c(fit.tmv   (X.d, symmetric = TRUE,  silent = TRUE), # symmetric t
            fit.NIGmv (X.d, symmetric = TRUE,  silent = TRUE), # symmetric NIG
            fit.hypmv (X.d, symmetric = TRUE,  silent = TRUE), # symmetric H
            fit.ghypmv(X.d, symmetric = TRUE,  silent = TRUE), # symmetric GH
            fit.VGmv  (X.d, symmetric = TRUE,  silent = TRUE), # symmetric VG
            fit.tmv   (X.d, symmetric = FALSE, silent = TRUE), # asymmetric t
            fit.NIGmv (X.d, symmetric = FALSE, silent = TRUE), # asymmetric NIG
            fit.hypmv (X.d, symmetric = FALSE, silent = TRUE), # asymmetric H
            fit.ghypmv(X.d, symmetric = FALSE, silent = TRUE), # asymmetric GH
            fit.VGmv  (X.d, symmetric = FALSE, silent = TRUE)) # asymmetric VG
## Note:
## 1) Warnings are due to deprecated recycling of an array of length 1 ('ghyp' problem)
## 2) ...@gamma would be the fitted skewness parameter

## Summary of results
res.mv <- sapply(fit.mv, function(x) x@aic)
names(res.mv) <- model.names
which.min(res.mv) # => the symmetric t fits best according to AIC


## c)

## Extract fitted parameters
param <- coef(fit.mv[[which.min(res.mv)]], type = "chi.psi")
lambda <- param$lambda
chi <- param$chi
psi <- param$psi
mu <- param$mu
Sig <- param$sigma
gamma <- param$gamma

## We know that if X ~ GH_d(lambda, chi, psi, mu, Sigma, gamma) then
## b^T X ~ GH_1(lambda, chi, psi, b^T mu, b^T Sigma b, b^T gamma);
## see MFE (2015, Proposition 6.13). Let us define this implied distribution
## of daily portfolio returns.

## Parameters
b <- rep(1, d) # vector of weights
mu. <- drop(t(b) %*% mu)
sig. <- sqrt(drop(t(b) %*% Sig %*% b))
gamma. <- drop(t(b) %*% gamma)

## Implied distribution of daily portfolio returns
ret.distr <- ghyp(lambda = lambda, chi = chi, psi = psi,
                  mu = mu., sigma = sig., gamma = gamma.)

## Evaluate fitted GH density and a non-parametric density estimate (as a comparison)
ret <- drop(X.d %*% b) # returns
ran <- range(ret) # range where to plot
x <- seq(ran[1], ran[2], length.out = 201) # points where to evaluate fitted GH density
dens.fit <- dghyp(x, object = ret.distr) # evaluate fitted GH density there
dens.emp <- density(ret) # (empirical) density estimate

## Plot
ylim <- c(0, max(dens.emp$y, dens.fit))
plot(x, dens.fit, type = "l", ylim = ylim, col = "royalblue3",
     ylab = "Implied distribution of daily portfolio returns")
lines(dens.emp)
legend("topright", bty = "n", lty = c(1, 1), col = c("black", "royalblue3"),
       legend = c("Non-parametric estimate", "Fitted GH density"))


### Exercise 6.24 (Fitting a one-factor model to equity return data) ###########

## a)

## Data preparation
data(DJ_const)
tickers <- c("AAPL", "CSCO", "DIS", "IBM", "INTC", "MCD", "MSFT", "NKE", "PG", "WMT")
Sc <- DJ_const['2005-01-03/2012-12-31', tickers]
d <- ncol(Sc)
Xc <- returns(Sc) # build daily log-returns of the constituents

## Index
data(DJ) # index data
Si <- DJ['2005-01-01/2012-12-31']
Xi <- returns(Si) # build daily log-returns of the index

## Fit a multivariate macroeconomic regression model X_t = a + B * F_t + eps_t;
## see MFE (2015, 6.54), where B is the single-column matrix in our single-index
## model here (F_t is univariate).
X <- Xc # constituents' returns
F <- Xi # (single) index returns
(res <- lm(X ~ F)) # more details via summary()

## Get parameter estimates
par.ests <- coefficients(res)
a <- par.ests["(Intercept)",]
B <- par.ests["F",]

## High-beta stocks
sort(B, decreasing = TRUE) # => high-beta stocks (in this order)

## R^2 values
R2 <- apply(X, 2, function(x) cor(x, F)^2) # manually computed
R2. <- sapply(summary(res), function(res.) res.$r.squared) # extract R^2 values
names(R2.) <- tickers
stopifnot(all.equal(R2, R2.))
sort(R2, decreasing = TRUE) # => high R^2 stocks (in this order)x


## b)

## We investigate the residuals
eps <- resid(res) # (2012, 10)-matrix

## Correlation matrix
cor.eps <- cor(eps)
matrix_plot(cor.eps, ran = c(-1, 1)) # is it diagonal (as required)? => yes (roughly)

## Are the errors uncorrelated with the factors?
cor.eps.F <- cor(eps, F) # 10 correlations (idiosyncratic risk)
summary(cor.eps.F) # => yes (as required)

## Construct the implied covariance and correlation matrix
Ups <- cov(eps) # Upsilon (covariance matrix of epsilon)
Omega <- as.matrix(cov(F)) # Omega (covariance matrix of F; systematic risk)
Sigma <- B %*% Omega %*% t(B) + diag(diag(Ups)) # Cov(X); (10, 10)-matrix
rownames(Sigma) <- colnames(Sigma)
P <- cov2cor(Sigma) # Cor(X)

## Look at discrepancies between the factor model correlation matrix and the
## sample correlation matrix
err <- P-cor(X)
matrix_plot(err, ran = c(-2, 2)) # => fine

## The macroeconomic one-factor model is overall fine.


### Exercise 6.25 (Principal components analysis of equity return data) ########

## Data preparation
data(DJ_const)
tickers <- c("AAPL", "CSCO", "DIS", "IBM", "INTC", "MCD", "MSFT", "NKE", "PG", "WMT")
S <- DJ_const['2005-01-01/2012-12-31', tickers]
d <- ncol(S)
X.d <- returns(S) # build daily log-returns of the constituents


## a)

## Principal component analysis
PCA <- prcomp(X.d)
summary(PCA)
plot(PCA, xlab = "Principal component"); box() # show important components

## Extracting all information
Gamma <- PCA$rotation # principal axes (jth column is orthonormal eigenvector of cov(X) corresponding to jth largest eigenvalue) or 'loadings'
mu <- PCA$center # estimated centers
Y <- PCA$x # estimated principal components of X or 'scores'; (2012, 10)-matrix
var <- PCA$sdev^2 # explained variances per principal component
## Note: var equals the sorted eigenvalues of Cov(X) since diag(<sorted sigma^2>)
##       = Cov(Y) = Cov(Gamma^T (X - mu)) = Gamma^T Cov(X) Gamma = diag(<sorted lambda>)

## Sanity checks
stopifnot(max(abs(Gamma %*% t(Gamma) - diag(d))) < 1e-15, # Gamma Gamma^T = I_10 (orthogonality)
          max(abs(t(Gamma) %*% Gamma - diag(d))) < 1e-14, # Gamma^T Gamma = I_10 (orthogonality)
          apply(Gamma, 2, function(x) all.equal(sum(x * x), 1))) # orthonormality check
matrix_plot(cor(Y), ran = c(-1, 1)) # => as good as uncorrelated
Y. <- (X.d - rep(mu, each = nrow(X.d))) %*% Gamma # manually computing Y = Gamma^T (X - mu); see MFE (2015, (6.60))
stopifnot(all.equal(Y, Y., check.attributes = FALSE)) # compare Y and Y.
stopifnot(all(diff(var) < 0)) # check that variances are decreasing

## Proportion of variability explained by the first two and first three principal axes/loadings
prop <- cumsum(var)/sum(var)
prop[2]
prop[3]


## b)

## Loadings (= principal axes) of the first two principal axes (= first two eigenvectors)
Gamma[,1:2] # (could be plotted as component samples, so each column)

## Comment:
## The first vector of loadings has weights of the same sign (-ve) for
## all stocks. The negative signs are an artefact of the linear algebra and all
## the loading vectors could in principle be multiplied by -1. The key
## observation is that the signs are the same. It can be seen as a kind of index
## portfolio (if the weights in the loading vector are scaled to sum to 1 =>
## 'principal-component-mimicking portfolio').  The second vector has negative
## weights for all titles except for the Apple stock; as a portfolio it can be
## thought of as prescribing a programme of short selling of all titles to buy
## Apple, or vice versa.


## c)

npr <- 2 # number of principal components
Y. <- xts(Y[,seq_len(npr)], time(X.d)) # grab out the first so-many principal components of X
plot.zoo(Y., type = "h", xlab = "Time", ylab = paste("Component", seq_len(npr)),
         main = "Principal components of X") # plot of the first so-many principal components of X
## Using the first two principal components as factors, the corresponding
## model is X_t = mu + B F_t + eps_t, where B = Gamma_1 and F_t = Y[,1:2].
