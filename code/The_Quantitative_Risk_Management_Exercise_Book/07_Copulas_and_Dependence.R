## By Marius Hofert

## R script for Chapter 7 of The QRM Exercise Book


### Setup ######################################################################

library(copula)
library(xts)
library(qrmtools)
library(qrmdata)
doPDF <- require(crop)
options(digits = 10)


### Exercise 7.8 (Distinguishing meta-distributions from scatterplots) #########

## Reproducing code

## Setup
n <- 1000 # sample size
d <- 2 # dimension
tau <- 0.5 # Kendall's tau

## a) Normal copula with N(0,1) margins, so N(0, P) sample
family <- "normal"
rho <- iTau(ellipCopula(family), tau = tau)
cop <- ellipCopula(family, param = rho, dim = d)
set.seed(271)
U.N <- rCopula(n, copula = cop)
X.N <- qnorm(U.N)

## b) t_nu copula with N(0,1) margins
##    Note: same 'rho' as above
family <- "t"
nu <- 3.5
cop <- ellipCopula(family, param = rho, dim = d, df = nu)
set.seed(271)
U.t <- rCopula(n, copula = cop)
X.t <- qnorm(U.t)

## c) t_nu(0, P) copula with negative P_{12} and N(0,1) margins
cop <- ellipCopula(family, param = -rho, dim = d, df = nu)
set.seed(271)
U.t.m <- rCopula(n, copula = cop)
X.t.m <- qnorm(U.t.m)

## d) t_nu(0, I) copula with N(0,1) margins
cop <- ellipCopula(family, param = 0, dim = d, df = nu)
set.seed(271)
U.t.I <- rCopula(n, copula = cop)
X.t.I <- qnorm(U.t.I)

## e) Gumbel copula with N(0,1) margins
family <- "Gumbel"
th <- iTau(archmCopula(family), tau)
cop <- archmCopula(family, param = th, dim = d)
set.seed(271)
U.G <- rCopula(n, copula = cop)
X.G <- qnorm(U.G)

## f) Clayton copula with N(0,1) margins
family <- "Clayton"
th <- iTau(archmCopula(family), tau)
cop <- archmCopula(family, param = th, dim = d)
set.seed(271)
U.C <- rCopula(n, copula = cop)
X.C <- qnorm(U.C)

## g) Copula sample of an elliptical distribution with bounded radial part (with N(0,1) margins)
set.seed(271)
Z <- matrix(rnorm(n * d), ncol = d)
U <- Z / sqrt(rowSums(Z^2))
P <- p2P(rho)
A <- t(chol(P))
R <- 1 + runif(n) # R ~ U(1,2)
U.bdd.R <- pobs(R * t(A %*% t(U)))
X.bdd.R <- qnorm(U.bdd.R)

## h) Marshall--Olkin copula sample with N(0,1) margins
alpha <- c(0.2, 0.6)
set.seed(271)
V <- matrix(runif(n * 3), ncol = 3)
U.MO <- cbind(pmax(V[,1]^(1/(1 - alpha[1])), V[,3]^(1/alpha[1])),
              pmax(V[,2]^(1/(1 - alpha[2])), V[,3]^(1/alpha[2])))
X.MO <- qnorm(U.MO)

## i) Copula sample of a normal mean-variance mixture with N(0,1) margins
set.seed(271)
W <- rexp(n) # too extreme: nu / rchisq(n, df = nu)
mu <- c(1, 1) # need a mu here (mu = 0 => not a normal mean-variance mixture)
X <- rep(mu, each = n) * W + sqrt(W) * t(A %*% t(Z)) # of the form m(W) + sqrt(W) * A * Z
U.nmvmix <- pobs(X)
X.nmvmix <- qnorm(U.nmvmix)

## Plot
r <- range(X.N, X.t, X.t.m, X.t.I, X.G, X.C, X.bdd.R, X.MO, X.nmvmix)
rmx <- max(abs(r))
ran <- c(-rmx, rmx) # symmetricize
if(doPDF) pdf(file = (file <- "fig_07_recognizing_features.pdf"), width = 9, height = 9)
opar <- par(mar = c(5, 4, 4, 2) + 0.1 - c(1, 0, 3, 1)) # # set margins (bottom, left, top, right)
lay <- matrix(1:9, ncol = 3, byrow = TRUE) # layout matrix
layout(lay) # layout
## Note: Order obtained from set.seed(314); c("N", "t", "t.m", "t.I", "G", "C", "bdd.R", "MO", "nmvmix")[sample(1:9)]
plot(X.N,      cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
plot(X.t.m,    cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
plot(X.C,      cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
plot(X.t,      cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
plot(X.bdd.R,  cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
plot(X.G,      cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
plot(X.nmvmix, cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
plot(X.MO,     cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
plot(X.t.I,    cex = 0.7, xlab = expression(X[1]), ylab = expression(X[2]), xlim = ran, ylim = ran)
par(opar)
if(doPDF) dev.off.crop(file)


### Exercise 7.12 (Minimum correlation for exponential random variables) #######

## b)
integrate(function(u) log(1-u)*log(u), lower = 0, upper = 1)$value - 1 # minimal correlation


### Exercise 7.13 (Maximal correlation for Pareto Type I random variables) #####

## Reproducing code

##' @title Correlation bounds for Par(theta_.) margins with theta_. > 2
##' @param th 2-column matrix or 2-vector of sigma's
##' @param method character string indicating the minimal or maximal correlation
##'        bound to be plotted
##' @return correlation bound
##' @author Marius Hofert
cor_bound_ParT1 <- function(th, method = c("max", "min"), N = 1e4)
{
    ## th = (theta_1, theta_2)
    if(!is.matrix(th)) th <- rbind(th, deparse.level = 0)
    E1 <- 1/(1-1/th[,1]) # E(X_1)
    E2 <- 1/(1-1/th[,2]) # E(X_1)
    E1.2 <- 1/(1-2/th[,1]) # E(X_1^2)
    E2.2 <- 1/(1-2/th[,2]) # E(X_2^2)
    method <- match.arg(method)
    E12 <- switch(method,
    "min" = {
        stopifnot(N >= 0)
        k <- 0:N
        apply(th, 1, function(th.) sum(choose(-1/th.[1], k) * (-1)^k / (k-1/th.[2]+1)))
    },
    "max" = {
        1/(1-1/th[,1]-1/th[,2])
    },
    stop("Wrong 'method'"))
    (E12 - E1 * E2) / sqrt((E1.2-E1^2) * (E2.2-E2^2))
}

## Extreme correlations for Pareto Type 1 margins
n.grid <- 26 # number of grid points in each dimension
th <- seq(2.01, 5, length.out = n.grid) # subdivision points in each dimension
grid <- expand.grid("theta[1]" = th, "theta[2]" = th) # build a grid
## Maximal
val.max <- cbind(grid, "rho[max](theta[1],theta[2])" =
                 cor_bound_ParT1(grid))
if(doPDF) pdf(file = (file <- "fig_07_cor_bounds_ParT1_up.pdf"))
wireframe2(val.max) # maximal correlation bound
if(doPDF) dev.off.crop(file)
## Minimal
val.min <- cbind(grid, "rho[min](theta[1],theta[2])" =
                 cor_bound_ParT1(grid, method = "min"))
if(doPDF) pdf(file = (file <- "fig_07_cor_bounds_ParT1_low.pdf"))
wireframe2(val.min) # minimal correlation bound (also included here for comparison)
if(doPDF) dev.off.crop(file)


### Exercise 7.24 (Sampling copulas and meta-distributions) ####################

## a)

## Define parameters
th.N <- 0.7
th.G <- 2
th.C <- 2.2
th.t <- 0.71

## Define copulas
Ncop <- normalCopula(th.N)
Gcop <- gumbelCopula(th.G)
Ccop <- claytonCopula(th.C)
tcop <- tCopula(th.t) # note: df = 4 is the default (and thus omitted)

## Generate samples
n <- 2000
set.seed(271) # for reproducibility
U.N <- rCopula(n, copula = Ncop)
U.G <- rCopula(n, copula = Gcop)
U.C <- rCopula(n, copula = Ccop)
U.t <- rCopula(n, copula = tcop)

##' @title Plot Copula Samples One by One
##' @param ... named copula samples
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ran x- and y-range
##' @param cex character expansion factor
##' @return invisible()
##' @author Marius Hofert
splot <- function(..., xlab, ylab, ran = NULL, cex = 0.7)
{
    opar <- par(pty = "s", ask = interactive())
    arg <- list(...)
    if(is.null(ran))
        ran <- do.call(range, arg)
    for(i in seq_len(length(arg))) {
        plot(arg[[i]], xlim = ran, ylim = ran, cex = cex,
             xlab = xlab, ylab = ylab)
        mtext(substitute("Sample of size"~n==n.~~"from a"~nm,
                         list(n. = nrow(arg[[i]]),
                              nm = parse(text = names(arg)[i])[[1]])),
              side = 4, line = 1, adj = 0)
    }
    par(opar)
    invisible()
}

## Plots
splot("normal~copula" = U.N, "Gumbel~copula" = U.G,
      "Clayton~copula" = U.C, "t[4]~copula" = U.t,
      xlab = expression(U[1]), ylab = expression(U[2]))


## b)

## Generate data
## Note: We recycle the copula samples, just map to N(0,1) margins.
X.N <- qnorm(U.N)
X.G <- qnorm(U.G)
X.C <- qnorm(U.C)
X.t <- qnorm(U.t)

## Plots
splot("meta-normal~distribution" = X.N, "meta-Gumbel~distribution" = X.G,
      "meta-Clayton~distribution" = X.C, "meta-t[4]~distribution" = X.t,
      xlab = expression(X[1]), ylab = expression(X[2]))


## c)

## Define parameters
(tau <- tau(Ncop)) # compute Kendall's tau for the normal copula (~= 0.5)
th.G. <- iTau(gumbelCopula(), tau = tau) # compute parameter of Gumbel copula such that Kendall's tau matches 'tau'
th.C. <- iTau(claytonCopula(), tau = tau)
th.t. <- iTau(tCopula(), tau = tau)

## Define (new) Gumbel, Clayton and t copulas
Gcop. <- gumbelCopula(th.G.)
Ccop. <- claytonCopula(th.C.)
tcop. <- tCopula(th.t.)

## Generate samples
n <- 2000
set.seed(271) # for reproducibility
X.N. <- qexp(U.N) # here we can recycle the copula sample from before
X.G. <- qexp(rCopula(n, copula = Gcop.))
X.C. <- qexp(rCopula(n, copula = Ccop.))
X.t. <- qexp(rCopula(n, copula = tcop.))

## Plots
splot("meta-normal~distribution" = X.N., "meta-Gumbel~distribution" = X.G.,
      "meta-Clayton~distribution" = X.C., "meta-t[4]~distribution" = X.t.,
      xlab = expression(X[1]), ylab = expression(X[2]))


### Exercise 7.25 (Fitting copulas to equity return data) ######################

## Data preparation
data(DJ_const)
tickers <- c("AAPL", "CSCO", "DIS", "IBM", "INTC", "MCD", "MSFT", "NKE", "PG", "WMT")
S <- DJ_const['2005-01-01/2012-12-31', tickers]
d <- ncol(S)


## a)

## Build daily log-returns and their pseudo-observations
X.d <- returns(S)
U.d <- pobs(X.d)
## Note:
## 1) We normally work with negative log-returns but it doesn't matter here as
##    we only fit radially symmetric copula models.
## 2) A more in-depth modeling procedure would need to take marginal non-stationarity
##    into account (for example, by fitting ARMA-GARCH models and then building
##    the pseudo-observations of the standardized residuals).

## Fitting
## Normal copula
Ncop <- normalCopula(dim = d, dispstr = "un")
system.time(fit.N.d <- fitCopula(Ncop, data = as.matrix(U.d))) # uses the default method "mpl" (for maximum pseudo-likelihood based estimation); ~= 35s
## t copula
tcop <- tCopula(dim = d, dispstr = "un")
system.time(fit.t.d <- fitCopula(tcop, data = as.matrix(U.d))) # ~= 3min

## Fitted parameters
p2P(coef(fit.N.d)) # correlation matrix of the fitted normal copula
p2P(tail(coef(fit.t.d), n = -1)) # correlation matrix of the fitted t copula
tail(coef(fit.t.d), n = 1) # degrees of freedom of the fitted t copula


## b)

## Comparison of AICs
k <- choose(d, 2) # number of estimated parameters (also for t as d.o.f. is treated as 'fixed')
AIC.N.d <- 2*k - 2*fit.N.d@loglik # AIC for normal copula
AIC.t.d <- 2*k - 2*fit.t.d@loglik # AIC for t copula
(AIC.N.d - AIC.t.d) / abs(AIC.N.d) # (> 0 => AIC of t is smaller)
## => AIC of the t copula is about 10% smaller


## c)

## Build monthly log-returns and their pseudo-observations
X.m <- apply.monthly(X.d, FUN = colSums)
U.m <- pobs(X.m)

## Fitting (faster here as the sample size is smaller)
## Normal copula
system.time(fit.N.m <- fitCopula(Ncop, data = as.matrix(U.m))) # ~= 3s
## t copula
system.time(fit.t.m <- fitCopula(tcop, data = as.matrix(U.m))) # ~= 10s

## Fitted parameters
p2P(coef(fit.N.m)) # correlation matrix of the fitted normal copula
p2P(tail(coef(fit.t.m), n = -1)) # correlation matrix of the fitted t copula
tail(coef(fit.t.m), n = 1) # degrees of freedom of the fitted t copula
## => quite a bit larger (data more 'normal' due to Central Limit Theorem effect)

## Comparison of AICs
AIC.N.m <- 2*k - 2*fit.N.m@loglik # AIC for normal copula
AIC.t.m <- 2*k - 2*fit.t.m@loglik # AIC for t copula => smaller
(AIC.N.m - AIC.t.m) / abs(AIC.N.m)
## => AIC of the t copula is still smaller by 4% but not as much anymore as
##    the monthly data is less non-normal (Central Limit Theorem effect)
