## By Marius Hofert

## Copula estimation and goodness-of-fit (5d example, 2 copulas)


### 0 Setup ####################################################################

library(copula)
library(qrmtools)
library(qrmdata)
library(rugarch)
library(xts)
library(ADGofTest)
library(qqtest)


### 1 Working with the data ####################################################

## Select the data we work with
data("SP500_const") # load the constituents data of the S&P 500
stocks <- c("INTC", "QCOM", "GOOGL", "AAPL", "MSFT") # Intel, Qualcomm, Google, Apple, Microsoft
time <- c("2007-01-03", "2009-12-31") # time period
S <- SP500_const[paste0(time, collapse = "/"), stocks] # data

## Check for missing data
stopifnot(all(!is.na(S))) # na.fill(, fill = "extend") is often helpful


### 2 Fitting marginal ARMA(1,1)-GARCH(1,1) models with standardized t residuals

## Build negative log-returns
X <- -log_returns(S) # -log-returns

## Basic plot
plot.zoo(X, main = "Negative log-returns")

## Fit marginal time series
uspec <- rep(list(ugarchspec(distribution.model = "std")), ncol(X))
fit.ARMA.GARCH <- fit_ARMA_GARCH(X, ugarchspec.list = uspec)
stopifnot(sapply(fit.ARMA.GARCH$error, is.null)) # NULL = no error
if(FALSE)
    fit.ARMA.GARCH$warning
    ## => Warning comes from finding initial values and can be ignored here
fits <- fit.ARMA.GARCH$fit # fitted models
Z <- as.matrix(do.call(merge, lapply(fits, residuals, standardize = TRUE))) # grab out standardized residuals
colnames(Z) <- colnames(S)
(nu.mar <- vapply(fits, function(x) x@fit$coef[["shape"]], NA_real_)) # vector of estimated df
n <- nrow(X) # sample size
d <- ncol(X) # dimension


### 3 Fitting copulas ##########################################################

## Compute pseudo-observations from the standardized t residuals
U <- pobs(Z)
pairs2(U, cex = 0.4, col = adjustcolor("black", alpha.f = 0.5))

## Fitting a Gumbel copula
fit.gc <- fitCopula(gumbelCopula(dim = d),
                    data = U, method = "mpl")
fit.gc@estimate # estimated copula parameter
gc <- fit.gc@copula # fitted copula

## Compute matrices of pairwise Kendall's tau and upper tail-dependence coefficients
p2P(tau(gc), d = d)
p2P(lambda(gc)["upper"], d = d)

## Fitting a t copula
fit.tc <- fitCopula(tCopula(dim = d, dispstr = "un"),
                    data = U, method = "itau.mpl")
(nu <- tail(fit.tc@estimate, n = 1)) # estimated degrees of freedom nu
(P <- p2P(head(fit.tc@estimate, n = -1))) # estimated correlation matrix
tc <- fit.tc@copula # fitted copula

## Compute matrices of pairwise Kendall's tau and upper tail-dependence coefficients
p2P(tau(tc))
p2P(lambda(tc)[(choose(d,2)+1):(d*(d-1))])


### 4 Goodness-of-fit ##########################################################

## We use the parametric bootstrap here, based on the maximum pseudo-likelihood
## estimator.
set.seed(271) # for reproducibility
N <- 100 # this is to save run time; it should be larger!

## Check the Gumbel copula
gof.gc <- gofCopula(gc, x = U, N = N)
stopifnot(gof.gc$p.value < 0.05) # => rejection

## Check the t copula
## Note: - This can currently only be done for fixed and integer degrees of
##         freedom as there is no algorithm to evaluate the multivariate t df for
##         non-integer degrees of freedom.
##       - ... yet it's still quite slow. We thus check the model graphically
##         after mapping the variates to a U(0,1)^d distribution with the
##         Rosenblatt transform.
U.Rsnbl <- cCopula(U, copula = tc)
pairs2(U.Rsnbl, cex = 0.4, col = adjustcolor("black", alpha.f = 0.5)) # looks ok

## Map it to a K_d distribution and do a Q-Q plot
U.Rsnbl.K <- sqrt(rowMeans(qnorm(U.Rsnbl)^2)) # map to a K_d
pK <- function(q, d) pchisq(d*q*q, df = d) # df of a K_d distribution
AD <- ad.test(U.Rsnbl.K, distr.fun = pK, d = d) # compute an AD test
stopifnot(AD$p.value >= 0.05)
## Note: The AD test here does not take into account parameter estimation
##       (that would require a parametric bootstrap, for example)

## A (sophisticated) Q-Q plot
qqtest(U.Rsnbl.K, dist = "kay", df = d, nreps = 1000, pch = 1,
       col = adjustcolor("black", alpha.f = 0.5), main = "",
       xlab = substitute(K[dof]~"quantiles", list(dof = d))) # => looks ok


### 5 Simulating paths from the full model #####################################

## Simulate paths
B <- 200
m <- ceiling(n/10) # length of the simulates paths
X.lst <- lapply(1:B, function(b) {
    ## 1) Simulate from the fitted copula
    U. <- rCopula(m, copula = tc)
    ## 2) Quantile-transform to standardized t distributions (for ugarchsim())
    Z. <- sapply(1:d, function(j) sqrt((nu.mar[j]-2)/nu.mar[j]) * qt(U.[,j], df = nu.mar[j]))
    ## 3) Use these multivariate dependent t innovations to sample from the
    ##    time series
    sim <- lapply(1:d, function(j)
        ugarchsim(fits[[j]], n.sim = m, m.sim = 1,
                  custom.dist = list(name = "sample",
                                     distfit = Z.[,j, drop = FALSE])))
    ## 4) Extract simulated series
    sapply(sim, function(x) fitted(x)) # simulated multivariate series X_t (= x@simulation$seriesSim)
})
## => List of length B containing (n x d)-matrices


### 6 Predict the aggregated loss and VaR_0.99 #################################

## Note: - This is merely a demo of what can be done with the simulated data.
##       - See also the vignette "ARMA_GARCH_VaR" in qrmtools

## Predict the aggregated loss and VaR (nonparametrically)
Xs <- rowSums(X) # aggregated loss; n-vector
Xs. <- sapply(X.lst, rowSums) # simulated aggregated losses; (m, B)-matrix
Xs.mean <- rowMeans(Xs.) # predicted aggregated loss; m-vector
Xs.CI <- apply(Xs., 1, function(x) quantile(x, probs = c(0.025, 0.975))) # CIs; (2, m)-matrix
alpha <- 0.99 # confidence level
VaR <- apply(Xs., 1, function(x) quantile(x, probs = alpha)) # VaR_alpha; m-vector

## Plot
tm <- index(SP500_const)
start <- match(time[1], as.character(tm))
past <- tm[start:(start+n-1)]
future <- tm[(start+n):(start+n+m-1)]
plot(past, Xs, type = "l", xlim = range(c(past, future)), xlab = "", ylab = "") # actual (past) losses
polygon(c(future, rev(future)), c(Xs.CI[1,], rev(Xs.CI[2,])),
        border = NA, col = "grey80") # CI region
lines(future, Xs.mean, col = "royalblue3") # predicted aggregated loss
lines(future, Xs.CI[1,], col = "grey50") # lower CI
lines(future, Xs.CI[2,], col = "grey50") # upper CI
lines(future, VaR, col = "maroon3") # VaR_alpha
legend("bottomright", bty = "n", lty = rep(1, 4),
       col = c("black", "royalblue3", "grey50", "maroon3"),
       legend = c("(Aggregated) loss", "(Simulated) predicted loss",
                  "95% CIs", as.expression(substitute("Simulated"~VaR[a], list(a = alpha)))))


