## By Alexander McNeil and Marius Hofert

## Confidence intervals for the GEV parameters in 05_GEV_BMM_SP500.R
## based on the profile likelihood method


### Background information (see Davison (2003, p. 126)) ########################

## 1) Log-likelihood theory:
##    - Likelihood ratio statistic:
##
##          W(th) = 2(LL(MLE) - LL(th))
##
##    - True underlying unknown parameter: th_0 (p-dimensional)
##    - Under regularity conditions,
##
##          W(th_0) -> chi_p^2 (in distribution)
##
##    - An asymptotic (1-alpha)-CI for th_0 is thus:
##
##          {th :  W(th) <= (chi_p^2)^-(1-alpha)}
##        = {th : LL(th) >= LL(MLE) - (chi_p^2)^-(1-alpha) / 2}
##
##    - Replacing ">=" by "=", we have to find the two roots
##
##          LL(th) - (LL(MLE) - (chi_p^2)^-(1-alpha) / 2) = 0
##
##      (corresponding to the lower and upper CI endpoints).
##
## 2) Idea profile log-likelihoods:
##    - Assume only (CIs for) some parameters th.i are of interest
##    - th.i = 'parameters of interest';
##      th.n = 'nuisance parameters' (remaining ones)
##    - One can define the 'generalized likelihood ratio statistic'
##
##          W(th.i) = 2(LL(MLE) - LL(th.i, hat(th.n)))
##
##      where hat(th.n) is the maximizer of the 'profile log-likelihood'
##
##          pLL(th.i) := LL(th.i, th.n) in th.n for given th.i.
##
##    - Under regularity conditions,
##
##          W(th.i) -> chi_{p.i}^2 (in distribution)
##
##      where p.i = length(th_i).
##
##    - An asymptotic (1-alpha)-CI for the true underlying, unknown th_{i,0} is
##      thus
##
##          {th_i :             W(th_i) <= (chi_{p.i}^2)^-(1-alpha)}
##        = {th_i : LL(th.i, hat(th.n)) >= LL(MLE) - (chi_{p.i}^2)^-(1-alpha) / 2}.
##
##    - Replacing ">=" by "=", we have to find the two roots
##
##          LL(th.i, hat(th.n)) - (LL(MLE) - (chi_{p.i}^2)^-(1-alpha) / 2) = 0
##
##      (corresponding to the lower and upper CI endpoints).


### Setup ######################################################################

library(xts) # for functions around time series objects
library(qrmdata) # for the S&P 500 data
library(qrmtools) # for returns()


##' @title Profile Log-Likelihood
##' @param r return level (parameter of interest)
##' @param k return period (parameter of interest)
##' @param x maxima
##' @param control see ?optim
##' @param ... additional arguments passed to optim()
##' @return return value of optim() with the minimizer of the -log-likelihood
##'         in the nuisance parameters theta = (xi, sigma) for fixed r or k.
##'         In other words, profLogLik()$value gives the profile log-likelihood
##'         in r or k.
##' @author Alexander McNeil and Marius Hofert
##' @note 1) r = H_{xi,mu,sig}^-(1-1/k) => mu implied by xi, sig, r, k
##'       2) reltol = 0 improves the smoothness of the objective function below
profLogLik <- function(r, k, x, control = list(reltol = 0), ...)
{
    ## log-likelihood as a function of the nuisance parameters xi, sigma
    ## (mu is implied from the given r, k)
    ## Note: Could very well be -Inf due to choice of r or k and thus implied mu
    logLik.xi.sig <- function(theta) { # note: r, k and x are locally 'seen' here
        xi <- theta[1]
        sig <- theta[2]
        imu <- if(xi == 0) { # mu implied by xi, sigma, r, k
                   r + sig * log(-log1p(-1/k))
               } else {
                   r - (sig/xi) * ((-log1p(-1/k))^(-xi) - 1)
               }
        stopifnot(is.finite(imu)) # sanity check
        sum(dGEV(x, xi = xi, mu = imu, sigma = sig, log = TRUE)) # log-likelihood (correctly deals with sigma <= 0)
    }
    ## Construct an initial value which guarantees that logLik.xi.sig is finite
    ## Note:
    ## 1) The case 'xi = 0' guarantees that 1 + xi (x - mu) / sigma > 0 but
    ##    since we use the implied mu here (to force a given r or k), y = (x - mu) / sigma
    ##    can still be so negative (e.g. ~= -1e4) such that exp(-y) = Inf and thus
    ##    dGEV(, log = TRUE) = -Inf; see code of dGEV(). We thus adjust sigma
    ##    here to guarantee that we get a finite initial log-likelihood.
    ## 2) y = (x - mu) / sigma should be >= -log(.Machine$double.xmax) for all y
    ##    so that exp(-y) < Inf for all y. So sigma >= (x - mu) / -log(.Machine$double.xmax)
    ##    for all x (so min(x) [since -log(..) < 0]).
    ## 3) ... but we can't take .Machine$double.xmax as sum() (over sample size)
    ##    would still be -Inf.
    ##    Version 1: take .Machine$double.xmax / sample size
    ##    imu <- r + init[2] * log(-log1p(-1/k)) # implied mu
    ##    const <- -log(.Machine$double.xmax/1e8) # good up to sample size 1e8
    ##    xmin <- min(x)
    ##    if((xmin - imu) / init[2] < const)
    ##         init[2] <- (xmin - imu) / const
    ##    Version 2: simply double sigma
    init <- c(0, sqrt(6 * var(x)) / pi) # case 'xi = 0' (xi, sigma)
    while(!is.finite(logLik.xi.sig(init))) init[2] <- init[2] * 2
    ## Maximize the log-likelihood in the nuisance parameters theta = (xi, sigma)
    ## for our given parameters of interest r and k (and the data x)
    control <- c(as.list(control), fnscale = -1) # maximization (overwrites possible additionally passed 'fnscale')
    optim(init, fn = logLik.xi.sig, control = control, ...) # maximize logLik.xi.sig
}

##' @title Objective Function for Finding (1-alpha)-Confidence Intervals for the
##'        Return Level r or the Return Period k
##' @param r return level
##' @param k return period
##' @param x maxima
##' @param maxLogLik maximal (overall) log-likelihood (when maximizing GEV in all parameters)
##' @param alpha significance level
##' @param control see ?optim
##' @param ... additional arguments passed to profLogLik()
##' @return One side of the confidence level (which one is determined by
##'         the initial values)
##' @author Alexander McNeil and Marius Hofert
##' @note - Our parameter of interest is either r or k (fix the other).
##'       - The nuisance parameters are xi and sigma; mu is implied when
##'         fixing r, k, xi, sigma
obj <- function(r, k, x, maxLogLik, alpha = 0.05, control = list(reltol = 0), ...)
    (profLogLik(r = r, k = k, x = x, control = control, ...)$value) - # profile log-likelihood
        (maxLogLik - qchisq(1-alpha, df = 1) / 2)

##' @title Plot obj() as a Function of r or k
##' @param r return level (if of length > 1, obj() is evaluated as a function of r)
##' @param k return perid (if of length > 1, obj() is evaluated as a function of k)
##' @param x maxima
##' @param maxLogLik maximal (overall) log-likelihood (when maximizing GEV in all parameters)
##' @param alpha significance level
##' @param control see ?optim
##' @param ... additional arguments passed to profLogLik()
##' @return invisible(); plot as side-effect
##' @author Marius Hofert
plot_obj <- function(r, k, x, maxLogLik, alpha = 0.05, control = list(reltol = 0), ...)
{
    lr <- length(r)
    lk <- length(k)
    stopifnot(xor(lr == 1, lk == 1)) # one of them has to have length > 1 (otherwise no plot)
    if(lr > 1) {
        x. <- r
        y. <- sapply(r, function(r.) obj(r = r., k = k, x = x,
                                         maxLogLik = maxLogLik, alpha = alpha,
                                         control = control, ...))
        xlab <- "Return level r"

    } else { # length(k) > 1
        x. <- k
        y. <- sapply(k, function(k.) obj(r = r, k = k., x = x,
                                         maxLogLik = maxLogLik, alpha = alpha,
                                         control = control, ...))
        xlab <- "Return period k"
    }
    plot(x., y., type = "l", xlab = xlab,
         ylab = "Objective function for finding root")
    abline(h = 0, lty = 2)
}


### 1 Working with the data ####################################################

## Load the data and compute the negative log-returns (risk-factor changes X)
data(SP500)
S <- SP500 # 'xts'/'zoo' object
X <- -returns(S) # -log-returns X_t = -log(S_t/S_{t-1})
X. <- X['1960-01-01/1987-10-16'] # grab out data we work with

## Extract (half-)yearly maxima method
M.y <- period.apply(X., INDEX = endpoints(X., "years"), FUN = max) # yearly maxima
endpts <- endpoints(X., "quarters") # end indices for quarters
endpts <- endpts[seq(1, length(endpts), by = 2)] # end indices for half-years
M.hy <- period.apply(X., INDEX = endpts, FUN = max) # half-yearly maxima


### 3 Profile-likelihood-based CIs for return levels and return periods ########

### 3.1 Based on annual maxima #################################################

## Fit GEV to yearly maxima
fit.y <- fit_GEV(M.y) # likelihood-based estimation of the GEV; see ?fit_GEV
(xi.y <- fit.y$par[1]) # ~= 0.2972 => Frechet domain with infinite ceiling(1/xi.y) = 4th moment
(mu.y  <- fit.y$par[2])
(sig.y <- fit.y$par[3])
sqrt(diag(fit.y$Cov)) # standard errors
mLL <- fit.y$value # ~= 88.5288; maximum log-likelihood (at MLE)


### 3.1.1 CI for 10-year and 50-year return level ##############################

## Fix the return period k = 10 (in years) and estimate the corresponding
## return level r
k <- 10 # fix return period
r <- qGEV(1-1/k, xi = xi.y, mu = mu.y, sigma = sig.y) # corresponding estimated return level r
## Find a suitable initial interval for computing the CI
I <- c(0.02, 0.1) # found by experimenting with the following plot
plot_obj(r = seq(I[1], I[2], length.out = 128), k = k, x = M.y, maxLogLik = mLL)
abline(v = r) # estimated return level
## Compute 95%-CI for the 10-year return level r
CI.low <- uniroot(obj, lower = I[1], upper = r,    k = k, x = M.y, maxLogLik = mLL)
CI.up  <- uniroot(obj, lower = r,    upper = I[2], k = k, x = M.y, maxLogLik = mLL)
(CI.r10 <- c(CI.low$root, CI.up$root)) # 95%-CI for 10-year return level r; [0.0346, 0.0757]
## Does it contain a drop as large as on Black Monday?
rBM <- as.numeric(X['1987-10-19']) # return level on Black Monday
CI.r10[1] <= rBM && rBM <= CI.r10[2] # => no!

## The same for k = 50
k <- 50 # fix return period
r <- qGEV(1-1/k, xi = xi.y, mu = mu.y, sigma = sig.y) # corresponding estimated return level r
## Find a suitable initial interval for computing the CI
I <- c(0.02, 0.3) # found by experimenting with the following plot
plot_obj(r = seq(I[1], I[2], length.out = 128), k = k, x = M.y, maxLogLik = mLL)
abline(v = r) # estimated return level
## Compute 95%-CI for the 50-year return level r
CI.low <- uniroot(obj, lower = I[1], upper = r,    k = k, x = M.y, maxLogLik = mLL)
CI.up  <- uniroot(obj, lower = r,    upper = I[2], k = k, x = M.y, maxLogLik = mLL)
(CI.r50 <- c(CI.low$root, CI.up$root)) # 95%-CI for 50-year return level r; [0.0488, 0.2483]
## Does it contain a drop as large as on Black Monday?
rBM <- as.numeric(X['1987-10-19']) # return level on Black Monday
CI.r50[1] <= rBM && rBM <= CI.r50[2] # => yes!

## Note: If obj() is called *without* smaller relative tolerance
##       (e.g., the default control = list()), then the upper CI differs!
uniroot(obj, lower = r, upper = I[2], k = k, x = M.y, maxLogLik = mLL,
        control = list())$root
## This would imply that the return level on Black Monday is *not* contained in the CI!
## Why does this happen? Watch this...
plot_obj(r = seq(I[1], I[2], length.out = 128), k = k, x = M.y, maxLogLik = mLL,
         control = list())
## => Always watch out for numerical issues!


### 3.1.2 CI for the return period of a loss as on Black Monday ################

## Fix the return level r (as on Black Monday) and estimate the corresponding
## return period k
r <- rBM # fix return level
k <- 1/(1-pGEV(r, xi = xi.y, mu = mu.y, sigma = sig.y)) # corresponding estimated return period k
## Find a suitable initial interval for computing the CI
I <- c(10, 1e5) # found by experimenting with the following plot
plot_obj(r = r, k = seq(I[1], I[2], length.out = 128), x = M.y, maxLogLik = mLL)
abline(v = k) # estimated return level
## => As good as impossible to find the upper CI (even for larger I[2])
##    Most likely requires to rescale the objective function etc.
## Compute (at least the lower end of) the 95%-CI for the return period (in years)
## of an event of size as Black Monday
CI.low <- uniroot(obj, lower = I[1], upper = k, r = r, x = M.y, maxLogLik = mLL)
CI.low$root


### 3.2 Based on biannual maxima ###############################################

## Fit GEV to half-yearly maxima
fit.hy <- fit_GEV(M.hy) # likelihood-based estimation of the GEV; see ?fit_GEV
(xi.hy <- fit.hy$par[1]) # ~= 0.3401 => Frechet domain with infinite ceiling(1/xi.hy) = 3rd moment
(mu.hy  <- fit.hy$par[2])
(sig.hy <- fit.hy$par[3])
sqrt(diag(fit.hy$Cov)) # standard errors
mLL <- fit.hy$value # ~= 191.3122; maximum log-likelihood (at MLE)


### 3.2.1 CI for 20-half-year and 100-half-year return level ###################

## Fix the return period k = 20 (in half-years) and estimate the corresponding
## return level r
k <- 20 # fix return period
r <- qGEV(1-1/k, xi = xi.hy, mu = mu.hy, sigma = sig.hy) # corresponding estimated return level r
## Find a suitable initial interval for computing the CI
I <- c(0.03, 0.1) # found by experimenting with the following plot
plot_obj(r = seq(I[1], I[2], length.out = 128), k = k, x = M.hy, maxLogLik = mLL)
abline(v = r) # estimated return level
## Compute 95%-CI for the 20-half-year return level r
CI.low <- uniroot(obj, lower = I[1], upper = r,    k = k, x = M.hy, maxLogLik = mLL)
CI.up  <- uniroot(obj, lower = r,    upper = I[2], k = k, x = M.hy, maxLogLik = mLL)
(CI.r20 <- c(CI.low$root, CI.up$root)) # 95%-CI for 20-half-year return level r; [0.0355, 0.0728]
## Does it contain a drop as large as on Black Monday?
rBM <- as.numeric(X['1987-10-19']) # return level on Black Monday
CI.r20[1] <= rBM && rBM <= CI.r20[2] # => no!

## The same for k = 100
k <- 100 # fix return period
r <- qGEV(1-1/k, xi = xi.hy, mu = mu.hy, sigma = sig.hy) # corresponding estimated return level r
## Find a suitable initial interval for computing the CI
I <- c(0.04, 0.25) # found by experimenting with the following plot
plot_obj(r = seq(I[1], I[2], length.out = 128), k = k, x = M.hy, maxLogLik = mLL)
abline(v = r) # estimated return level
## Compute 95%-CI for the 100-half-year return level r
CI.low <- uniroot(obj, lower = I[1], upper = r,    k = k, x = M.hy, maxLogLik = mLL)
CI.up  <- uniroot(obj, lower = r,    upper = I[2], k = k, x = M.hy, maxLogLik = mLL)
(CI.r100 <- c(CI.low$root, CI.up$root)) # 95%-CI for 100-half-year return level r; [0.0504, 0.1912]
## Does it contain a drop as large as on Black Monday?
rBM <- as.numeric(X['1987-10-19']) # return level on Black Monday
CI.r100[1] <= rBM && rBM <= CI.r100[2] # => no!

## Note: Also here, note that this looks better than it actually is. Watch this:
I <- c(0.05, 0.25) # ... just a different initial interval
plot_obj(r = seq(I[1], I[2], length.out = 128), k = k, x = M.hy, maxLogLik = mLL)
abline(v = r)
## ... and problems of this sort appear frequently => always good to (at least)
## check with a plot whether the found roots make sense.


### 3.2.2 CI for the return period of a loss as on Black Monday ################

## Fix the return level r (as on Black Monday) and estimate the corresponding
## return period k
r <- rBM # fix return level
k <- 1/(1-pGEV(r, xi = xi.hy, mu = mu.hy, sigma = sig.hy)) # corresponding estimated return period k
## Find a suitable initial interval for computing the CI
I <- c(10, 1e5) # found by experimenting with the following plot
plot_obj(r = r, k = seq(I[1], I[2], length.out = 128), x = M.hy, maxLogLik = mLL)
abline(v = k) # estimated return level
## => As before
## Compute (at least the lower end of) the 95%-CI for the return period (in years)
## of an event of size as Black Monday
CI.low <- uniroot(obj, lower = I[1], upper = k, r = r, x = M.hy, maxLogLik = mLL)
CI.low$root
