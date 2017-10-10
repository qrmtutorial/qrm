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
library(QRM) # for fit.GEV()
library(qrmdata) # for the S&P 500 data
library(qrmtools) # for returns()


##' @title Negative Profile Log-Likelihood
##' @param r return level (parameter of interest)
##' @param k return period (parameter of interest)
##' @param theta initial values for the nuisance parameters (xi, sigma)
##' @param x maxima
##' @param ... additional arguments passed to optim()
##' @return return value of optim() with the minimizer of the -log-likelihood
##'         in the nuisance parameters theta = (xi, sigma) for fixed r or k.
##'         In other words, npLL()$value gives the negative profile
##'         log-likelihood in r or k.
##' @author Alexander McNeil and Marius Hofert
##' @note r = H^-(1-1/k) = mu + (sig/xi)((-log(1-1/k))^xi - 1) so the mu
##'       implied by a fixed xi, sigma, r, k is
##'       mu = r - (sig/xi)((-log(1-1/k))^xi - 1)
npLL <- function(r, k, theta = c(xi = 0.01, sigma = sqrt(6*var(x)/pi)), x, ...)
{
    nLL <- function(theta) { # note: r, k and x are locally 'seen' here
        xi <- theta[1]
        sig <- theta[2]
        implied.mu <- r - (sig/xi) * ((-log1p(-1/k))^(-xi) - 1) # mu implied by xi, sigma, r, k
        -sum(dGEV(x, xi = xi, mu = implied.mu, sigma = abs(sig), log = TRUE)) # -log-likelihood
    }
    ## Minimize the -log-likelihood in the nuisance parameters theta = (xi, sigma)
    ## for our given parameters of interest r and k (and the data x)
    optim(theta, fn = nLL, ...) # minimize nLL
}

##' @title Computing a Side of the (1-alpha)-Confidence Interval for the
##'        Return Level r or the Return Period k
##' @param r return level
##' @param k return period
##' @param x maxima
##' @param mLL maximal (overall) log-likelihood (when maximizing GEV in all parameters)
##' @param alpha significance level
##' @param ... additional arguments passed to npLL()
##' @return One side of the confidence level (which one is determined by
##'         the initial values)
##' @author Alexander McNeil and Marius Hofert
##' @note - Our parameter of interest is either r or k (fix the other).
##'       - The nuisance parameters are xi and sigma; mu is implied when
##'         fixing r, k, xi, sigma
root <- function(r, k, x, mLL, alpha = 0.05, ...)
    (-npLL(r, k = k, x = x, ...)$value) - # profile +log-likelihood
        (mLL - qchisq(1-alpha, df = 1) / 2)
## Note: One could also write a function (say, 'find_CI()') which calls uniroot()
##       based on root() two times to find the left and right CI endpoints.
##       find_CI() could also start from computing the overall MLE and then
##       extend the likelihood to the left and right (by the usual 'doubling')
##       to find the roots.


### 1 Working with the data ####################################################

## Load the data and compute the negative log-returns (risk-factor changes X)
data(SP500)
S <- SP500 # 'xts'/'zoo' object
X <- -returns(S) # -log-returns X_t = -log(S_t/S_{t-1})
X. <- X['1960-01-01/1987-10-16'] # grab out data we work with

## Extract (half-)yearly maxima method
M.year <- period.apply(X., INDEX = endpoints(X., "years"), FUN = max) # yearly maxima
endpts <- endpoints(X., "quarters") # end indices for quarters
endpts <- endpts[seq(1, length(endpts), by = 2)] # end indices for half-years
M.hyear <- period.apply(X., INDEX = endpts, FUN = max) # half-yearly maxima


### 3 Profile-likelihood-based CIs for return levels and return periods ########

### 3.1 Based on annual maxima #################################################

## Fit GEV to yearly maxima
fit.year <- fit.GEV(M.year) # likelihood-based estimation of the GEV; see ?fit.GEV
(xi.year <- fit.year$par.ests[["xi"]]) # => ~= 0.2971 => Frechet domain with infinite 4th moment
(mu.year  <- fit.year$par.ests[["mu"]])
(sig.year <- fit.year$par.ests[["sigma"]])

fit.year$par.ses # standard errors
mLL <- fit.year$llmax # maximum log-likelihood (at MLE)


### 3.1.1 CI for 10-year and 50-year return level ##############################

## Fix the return period k = 10 (in years) and estimate the corresponding
## return level r
k <- 10 # fix return period
r <- qGEV(1-1/k, xi = xi.year, mu = mu.year, sigma = sig.year) # corresponding estimated return level r
## Compute 95%-CI for the 10-year return level r
low <- uniroot(root, lower = 1e-6, upper = r,   k = k, x = M.year, mLL = mLL)
up  <- uniroot(root, lower = r,    upper = 1e2, k = k, x = M.year, mLL = mLL)
(CI.r10 <- c(low$root, up$root)) # 95%-CI for 10-year return level r
## Does it contain a drop as large as on Black Monday?
rBM <- as.numeric(X['1987-10-19']) # return level on Black Monday
CI.r10[1] <= rBM && rBM <= CI.r10[2] # => no!

## The same for k = 50
k <- 50 # fix return period
r <- qGEV(1-1/k, xi = xi.year, mu = mu.year, sigma = sig.year) # corresponding estimated return level r
## Compute 95%-CI for the 50-year return level r
low <- uniroot(root, lower = 1e-6, upper = r,   k = k, x = M.year, mLL = mLL)
up  <- uniroot(root, lower = r,    upper = 1e2, k = k, x = M.year, mLL = mLL)
(CI.r50 <- c(low$root, up$root)) # 95%-CI for 50-year return level r
## Does it contain a drop as large as on Black Monday?
CI.r50[1] <= rBM && rBM <= CI.r50[2] # => yes!


### 3.1.2 CI for the return period of a loss as on Black Monday ################

## Fix the return level r (as on Black Monday) and estimate the corresponding
## return period k
r <- rBM # fix return level
k <- 1/(1-pGEV(r, xi = xi.year, mu = mu.year, sigma = sig.year)) # corresponding estimated return period k
## Compute 95%-CI for the return period (in years) of an event of size as Black Monday
low <- uniroot(root, lower = 1+1e-2, # has to be > 1 (for log1p(-1/k) to not be NaN)
               upper = k, r = r, x = M.year, mLL = mLL)
## Upper bound of the 95%-CI difficult here because of numerical problems:
pwr <- seq(1, 17, by = 0.2) # powers of 10
profile <- sapply(pwr, function(p) root(r = r, k = 10^p, x = M.year, mLL = mLL,
                                        control = list(reltol = 0)))
## Plot the objective function for finding the root of for computing the upper
## 95%-CI for the return period of a loss with return level r as on Black Monday
plot(10^pwr, profile, # each function evaluation is a numerically found root (not stable)
     type = "l", log = "x", xlab = "Return period k",
     ylab = "Function for finding the upper 95%-CI bound for k (of event as on BM)")
abline(h = 0)
## Note: - Useless without reltol = 0
##       - Even with reltol = 0, the objective function is quite erratic
##         as it depends on computed optima.
##       - Problem probably requires a rescaling of the likelihood.


### 3.2 Based on biannual maxima ###############################################

## Fit GEV to half-yearly maxima
fit.hyear <- fit.GEV(M.hyear)
(xi.hyear <- fit.hyear$par.ests[["xi"]]) # => ~= 0.3401 => Frechet domain with infinite 3rd moment
(mu.hyear  <- fit.hyear$par.ests[["mu"]])
(sig.hyear <- fit.hyear$par.ests[["sigma"]])

fit.hyear$par.ses # standard errors
mLL <- fit.hyear$llmax # maximum log-likelihood (at MLE)


### 3.2.1 CI for 20-half-year and 100-half-year return level ###################

## Fix the return period k = 20 (in half-years) and estimate the corresponding
## return level r
k <- 20 # fix return period
r <- qGEV(1-1/k, xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear) # corresponding estimated return level r
## Compute 95%-CI for the 20-half-year return level r
low <- uniroot(root, lower = 1e-6, upper = r,   k = k, x = M.hyear, mLL = mLL)
up  <- uniroot(root, lower = r,    upper = 1e2, k = k, x = M.hyear, mLL = mLL)
(CI.r20 <- c(low$root, up$root)) # 95%-CI for 20-half-year return level r
## Does it contain a drop as large as on Black Monday?
CI.r20[1] <= rBM && rBM <= CI.r20[2] # => no!

## The same for k = 100
k <- 100
r <- qGEV(1-1/k, xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear) # corresponding estimated return level r
## Compute 95%-CI for the 100-half-year return level r
low <- uniroot(root, lower = 1e-6, upper = r,   k = k, x = M.hyear, mLL = mLL)
up  <- uniroot(root, lower = r,    upper = 1e2, k = k, x = M.hyear, mLL = mLL)
(CI.r100 <- c(low$root, up$root)) # 95%-CI for 100-half-year return level r
## Does it contain a drop as large as on Black Monday?
CI.r100[1] <= rBM && rBM <= CI.r100[2] # => no (not quite!)


### 3.2.2 CI for the return period of a loss as on Black Monday ################

## Fix the return level r (as on Black Monday) and estimate the corresponding
## return period k
r <- rBM # fix return level
k <- 1/(1-pGEV(r, xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear)) # corresponding estimated return period k
## Compute 95%-CI for the return period (in half-years) of an event of size as Black Monday
low <- uniroot(root, lower = 1+1e-2, # has to be > 1 (for log1p(-1/k) to not be NaN)
               upper = k, r = r, x = M.hyear, mLL = mLL)
up  <- uniroot(root, lower = k, # has to be > 1 (for log1p(-1/k) to not be NaN)
               upper = 10^8, r = r, x = M.hyear, mLL = mLL)
(CI.kBM <- c(low$root, up$root)) # 95%-CI (in half-years) for return period k of event as on Black Monday
