## By Alexander McNeil and Marius Hofert

## This script calculates confidence intervals for 05_GEV_BMM_SP500.R
## by the profile likelihood method.


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
##'         In other words, pLL()$value gives the profile log-likelihood
##'         in r or k.
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
##' @note 1) Background information (see Davison (2003, p. 126)):
##'          Under regularity conditions, the likelihood ratio statistic
##'          W(th) = 2(LL(MLE) - LL(th)) evaluated at the true underlying, unknown
##'          th_0 converges in distribution to a chi_p^2 where p = length(th_0).
##'          An asymptotic (1-alpha)-CI for th_0 is thus {th : W(th) <= (chi_p^2)^-(1-alpha)}
##'          = {th : LL(th) >= LL(MLE) - (chi_p^2)^-(1-alpha) / 2}. Replacing ">="
##'          by "=", we have to find the two roots
##'          LL(th) - (LL(MLE) - (chi_p^2)^-(1-alpha) / 2) = 0
##'          (corresponding to the lower and upper CI endpoints).
##'       2) Idea profile log-likelihoods:
##'          If only (CIs for) some parameters th.i are of interest (the
##'          'parameters of interest'; the remaining ones are the 'nuisance
##'          parameters' th.n), one can define the 'generalized likelihood ratio
##'          statistic' W(th.i) = 2(LL(MLE) - LL(th.i, hat(th.n))) where
##'          hat(th.n) is the maximizer of the 'profile log-likelihood'
##'          pLL(th.i) = LL(th.i, th.n) in th.n for given th.i. Under regularity
##'          conditions, W(th.i) converges in distribution to a chi_{p.i}^2 where
##'          p.i = length(th_i). An asymptotic (1-alpha)-CI for the true underlying,
##'          unknown th_{i,0} is thus {th_i : W(th_i) <= (chi_{p.i}^2)^-(1-alpha)}
##'          = {th_i : LL(th.i, hat(th.n)) >= LL(MLE) - (chi_{p.i}^2)^-(1-alpha) / 2}.
##'          Replacing ">=" by "=", we have to find the two roots
##'          LL(th.i, hat(th.n)) - (LL(MLE) - (chi_{p.i}^2)^-(1-alpha) / 2) = 0
##'          (corresponding to the lower and upper CI endpoints).
##'       3) Our parameter of interest is either r or k (by fixing one, the
##'          other is fixed, too). The nuisance parameters are xi and sigma,
##'          and mu is implied when fixing r, k, xi, sigma; see minLL_r_k()
root <- function(r, k, x, mLL, alpha = 0.05, ...)
    (-npLL(r, k = k, x = x, ...)$value) - # profile log-likelihood
        (mLL - qchisq(1-alpha, df = 1) / 2)
## Note: One could also write a function (say, 'find_CI()') which calls uniroot()
##       based on root() two times to find the left and right CI endpoints.
##       find_CI() could also start by computing the overall MLE and then
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

## Fit GEV to yearly maxima and compute MLE
fit.year <- fit.GEV(M.year) # likelihood-based estimation of the GEV; see ?fit.GEV
(xi.year <- fit.year$par.ests[["xi"]]) # => ~= 0.2971 => Frechet domain with infinite 4th moment
(mu.year  <- fit.year$par.ests[["mu"]])
(sig.year <- fit.year$par.ests[["sigma"]])
fit.year$par.ses # standard errors
mLL <- fit.year$llmax # maximum log-likelihood (at MLE)


### 3.1.1 CI for 10-year and 50-year return level ##############################

## Fix return period k = 10 (in years) and estimate the corresponding return
## level r
k <- 10 # return period
r <- qGEV(1-1/k, xi = xi.year, mu = mu.year, sigma = sig.year) # corresponding estimated return level r
## Compute 95%-CI for 10-year return level r
low <- uniroot(root, lower = 1e-6, upper = r,   k = k, x = M.year, mLL = mLL)
up  <- uniroot(root, lower = r,    upper = 1e2, k = k, x = M.year, mLL = mLL)
c(low$root, up$root)

## Fix return period k = 50 (in years) and estimate the corresponding return
## level r
k <- 50 # return period
r <- qGEV(1-1/k, xi = xi.year, mu = mu.year, sigma = sig.year) # corresponding estimated return level r
## Compute 95%-CI for 50-year return level r
low <- uniroot(root, lower = 1e-6, upper = r,   k = k, x = M.year, mLL = mLL)
up  <- uniroot(root, lower = r,    upper = 1e2, k = k, x = M.year, mLL = mLL)
c(low$root, up$root)


### 3.1.2 CI for the return period of a loss as on Black Monday

## Fix return level r (as on Black Monday) and estimate the corresponding return period k
r <- as.numeric(X['1987-10-19']) # return level
k <- 1/(1-pGEV(r, xi = xi.year, mu = mu.year, sigma = sig.year))
## Compute 95%-CI for return period (in years) of an event of size as Black Monday
low <- uniroot(root, lower = 1+1e-2, # has to be > 1 (for log1p(-1/k) to not be NaN)
               upper = k, r = r, x = M.year, mLL = mLL)
## Upper bound of the 95%-CI not easy here because of numerical problems:
pwr <- seq(1, 17, by = 0.2) # powers of 10
profile <- sapply(pwr, function(p) root(r = r, k = 10^p, x = M.year, mLL = mLL,
                                        control = list(reltol = 0)))
plot(10^pwr, profile, type = "l", log = "x", xlab = "return period")
abline(h = 0)
## Note: Useless without reltol = 0. Even with, the objective function
##       is quite erratic as it depends on computed optima. Problem probably
##       requires a rescaling of the likelihood etc. to be more tractable.


### 3.2 Based on biannual maxima ###############################################

## Fit GEV to half-yearly maxima
fit.hyear <- fit.GEV(M.hyear)
(xi.hyear <- fit.hyear$par.ests[["xi"]]) # => ~= 0.3401 => Frechet domain with infinite 3rd moment
(mu.hyear  <- fit.hyear$par.ests[["mu"]])
(sig.hyear <- fit.hyear$par.ests[["sigma"]])
fit.hyear$par.ses # standard errors
mLL <- fit.hyear$llmax # maximum log-likelihood (at MLE)


### 3.2.1 CI for 20-hyear return level #########################################

## Fix return period k = 20 (in half-years) and estimate the corresponding
## return level r
k <- 20 # return period
r <- qGEV(1-1/k, xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear) # corresponding estimated return level r
## Compute 95%-CI for 20-year return level r
low <- uniroot(root, lower = 1e-6, upper = r,   k = k, x = M.hyear, mLL = mLL)
up  <- uniroot(root, lower = r,    upper = 1e2, k = k, x = M.hyear, mLL = mLL)
c(low$root, up$root)


### 3.2.2 CI for the return period of a loss as on Black Monday ################

## Fix return level r (as on Black Monday) and estimate the corresponding return period k
r <- as.numeric(X['1987-10-19']) # return level
k <- 1/(1-pGEV(r, xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear))
## Compute 95%-CI for return period (in half-years) of an event of size as Black Monday
low <- uniroot(root, lower = 1+1e-2, # has to be > 1 (for log1p(-1/k) to not be NaN)
               upper = k, r = r, x = M.hyear, mLL = mLL)
up  <- uniroot(root, lower = k, # has to be > 1 (for log1p(-1/k) to not be NaN)
               upper = 10^8, r = r, x = M.hyear, mLL = mLL)
c(low$root, up$root)
