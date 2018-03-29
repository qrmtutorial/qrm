## By Marius Hofert

## 1) Compute the probability of next year's maximal risk-factor change
##    exceeding all previous ones via the BBM based on S&P 500 data;
##    see McNeil et al. (2015, Example 5.12).
## 2) Compute a return level;  see McNeil et al. (2015, Example 5.15)
## 3) Compute a return period; see McNeil et al. (2015, Example 5.15)

## Note: Different databases can (and indeed do) provide different stock values
##       for the S&P 500. The data below (from qrmdata) is from finance.yahoo.com
##       and thus our computed values differ from the values reported
##       in McNeil et al. (2015). Also, we (mostly) work with log-returns
##       instead of classical returns here.


### Setup ######################################################################

library(xts) # for functions around time series objects
library(qrmdata) # for the S&P 500 data
library(qrmtools) # for returns()


### 1 Working with the data ####################################################

## Load the data and compute the negative log-returns (risk-factor changes X)
data(SP500) # load the S&P 500 data
S <- SP500 # 'xts'/'zoo' object
X <- -returns(S) # -log-returns X_t = -log(S_t/S_{t-1})
stopifnot(all.equal(X, -diff(log(S))[-1], check.attributes = FALSE))

## Let's briefly work out some numbers around next Monday (= Black Monday!)
X['1987-10-16'] # ~=  5.3%; risk-factor change on the Friday before Black Monday
X['1987-10-19'] # ~= 22.9%; risk-factor change on Black Monday!

## Let's briefly consider negative classical instead of -log-returns
## Note: A change of beta from yesterday's value to today's satisfies
##       S_t = (1+beta) * S_{t-1} => Y_t = -(S_t/S_{t-1}-1) = -beta
##       => The negative classical returns Y_t give exactly the drop beta (= -beta)
Y <- -returns(S, method = "simple") # classical negative returns
stopifnot(all.equal(Y, -diff(S)[-1]/as.numeric(S[-length(S)]),
                    check.attributes = FALSE))
Y['1987-10-16'] # ~=  5.16% (drop)
Y['1987-10-19'] # ~= 20.47% (drop)

## To see the same change from -log-returns, note that X_t = -log(S_t/S_{t-1})
## = -log(1+beta) => -beta = -(exp(-X_t)-1) = -expm1(-X_t),
## so negative classical returns can be obtained from -log-returns via -expm1(-.)
stopifnot(all.equal(-expm1(-X['1987-10-16']), Y['1987-10-16'], check.attributes = FALSE))
stopifnot(all.equal(-expm1(-X['1987-10-19']), Y['1987-10-19'], check.attributes = FALSE))
## ... and over a time period: The drop (= loss) from (end of) Mon 1987-10-12
## to (end of) Fri 1987-10-16 can be obtained via -expm1(-sum(.)):
## Note: S_t/S_{t-4} = S_t/S_{t-1} * S_{t-1}/S_{t-2} ... * S_{t-3}/S_{t-4}
##       = exp(-X_t) * exp(-X_{t-1}) * ... * exp(-X_{t-3}) = exp(-sum(X_i, i=t-3,..,t))
##       => (Positive drop) -beta = -(S_t/S_{t-4}-1) = -(exp(-sum(X_i, i=t-3,..,t))-1)
-expm1(-sum(X['1987-10-12/1987-10-16'])) # ~= 9.12% (drop)

## Does working with either notion of (classical/log-)returns matter?
## The tangent to the curve log(x) in 1 is x - 1.
## 1) Near 1, log(x) ~ x - 1. So if S_t/S_{t-1} ~= 1 (so beta ~= 0),
##    then log(S_t/S_{t-1}) is roughly S_t/S_{t-1} - 1
##    => It does not matter much with which version one works
## 2) However, "Friday before Black Monday" S_t differs substantially from S_{t-1}
##    so it can/does matter:
(x <- S["1987-10-19"]/as.numeric(S["1987-10-16"])) # S_t/S_{t-1} ~= 0.7953
stopifnot(all.equal(-(x - 1), Y['1987-10-19'], # classical -return; ~= 0.2047
                    check.attributes = FALSE))
stopifnot(all.equal(-log(x),  X['1987-10-19'], # -log-return;       ~= 0.2290
                    check.attributes = FALSE))
## => Difference of 2.43%.
##    Since log(x) <= x - 1, classical returns are larger than log-returns
##    and so -log-returns are larger than negative classical returns.

## From now on we only consider the -log-returns from 1960-01-01 until
## the evening of 1987-10-16
X. <- X['1960-01-01/1987-10-16']

## Plot the S&P 500 -log-returns
plot.zoo(X., main = "S&P 500 risk-factor changes (-log-returns)",
         xlab = "Time t", ylab = expression(X[t] == -log(S[t]/S[t-1])))
## One would need to fit a time series model to X. (e.g., GARCH process), but
## we omit that here.


### 2 Block Maxima Method (BMM) and fitting the GEV ############################

## Extract (half-)yearly maxima method
M.y <- period.apply(X., INDEX = endpoints(X., "years"), FUN = max) # yearly maxima
endpts <- endpoints(X., "quarters") # end indices for quarters
endpts <- endpts[seq(1, length(endpts), by = 2)] # end indices for half-years
M.hy <- period.apply(X., INDEX = endpts, FUN = max) # half-yearly maxima

## Fit the GEV distribution H_{xi,mu,sigma} to the (half-)yearly maxima

## Yearly maxima
fit.y <- fit_GEV(M.y) # likelihood-based estimation of the GEV; see ?fit.GEV
stopifnot(fit.y$convergence == 0) # => converged
(xi.y <- fit.y$par[1]) # ~= 0.2972 => Frechet domain with infinite ceiling(1/xi.y) = 4th moment
(mu.y  <- fit.y$par[2])
(sig.y <- fit.y$par[3])
sqrt(diag(fit.y$Cov)) # standard errors

## Half-yearly maxima
fit.hy <- fit_GEV(M.hy)
stopifnot(fit.hy$convergence == 0) # => converged
(xi.hy <- fit.hy$par[1]) # ~= 0.3401 => Frechet domain with infinite ceiling(1/xi.hy) = 3rd moment
(mu.hy  <- fit.hy$par[2])
(sig.hy <- fit.hy$par[3])
sqrt(diag(fit.hy$Cov)) # standard errors


### 3 Compute exceedance probabilities, return levels and return periods #######

## Q: What is the probability that next year's maximal risk-factor change
##    exceeds all previous ones?
1-pGEV(max(head(M.y, n = -1)),  xi = xi.y,  mu = mu.y,  sigma = sig.y) # exceedance prob. ~= 2.58%
1-pGEV(max(head(M.hy, n = -1)), xi = xi.hy, mu = mu.hy, sigma = sig.hy) # exceedance prob. ~= 1.48%
## Note: mu and sig also differ for half-yearly vs yearly data; if it was only xi,
##       the exceedance probability based on half-yearly data would be estimated
##       larger than the one based on yearly data (as xi is larger => heavier tailed GEV)

## Q: What is the 10-year and 50-year return level?
##    Recall: k n-block return level = r_{n,k} = H^-(1-1/k) = level which is
##            expected to be exceeded in one out of every k n-blocks.
qGEV(1-1/10, xi = xi.y, mu = mu.y, sigma = sig.y) # r_{n = 260, k = 10} ~= 4.42%; n ~ 1y
qGEV(1-1/50, xi = xi.y, mu = mu.y, sigma = sig.y) # r_{n = 260, k = 50} ~= 7.49%
## 20-half-year and 100-half-year return levels
qGEV(1-1/20,  xi = xi.hy, mu = mu.hy, sigma = sig.hy) # r_{n = 130, k = 20}  ~= 4.56%; n ~ 1/2y
qGEV(1-1/100, xi = xi.hy, mu = mu.hy, sigma = sig.hy) # r_{n = 130, k = 100} ~= 7.90%
## => Close to r_{n = 260, k = 10} and r_{n = 260, k = 50}, respectively (reassuring).

## Q: What is the return period of a risk-factor change at least as large as
##    on Black Monday?
##    Recall: k_{n,u} = 1/\bar{H}(u) = period (= number of n-blocks) in which we
##            expect to see a single n-block exceeding u (= risk-factor change
##            as on Black Monday)
1/(1-pGEV(as.numeric(X['1987-10-19']),
          xi = xi.y, mu = mu.y, sigma = sig.y)) # ~= 1873 years
1/(1-pGEV(as.numeric(X['1987-10-19']),
          xi = xi.hy, mu = mu.hy, sigma = sig.hy)) # ~= 2302 half-years = 1151 years


### 4 Behind the scenes ########################################################

## Note:
## 1) Fitting a GEV via MLE needs initial values. One (of multiple) problem(s) is
##    that for xi != 0, we need 1 + xi * (x - mu) / sigma > 0 for *all* data points
##    x in order for the density to be positive and thus the log-likelihood to be
##    finite.
## 2) A cheap solution is to use the case xi = 0 as initial values. So xi = 0,
##    sigma = sqrt(6 * var(x)) / pi and mu = mean(x) - sigma * gam in this case,
##    where gam = 0.5772157 (= -digamma(1)) denotes the Euler--Mascheroni constant.
## 3) This is no problem for x = M.y, however, it is a problem for x = M.hy as we
##    now explore a bit.
## 4) These types of numerical explorations are important for many practical
##    problems. It is always good to keep an eye on numerics...

## Fit based on half-yearly maxima with the case 'xi = 0' as initial values
x <- M.hy # for simplicity
sig.init <- sqrt(6 * var(x)) / pi # initial sigma for the 'xi = 0' case
mu.init <- mean(x) - sig.init * 0.5772157
init <- c(0, mu.init, sig.init) # initial xi, mu, sigma
fit.hy. <- fit_GEV(M.hy, init = init) # likelihood-based estimation of the GEV; see ?fit.GEV
stopifnot(fit.hy.$convergence == 0) # => converged, but...
(xi.hy. <- fit.hy.$par[1]) # => ~= -7e-17, suspicious
if(FALSE) {
    ## Let's look at a trace
    fit_GEV(M.hy, init = init, control = list(trace = TRUE)) # => nothing obvious
    ## On a second thought (comparing the log-likelihoods below (~ 191) with the
    ## ones of this trace), Nelder--Mead seems to 'jump' on the second reflection.
}

## Check a different optimizer
fit_GEV(M.hy, init = init, method = "BFGS") # ~= 0.3265 => ok

## Check a whole sequence of xi's (for the same mu, sigma)
s <- 2^(-seq(3, 20, by = 0.5)) # sequence concentrating near 0
xi.init <- c(-s, 0, rev(s)) # sequence concentrating around 0
fit.hy.xis <- lapply(xi.init, function(xi) fit_GEV(M.hy, init = c(xi, init[2:3])))
stopifnot(sapply(fit.hy.xis, function(f) f$convergence) == 0) # => converged
fitted.xi <- lapply(fit.hy.xis, function(f) f$par[1]) # fitted xi's
plot(xi.init, fitted.xi, type = "b",
     xlab = expression("initial"~xi), ylab = expression("fitted"~xi))
## => The optimization method Nelder--Mead fails exactly for xi = 0...

## ... with no obvious reason (at least not when looking at the log-likelihood
## as a function of xi
logL <- sapply(xi.init, function(xi) logLik_GEV(c(xi, init[2:3]), M.hy))
plot(xi.init, logL, type = "b") # => ok; see also 'logL'

## For initial |xi| small but nonzero, fitting works.
fit_GEV(M.hy, init = c( .Machine$double.eps, init[2:3])) # ~= 0.3401 => ok
fit_GEV(M.hy, init = c(-.Machine$double.eps, init[2:3])) # ~= 0.3401 => ok
## ... and this and another trick are in fact used internally.
## In some other R packages, initial xi are more brutally set as 0.01.

## Note:
## 1) fit_GEV() internally also provides more algorithms for finding
##    initial values, and they also (roughly) lead to xi ~= 0.3402.
## 2) The problem overall is significantly more severe and challenging
##    when computing profile log-likelihood based confidence intervals
##    for return levels and return periods; see the script
##    05_GEV_BMM_SP500_CIs_rlevels_rperiods.R
