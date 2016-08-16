## By Marius Hofert

## Compute the probability of next year's maximal risk-factor change
## exceeding all previous ones via the block maxima method based on
## S&P 500 data; see McNeil et al. (2015, Example 5.12). Compute a
## return level and return period; see McNeil et al. (2015, Example 5.15)

## Note: Different databases can (and indeed do) provide different stock values
##       for the S&P 500. The data below (from qrmdata) is from finance.yahoo.com
##       and thus our computed values differ from the values reported
##       in McNeil et al. (2015). Also, we (mostly) work with log-returns
##       instead of classical returns here.


### Setup ######################################################################

library(qrmdata) # for the S&P 500 data
library(qrmtools) # for log_returns
library(xts) # for functions around time series objects
library(QRM) # for fit.GEV()


### 1 Working with the data ####################################################

## Load the data and compute the negative log-returns (risk-factor changes X)
data(SP500)
S <- SP500 # 'xts'/'zoo' object
X <- -log_returns(S) # -log-returns X_t = -log(S_t/S_{t-1})
stopifnot(identical(X, -diff(log(S))[-1]))
date <- index(X)

## Let's briefly work out some numbers around next Monday (= Black Monday!)
X[date == "1987-10-16"] # ~= 5.3%; risk-factor change on the Friday before Black Monday
X[date == "1987-10-19"] # ~= 22.9%; risk-factor change on Black Monday!

## Let's briefly consider negative classical instead of -log-returns
## Note: A change of beta from yesterday's value to today's satisfies
##       (1+beta) * S_{t-1} = S_t => Y_t = -(S_t/S_{t-1}-1) = -beta
##       so Y_t gives exactly the drop beta (= -beta)
Y <- -diff(S)[-1]/as.numeric(S[-length(S)]) # classical negative returns
Y[date == "1987-10-16"] # ~= 5.16% (drop)
Y[date == "1987-10-19"] # ~= 20.47% (drop)

## To see the same change from -log-returns, note that X_t = -log(S_t/S_{t-1})
## = -log(1+beta) => -beta = -(exp(-X_t)-1) = -expm1(-X_t)
stopifnot(all.equal(-expm1(-X[date == "1987-10-16"]), Y[date == "1987-10-16"]))
stopifnot(all.equal(-expm1(-X[date == "1987-10-19"]), Y[date == "1987-10-19"]))

## Does working with either notion of (classical/log-)returns matter?
## The tangent to the curve log(x) in 1 is x - 1.
## 1) Near 1, log(x) ~ x - 1. So if S_t/S_{t-1} is roughly 1 (beta ~= 0),
##    then log(S_t/S_{t-1}) is roughly S_t/S_{t-1} - 1
##    => It does not matter much with which version one works
## 2) However, "Friday before Black Monday" S_t differs substantially from S_{t-1}
##    so it can/does matter:
x <- S["1987-10-19"]/as.numeric(S["1987-10-16"]) # S_t/S_{t-1} ~= 0.7953
stopifnot(all.equal(x - 1,  -Y[date == "1987-10-19"])) # classical return; ~= -0.2047
stopifnot(all.equal(log(x), -X[date == "1987-10-19"])) # log-return; ~= -0.2290
## => Since the tangent provides an upper bound for the log, classical returns
##    are larger than log-returns

## Drop from (end of) Mon 1987-10-12 to (end of) Fri 1987-10-16:
## Note: S_t/S_{t-4} = S_t/S_{t-1} * S_{t-1}/S_{t-2} ... * S_{t-3}/S_{t-4}
##       = exp(-X_t) * exp(-X_{t-1}) * ... * exp(-X_{t-3}) = exp(-sum(X_i, i=t-3,..,t))
##       => (Positive drop) -beta = -(S_t/S_{t-4}-1) = -(exp(-sum(X_i, i=t-3,..,t))-1)
-expm1(-sum(X["1987-10-12" <= date & date <= "1987-10-16"])) # ~= 9.12% (drop)

## From now on we only consider the risk-factor changes from 1960-01-01 until
## the evening of 1987-10-16
X. <- X["1960-01-01" <= date & date <= "1987-10-16"]

## Plot the S&P 500 log-returns
plot.zoo(X., main = "S&P 500 risk-factor changes (-log-returns)",
         xlab = "Time t", ylab = expression(X[t] == -log(S[t]/S[t-1])))


### 2 Block Maxima Method (BMM) and fitting the GEV ############################

## Extract (half-)yearly maxima method
M.year <- period.apply(X., INDEX = endpoints(X., "years"), FUN = max) # yearly maxima
endpts <- endpoints(X., "quarters") # end indices for quarters
endpts <- endpts[seq(1, length(endpts), by = 2)] # end indices for half-years
M.hyear <- period.apply(X., INDEX = endpts, FUN = max) # half-yearly maxima

## Fit the GEV distribution H_{xi,mu,sigma} to the (half-)yearly maxima
## Yearly maxima
fit.year <- fit.GEV(M.year)
(xi.year <- fit.year$par.ests[["xi"]]) # => ~= 0.2971 => Frechet domain with infinite 4th moment
(mu.year  <- fit.year$par.ests[["mu"]])
(sig.year <- fit.year$par.ests[["sigma"]])
fit.year$par.ses # standard errors
## Half-yearly maxima
fit.hyear <- fit.GEV(M.hyear)
(xi.hyear <- fit.hyear$par.ests[["xi"]]) # => ~= 0.3401 => Frechet domain with infinite 3rd moment
(mu.hyear  <- fit.hyear$par.ests[["mu"]])
(sig.hyear <- fit.hyear$par.ests[["sigma"]])
fit.hyear$par.ses # standard errors


### 3 Compute exceedance probabilities, return levels and return periods #######

## Q: What is the probability that next year's maximal risk-factor change
##    exceeds all previous ones?
prev.max <- max(M.year[-length(M.year)]) # maximum over all previous years
1-pGEV(prev.max, xi = xi.year,  mu = mu.year,  sigma = sig.year) # exceedance probability ~= 2.58%
1-pGEV(prev.max, xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear) # exceedance probability ~= 1.49%
## Note: mu and sig also differ for half-yearly vs yearly data; if it was only xi,
##       the exceedance probability based on half-yearly data would be larger
##       than the one based on yearly data (as xi is larger => heavier tailed GEV)

## Q: What is the 10-year and 50-year return level?
##    Recall: k n-block return level = r_{n,k} = H^-(1-1/k) = level which is
##            expected to be exceeded in one out of every k n-blocks.
qGEV(1-1/10, xi = xi.year, mu = mu.year, sigma = sig.year) # r_{n=260, k=10} ~= 4.42%; n ~ 1y
qGEV(1-1/50, xi = xi.year, mu = mu.year, sigma = sig.year) # r_{n=260, k=50} ~= 7.49%
## 20-half-year and 100-half-year return levels
qGEV(1-1/20,  xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear) # r_{n=130, k=20}  ~= 4.56%; n ~ 1/2y
qGEV(1-1/100, xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear) # r_{n=130, k=100} ~= 7.90%
## => Close to r_{n = 260, k = 10} and r_{n = 260, k = 50}, respectively (reassuring).

## Q: What is the return period of a risk-factor change at least as large as
##    on Black Monday?
##    Recall: k_{n,u} = 1/\bar{H}(u) = period (= number of n-blocks) in which we
##            expect to see a single n-block exceeding u (= risk-factor change
##            as on Black Monday)
1/(1-pGEV(as.numeric(X[date == "1987-10-19"]),
          xi = xi.year, mu = mu.year, sigma = sig.year)) # ~= 1877 years
1/(1-pGEV(as.numeric(X[date == "1987-10-19"]),
          xi = xi.hyear, mu = mu.hyear, sigma = sig.hyear)) # ~= 2300 half-years = 1150 years

