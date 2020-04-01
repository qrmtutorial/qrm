## By Marius Hofert

## R script for Chapter 3 of The QRM Exercise Book


### Setup ######################################################################

library(zoo) # for na.fill()
library(xts)
library(tseries)
library(qrmtools)
library(qrmdata)
library(moments)
doPDF <- require(crop)
options(digits = 10)


### Exercise 3.2 (Negative log-returns and the correlogram) ####################

## Reproducing code for exercise

## Data preparation
data(SP500)
S <- SP500['2011-01-03/2015-12-31']
X <- -returns(S)
n <- nrow(X)

## Plot of the negative log-returns
if(doPDF) pdf(file = (file <- "fig_03_returns_correlograms_-log_returns.pdf"))
opar <- par(pty = "s")
plot.zoo(X, xlab = "Time", ylab = "-Log-returns of the S&P 500 index")
mtext(substitute("Trading days:"~n., list(n. = n)),
      side = 4, adj = 0, line = 0.5)
par(opar)
if(doPDF) dev.off.crop(file)

## Corresponding correlogram
if(doPDF) pdf(file = (file <- "fig_03_returns_correlograms_ACF.pdf"))
opar <- par(pty = "s")
acf(X, ci.col = "black", main = "", ylab = "ACF of -log-returns")
par(opar)
if(doPDF) dev.off.crop(file)


## c) Reproducing code

## Plot of the squared negative log-returns
if(doPDF) pdf(file = (file <- "fig_03_returns_correlograms_-log_returns_squared.pdf"))
opar <- par(pty = "s")
plot.zoo(X^2, xlab = "Time", ylab = "Squared -log-returns of the S&P 500 index")
mtext(substitute("Trading days:"~n., list(n. = n)),
      side = 4, adj = 0, line = 0.5)
par(opar)
if(doPDF) dev.off.crop(file)

## Corresponding correlogram
if(doPDF) pdf(file = (file <- "fig_03_returns_correlograms_ACF_squared.pdf"))
opar <- par(pty = "s")
acf(X^2, ci.col = "black", main = "", ylab = "ACF of squared -log-returns")
par(opar)
if(doPDF) dev.off.crop(file)


### Exercise 3.3 (Interpreting Q-Q plots) ######################################

## Reproducing code for exercise

## Data preparation
data(crypto)
S <- crypto['2014-05-29/2018-05-29',"BTC"]
X <- -returns(S)
plot.zoo(X)

## Fit a normal distribution
mu <- mean(X) # fitted mean
sig <- sd(X) # fitted standard deviation
n <- nrow(X) # sample size

## Q-Q plot
if(doPDF) pdf(file = (file <- "fig_03_S-shape_in_Q-Q_plots_BTC_vs_fitted_normal_Q-Q.pdf"))
opar <- par(pty = "s")
qq_plot(X, FUN = function(p) qnorm(p, mean = mu, sd = sig),
        xlab = "Quantiles of the fitted normal distribution",
        ylab = "Sample quantiles of -log-returns of Bitcoin")
mtext(substitute("Location:"~mu.*", scale:"~sig.*", sample size:"~n.,
                 list(mu. = round(mu, 4), sig. = round(sig, 4), n. = n)),
      side = 4, adj = 0, line = 0.5)
par(opar)
if(doPDF) dev.off.crop(file)


### Exercise 3.5 (Joint negative log-returns and cross-correlogram) ############

## Reproducing code for exercise

## Data preparation
data(SP500)
data(DAX)
S <- merge(SP500, DAX)['2011-01-03/2015-12-31']
numNA <- rowSums(is.na(S))
any(numNA > 1) # => at most one data point missing
sum(numNA == 1) # => 46 dates with missing data points
S <- na.fill(S, fill = "extend") # fill NAs (by extension)
colnames(S) <- c("S&P 500", "DAX")
X <- -returns(S)
n <- nrow(X)
X. <- as.matrix(X)

## Plot of the joint negative log-returns
if(doPDF) pdf(file = (file <- "fig_03_returns_cross-correlograms_-log_returns.pdf"))
opar <- par(pty = "s")
plot(X., xlab = "-Log-returns of the S&P 500 index", ylab = "-Log-returns of the DAX index")
mtext(substitute("Trading days:"~n., list(n. = n)),
      side = 4, adj = 0, line = 0.5)
par(opar)
if(doPDF) dev.off.crop(file)

## Plot of the corresponding scaled ranks
if(doPDF) pdf(file = (file <- "fig_03_returns_cross-correlograms_-log_returns_pobs.pdf"))
opar <- par(pty = "s")
plot(apply(X., 2, rank)/(n+1), xlab = "Scaled ranks of -log-returns of the S&P 500 index",
     ylab = "Scaled ranks of -log-returns of the DAX index")
mtext(substitute("Trading days:"~n., list(n. = n)),
      side = 4, adj = 0, line = 0.5)
par(opar)
if(doPDF) dev.off.crop(file)


## c) Reproducing code

## Corresponding cross-correlogram (showing Cor(X_{t+h}, Y_t); see ?acf -> 'Value')
if(doPDF) pdf(file = (file <- "fig_03_returns_cross-correlograms_CCF.pdf"))
opar <- par(pty = "s")
ccf(as.ts(X[,1]), as.ts(X[,2]), ci.col = "black", main = "", ylab = "CCF of -log-returns") # note: ccf(X[,1], X[,2]) fails
if(FALSE)
    acf(X) # => top-left/bottom-right overlaid at h = 0 are the CCF; this is what ccf() does after calling acf()
par(opar)
if(doPDF) dev.off.crop(file)

## Cross-correlogram for squared negative log-returns
if(doPDF) pdf(file = (file <- "fig_03_returns_cross-correlograms_CCF_squared.pdf"))
opar <- par(pty = "s")
ccf(as.ts(X[,1]^2), as.ts(X[,2]^2), ci.col = "black", main = "", ylab = "CCF of squared -log-returns")
par(opar)
if(doPDF) dev.off.crop(file)


### Exercise 3.8 (Experimenting with Q-Q plots) ################################

n <- 1000
set.seed(271)

## a)
qq_plot(rt(n, df = 6), method = "empirical")

## b)
qq_plot(runif(n), method = "empirical")

## c)
qq_plot(rlogis(n), method = "empirical")

## d)
qq_plot(exp(rnorm(n)), method = "empirical")

## Note:
## 1) Obviously, all Q-Q plots indicate a departure from normality.
##    Note that each Q-Q plot actually checks whether the provided sample comes
##    from a type of standard normal distribution, so a location-scale transformed
##    standard normal (and thus an ordinary normal) distribution. To check for an
##    exact standard normal, omit 'method = "empirical"'; this would only indicate
##    even bigger departures from the hypothesized (the standard normal)
##    distribution.
## 2) The 'inverted-S' shape in a) shows that the empirical distribution (of
##    the sample) has heavier tails than the theoretical distribution compared
##    against (the normal). The 'S' shape in b) shows the opposite. Both
##    observations are in line with what we expect given the distributions from
##    which we sampled in these cases.


### Exercise 3.9 (Ljung--Box tests) ############################################

## a)

## Data preparation
data(EURSTX_const)
S.raw <- EURSTX_const['2000-01-01/2015-12-31']

## Dealing with missing data (keep those not having missing data on first day
## and fill the remaining NAs by interpolation)
all(apply(is.na(S.raw), 2, any)) # => all stocks have some missing data
NA_plot(S.raw) # see missing data
keep <- !is.na(S.raw[1,]) # indicator which stocks to keep
S <- S.raw[,keep] # keep those stocks
S <- na.fill(S, fill = "extend") # 'fill' the missing data
NA_plot(S) # visual check
stopifnot(!is.na(S)) # no more missing data

## Log-returns
X.d <- returns(S) # daily
X.m <- apply.monthly(X.d, FUN = colSums) # monthly


## b)

##' @title Manual Implementation of the Ljung--Box Test
##' @param x univariate data
##' @param lag maximal lag considered
##' @return p-value of the Ljung--Box test
##' @author Marius Hofert
ljung_box <- function(x, lag)
{
    rho.n <- acf(x, lag.max = lag, plot = FALSE)$acf[-1,,] # autocorrelation estimates for lags 1 to 'lag'
    n <- length(x)
    T.n <- n * (n+2) * sum(rho.n^2 / (n - seq_len(lag))) # test statistic; asymptotically chi_lag^2 under the null
    pchisq(T.n, df = lag, lower.tail = FALSE) # (asymptotic) p-value
}

## Compute p-values of the Ljung--Box test for each stock manually
mlag <- 10
p.d   <- apply(X.d,      2, ljung_box, lag = mlag)
p.d.a <- apply(abs(X.d), 2, ljung_box, lag = mlag)
p.m   <- apply(X.m,      2, ljung_box, lag = mlag)
p.m.a <- apply(abs(X.m), 2, ljung_box, lag = mlag)

## Compute p-values of the Ljung--Box test for each stock with Box.test()
LBtest <- function(x) Box.test(x, lag = mlag, type = "Ljung-Box")[["p.value"]] # auxiliary function
p.d.   <- apply(X.d,      2, LBtest)
p.d.a. <- apply(abs(X.d), 2, LBtest)
p.m.   <- apply(X.m,      2, LBtest)
p.m.a. <- apply(abs(X.m), 2, LBtest)

## Compare results
stopifnot(all.equal(p.d, p.d.), all.equal(p.d.a, p.d.a.),
          all.equal(p.m, p.m.), all.equal(p.m.a, p.m.a.))

## Which p-values are 'small'?
plot(p.d, ylim = 0:1, xlab = "Stock",
     ylab = paste0("p-values of Ljung-Box tests with maximal lag ",mlag))
points(p.d.a, col = "royalblue3")
points(p.m, pch = 4)
points(p.m.a, pch = 4, col = "royalblue3")
abline(h = 0.05, lty = 2)
legend("topleft", bty = "n",
       pch = c(1, 1, 4, 4), col = rep(c("black", "royalblue3"), 2),
       legend = c("Daily returns", "Daily absolute returns",
                  "Monthly returns", "Monthly absolute returns"))
## Note:
## We see that more 'blue' points are found below 0.05 (p-values of absolute
## returns) and overall 'many' which indicates that absolute values (and thus
## the returns) are not iid.


### Exercise 3.10 (Jarque--Bera tests) #########################################

## a)

## We use the NA-cleaned data from 3.9 here, too.

## Log-returns
X.d <- returns(S) # daily (as before)
X.m <- apply.monthly  (X.d, FUN = colSums) # monthly (as before)
X.q <- apply.quarterly(X.d, FUN = colSums) # quarterly


## b)

##' @title Manual Implementation of the Jarque--Bera Test
##' @param x univariate data
##' @return p-value of the Jarque--Bera test
##' @author Marius Hofert
jarque_bera <- function(x)
{
    x <- as.numeric(x)
    n <- length(x)
    skew.n <- skewness(x)
    kurt.n <- kurtosis(x)
    T.n <- (n/6) * (skew.n^2 + ((kurt.n - 3)^2) / 4) # test statistic; asymptotically chi_2^2 under the null
    pchisq(T.n, df = 2, lower.tail = FALSE) # p-value
}

## Compute p-values of the Jarque--Bera test for each stock manually
p.d <- apply(X.d, 2, jarque_bera)
p.m <- apply(X.m, 2, jarque_bera)
p.q <- apply(X.q, 2, jarque_bera)

## Compute p-values of the Jarque--Bera test for each stock with jarque.bera()
JBtest <- function(x) jarque.test(x)$p.value
p.d. <- apply(X.d, 2, JBtest)
p.m. <- apply(X.m, 2, JBtest)
p.q. <- apply(X.q, 2, JBtest)

## Compare results
stopifnot(all.equal(p.d, p.d.), all.equal(p.m, p.m.), all.equal(p.q, p.q.))

## Which p-values are 'small'?
plot(p.d, ylim = 0:1, xlab = "Stock", ylab = "p-values of Jarque-Bera tests",
     col = "maroon3")
points(p.m, col = "royalblue3")
points(p.q)
abline(h = 0.05, lty = 2)
legend("topleft", bty = "n", pch = rep(1, 3),
       col = c("maroon3", "royalblue3", "black"),
       legend = c("Daily returns", "Monthly returns", "Quarterly returns"))
## Note:
## For daily/monthly/quarterly log-returns, the hypothesis of normality is
## always/sometimes/less often rejected. We thus see that the larger the time
## period, the less often the Jarque--Bera test of normality rejects (due to a
## Central Limit Theorem effect).


### Exercise 3.12 (Changing correlation over time) #############################

## a), b)

## Data preparation
data(SP500)
data(DAX)
S <- merge(SP500, DAX)['2011-01-03/2015-12-31']
numNA <- rowSums(is.na(S))
any(numNA > 1) # => at most one data point missing
sum(numNA == 1) # => 46 dates with missing data points
S <- na.fill(S, fill = "extend") # fill NAs (by extension)
colnames(S) <- c("S&P 500", "DAX")
X <- returns(S)

##' @title Estimating Cross-Correlations Over Time
##' @param X bivariate xts object
##' @param window.size window size for estimating cross-correlations over time
##' @param ... additional arguments passed to the underlying cor()
##' @return xts object containing the estimated cross-correlations over time
##' @author Marius Hofert
##' @note - The default cor() call is equivalent to (but faster than)
##'         as.numeric(ccf(as.ts(X[t,1]), as.ts(X[t,2]), plot = FALSE, lag.max = 0)$acf)
##'       - apply.monthly(X, FUN = cor)[,2] gives monthly non-overlapping
##'         cross-correlation estimates (but the overlapping estimates show
##'         the time-varying effect equally well)
CCF <- function(X, window.size, ...)
{
    n <- nrow(X)
    stopifnot(window.size <= n, ncol(X) == 2)
    n.res <- n - window.size + 1
    res <- numeric(n.res)
    for(i in seq_len(n.res)) {
        t <- i:(i+window.size-1) # cross-correlation estimated for data from index i to i + window.size - 1
        res[i] <- cor(X[t,1], X[t,2], ...)
    }
    xts(res, order.by = index(X)[window.size:n]) # date is end point of each window
}

## Compute cross-correlations over time
window.size <- c(5, 25, 125)
rho.n.lst <- lapply(window.size, function(ws) CCF(X, window.size = ws))
rho.n <- do.call(merge, rho.n.lst)
names(rho.n) <- paste0("WS.",window.size)

## Plot
time <- index(rho.n)
cols <- c("maroon3", "black", "royalblue3")
plot(rep(NA, length(time)) ~ time, type = "l", ylim = c(-1, 1),
     xlab = "Time", ylab = "S&P500-DAX rolling correlation estimates")
for(j in seq_len(ncol(rho.n)))
    lines(time, rho.n[,j], col = cols[j])
legend("bottomleft", bty = "n", lty = rep(1, 3), col = cols,
       legend = paste0("Window size ",window.size))
## Note:
## The plot changes drastically with different window sizes, being too erratic for
## the smallest window size and perhaps oversmoothed for the largest. In all cases
## we see that the cross-correlation varies in time, which confirms (M3).


### Exercise 3.14 (Sensitivity of correlation estimates to volatility) #########

## a)

## Block the data in blocks of size 25 and estimate non-overlapping correlations
## and volatilities
n <- nrow(X)
ws <- 25 # window size
X.blk <- split(X, f = rep(seq_len(ceiling(n/ws)), each = ws, length.out = n))
rho.n. <- sapply(X.blk, function(x) cor(x)[1,2]) # non-overlapping correlation estimates
vol.n. <- t(sapply(X.blk, function(x) apply(x, 2, sd))) # non-overlapping volatility estimates; 2-column matrix

## Plot
x.S <- vol.n.[,"S&P 500"] # volatility estimates of S&P 500
x.D <- vol.n.[,"DAX"] # volatility estimates of DAX
y <- atanh(rho.n.) # Fisher (z-)transformed correlation estimates
plot(x.S, y, xlim = range(x.S, x.D),
     xlab = "Estimated stock volatility", ylab = "Fisher-transformed correlation estimates")
points(x.D, y, col = "royalblue3")
legend("bottomright", bty = "n", pch = c(1, NA, 1, NA), lty = c(NA, 1, NA, 1),
       col = rep(c("black", "royalblue3"), each = 2),
       legend = c("S&P 500", "Regression line S&P 500", "DAX", "Regression line DAX"))
mtext(substitute("All estimates are computed from non-overlapping windows of size"~ws.,
      list(ws. = ws)), side = 4, adj = 0, line = 0.5)
## Add regression lines (of Fisher transformed non-overlapping correlations on the
## volatility estimates)
abline(lm(y ~ x.S))
abline(lm(y ~ x.D), col = "royalblue3")
## Note:
## As in MFE (2015, Section 3.2.2), we see that the regression lines are
## rather meaningless here. Again, this is not an argument against the view
## that correlations are higher when volatilites are higher, it is just difficult
## to show with our estimated correlations and volatilities in this ad hoc way.


## b)

##' @title Estimating Volatilities Over Time
##' @param X multivariate xts object
##' @param window.size window size for estimating volatilities over time
##' @return xts object containing the estimated volatilities over time
##' @author Marius Hofert
VF <- function(X, window.size)
{
    n <- nrow(X)
    stopifnot(window.size <= n)
    n.res <- n - window.size + 1
    res <- matrix(, nrow = n.res, ncol = ncol(X))
    colnames(res) <- colnames(X)
    for(i in seq_len(n.res)) {
        t <- i:(i+window.size-1) # cross-correlation estimated for data from index i to i + window.size - 1
        res[i,] <- apply(X[t,], 2, sd)
    }
    xts(res, order.by = index(X)[window.size:n]) # date is end point of each window
}

## Rolling estimates (for window size 25)
rho.n.. <- rho.n.lst[[2]] # from previous exercise
vol.n.. <- VF(X, window.size = ws) # marginal rolling volatility estimates

## Plot
x.S. <- as.numeric(vol.n..[,"S&P 500"]) # volatility estimates of S&P 500
x.D. <- as.numeric(vol.n..[,"DAX"]) # volatility estimates of DAX
y. <- atanh(rho.n..) # Fisher (z-)transformed correlation estimates
plot(x.S., y., xlim = range(x.S., x.D.),
     xlab = "Estimated stock volatility", ylab = "Fisher-transformed correlation estimates")
points(x.D., y., col = "royalblue3")
legend("bottomright", bty = "n", pch = c(1, NA, 1, NA), lty = c(NA, 1, NA, 1),
       col = rep(c("black", "royalblue3"), each = 2),
       legend = c("S&P 500", "Regression line S&P 500", "DAX", "Regression line DAX"))
mtext(substitute("All estimates are computed from rolling windows of size"~ws.,
      list(ws. = ws)), side = 4, adj = 0, line = 0.5)
abline(lm(y. ~ x.S.))
abline(lm(y. ~ x.D.), col = "royalblue3")
## Note: The same interpretation as in a) applies.


### Exercise 3.15 (Geometric spacings between largest values in iid samples) ###

## d)

## Setup
p <- c(0.01, 0.02, 0.04)
k <- 1:200
true <- sapply(p, function(p.) pgeom(k-1, prob = p.))
appr <- sapply(p, function(p.) pexp(k, rate = p.))

## Plot of F_Geo(p) and F_Exp(p)
if(doPDF) pdf(file = (file <- "fig_03_geometric_spacings_F.pdf"))
plot(k, true[,1], ylim = 0:1, type = "l",
     ylab = expression(F[Geo(p)](k)~"(solid) and"~F[Exp(p)](k)~"(dashed)"))
lines(k, appr[,1], lty = 2)
lines(k, true[,2], col = "royalblue3")
lines(k, appr[,2], col = "royalblue3", lty = 2)
lines(k, true[,3], col = "maroon3")
lines(k, appr[,3], col = "maroon3", lty = 2)
legend("bottomright", bty = "n", lty = rep(1, 3),
       col = c("black", "royalblue3", "maroon3"),
       legend = as.expression(lapply(1:3, function(i)
           substitute(p==p., list(p. = p[i])))))
if(doPDF) dev.off.crop(file)

## Plot of the relative error
if(doPDF) pdf(file = (file <- "fig_03_geometric_spacings_rel_error.pdf"))
rel.err <- abs((true-appr)/true)
plot(k, rel.err[,1], type = "l", ylim = range(rel.err),
     ylab = "Relative error at k")
lines(k, rel.err[,2], col = "royalblue3")
lines(k, rel.err[,3], col = "maroon3")
legend("topright", lty = rep(1, 3), bty = "n",
       col = c("black", "royalblue3", "maroon3"),
       legend = as.expression(lapply(1:3, function(i)
           substitute(p==p., list(p. = p[i])))))
if(doPDF) dev.off.crop(file)
