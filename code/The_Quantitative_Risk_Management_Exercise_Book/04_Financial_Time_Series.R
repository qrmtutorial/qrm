## By Marius Hofert

## R script for Chapter 4 of The QRM Exercise Book


### Setup ######################################################################

library(xts)
library(rugarch)
library(qrmdata)
library(qrmtools)
doPDF <- require(crop)
options(digits = 10)


### Exercise 4.2 (Identification of ARMA processes) ############################

## Reproducing code

## Note:
## - ARMA(p1, q1)-GARCH(p2, q2) in 'rugarch' is:
##   X_t = mu_t + eps_t for eps_t = sig_t * Z_t and
##   mu_t = mu + sum_{k=1}^{p1} phi_k (X_{t-k} - mu) + sum_{k=1}^{q1} theta_k eps_{t-k},
##   sig_t^2 = alpha_0 + sum_{k=1}^{p2} alpha_k eps_{t-k}^2 + sum_{k=1}^{q2} beta_k sig_{t-k}^2
## - rugarch's phi_k is ar<k>, theta_k is ma<k>, omega is alpha_0

## a) MA(2) with theta_1 = -0.3 and theta_2 = -0.4
meanModel <- list(armaOrder = c(0, 2))
varModel <- list(garchOrder = c(0, 0))
pars <- list(mu = 0, ma1 = -0.3, ma2 = -0.4, omega = 1) # need omega = alpha_0 = 1 (see above)
uspec <- ugarchspec(mean.model = meanModel, variance.model = varModel,
                    fixed.pars = pars)
n <- 500 # sample size
pth <- ugarchpath(uspec, n.sim = n, n.start = 1, rseed = 271)
X.MA2 <- pth@path$seriesSim # length n

## b) Non-stationary process
set.seed(271)
eps <- rnorm(n)
t <- 1:n
X.nons <- 10 + 20 * t + 2 * t * sin(t) + eps

## c) AR(2) with phi_1 = 0.3 and phi_2 = 0.4
meanModel <- list(armaOrder = c(2, 0))
varModel <- list(garchOrder = c(0, 0))
pars <- list(mu = 0, ar1 = 0.3, ar2 = 0.4, omega = 1) # need omega = alpha_0 = 1 (see above)
uspec <- ugarchspec(mean.model = meanModel, variance.model = varModel,
                    fixed.pars = pars)
pth <- ugarchpath(uspec, n.sim = n, n.start = 1, rseed = 271)
X.AR2 <- pth@path$seriesSim # length n

## d) MA(1) process with theta_1 = -0.7
meanModel <- list(armaOrder = c(0, 1))
varModel <- list(garchOrder = c(0, 0))
pars <- list(mu = 0, ma1 = -0.7, omega = 1) # need omega = alpha_0 = 1 (see above)
uspec <- ugarchspec(mean.model = meanModel, variance.model = varModel,
                    fixed.pars = pars, distribution.model = "norm")
n <- 500 # sample size
pth <- ugarchpath(uspec, n.sim = n, n.start = 1, rseed = 271)
X.MA1 <- pth@path$seriesSim # length n

## Check sign of theta_1 (ma1)
if(FALSE) {
    acf(X.MA1, type = "covariance", plot = FALSE) # should be theta_1 at lag 1 => correct
    uspec. <- ugarchspec(mean.model = meanModel, # variance.model = varModel,
                         distribution.model = "norm")
    ugarchfit(uspec., data = X.MA1, solver = "hybrid", solver.control = list(trace = 1))@fit$solver$sol$pars
}

## e) ARMA(1,1) with phi_1 = 0.3 and theta_1 = -0.6
meanModel <- list(armaOrder = c(1, 1))
varModel <- list(garchOrder = c(0, 0))
pars <- list(mu = 0, ma1 = 0.3, ar1 = -0.6, omega = 1) # need omega = alpha_0 = 1 (see above)
uspec <- ugarchspec(mean.model = meanModel, variance.model = varModel,
                    fixed.pars = pars)
pth <- ugarchpath(uspec, n.sim = n, n.start = 1, rseed = 271)
X.ARMA11 <- pth@path$seriesSim # length n

## Plot
if(doPDF) pdf(file = (file <- "fig_04_ARMA_order.pdf"), width = 8, height = 10)
opar <- par(no.readonly=TRUE) # save plot parameters
lay <- matrix(1:18, ncol = 3, byrow = TRUE) # layout matrix
layout(lay, heights = c(0.2, rep(1, 5))) # layout
if(FALSE)
    layout.show(18) # display the layout
## Column headings
par(mar = c(0, 2.4, 0, 0)) # (bottom, left, top, right)
plot.new()
text(0.5, 0.5, labels = "Time series", font = 2, cex = 1.5)
plot.new()
text(0.5, 0.5, labels = "ACF", font = 2, cex = 1.5)
plot.new()
text(0.5, 0.5, labels = "PACF", font = 2, cex = 1.5)
par(mar = c(3, 3, 0, 0.5)) # set margins
## X.MA1
plot(X.MA1, type = "l", main = "", xlab = "", ylab = "") # time series
acf (X.MA1,             main = "", xlab = "", ylab = "", ci.col = "black") # ACF
pacf(X.MA1,             main = "", xlab = "", ylab = "", ci.col = "black") # PACF
## X.nons
plot(X.nons, type = "l", main = "", xlab = "", ylab = "") # time series
acf (X.nons,             main = "", xlab = "", ylab = "", ci.col = "black") # ACF
pacf(X.nons,             main = "", xlab = "", ylab = "", ci.col = "black") # PACF
## X.AR2
plot(X.AR2, type = "l", main = "", xlab = "", ylab = "") # time series
acf (X.AR2,             main = "", xlab = "", ylab = "", ci.col = "black") # ACF
pacf(X.AR2,             main = "", xlab = "", ylab = "", ci.col = "black") # PACF
## X.MA2
plot(X.MA2, type = "l", main = "", xlab = "", ylab = "") # time series
acf (X.MA2,             main = "", xlab = "", ylab = "", ci.col = "black") # ACF
pacf(X.MA2,             main = "", xlab = "", ylab = "", ci.col = "black") # PACF
## X.ARMA11
plot(X.ARMA11, type = "l", main = "", xlab = "", ylab = "") # time series
acf (X.ARMA11,             main = "", xlab = "", ylab = "", ci.col = "black") # ACF
pacf(X.ARMA11,             main = "", xlab = "", ylab = "", ci.col = "black") # PACF
## Finalize
par(opar)
if(doPDF) dev.off.crop(file)
## Note:
## X1: ACF cuts off after leg 1, PACF decays exponentially => MA(1) as in d)
## X2: Non-stationary => b)
## X3: ACF decays exponentially, PACF cuts off after leg 2 => AR(2) as in c)
## X4/5: Not easy to distinguish (X4 comes from a), X5 from e), but seeing it
##       the other way around is equally fine)


### Exercise 4.16 (Fitting GARCH models to equity index return data) ###########

## Data preparation
data(SP500)
S <- SP500['2006-01-01/2019-12-31'] # pick out the data we work with
X.d <- -returns(S)


## a)

##' @title Scatterplot and Q-Q Plot of the Residuals for GARCH(1,1) Model Fit Assessment
##' @param x data
##' @param armaOrder ARMA (p,q) order
##' @param distribution.model innovation distribution (standardized)
##' @return fitted ARMA-GARCH(1,1) object (invisibly)
##' @author Marius Hofert
GARCH11_check <- function(x, armaOrder, distribution.model = c("norm", "std", "sstd"))
{
    ## Check
    distribution.model <- match.arg(distribution.model)
    ## Model fit
    uspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                        mean.model = list(armaOrder = armaOrder, include.mean = TRUE),
                        distribution.model = distribution.model)
    fit <- ugarchfit(spec = uspec, data = x) # fit
    resi <- as.numeric(residuals(fit, standardize = TRUE)) # standardized residuals (equivalent to fit@fit$z)
    ## Quantile function of fitted innovation distribution
    switch(distribution.model,
           "norm" = {
               qF <- function(p) qdist("norm", p = p)
               lab <- c(distr = "Standard normal", info = NA)
           },
           "std" = {
               nu <- coef(fit)[["shape"]] # fitted degrees of freedom
               qF <- function(p) qdist("std", p = p, shape = nu)
               ## Note: As ./src/distributions.c -> qstd() reveals, this t distribution
               ##       is already standardized to have variance 1 (so t_nu(0, (nu-2)/nu))
               ##       when sigma = 1 (the default).
               lab <- c(distr = "Standardized t",
                        info = paste("Fitted degrees of freedom:", round(nu, 2)))
           },
           "sstd" = {
               sk <- coef(fit)[["skew"]] # fitted skewness parameter
               nu <- coef(fit)[["shape"]] # fitted degrees of freedom
               qF <- function(p) qdist("sstd", p = p, skew = sk, shape = nu)
               ## Note: Also this (skewed t) distribution is already standardized.
               ##       This can also be quickly verified via:
               ##       set.seed(271); var(rdist("sstd", n = 1e6, skew = 3.5, shape = 5.5))
               lab <- c(distr = "Standardized skewed t",
                        info = paste0("Fitted skewness: ",round(sk, 2),
                                      "; fitted degrees of freedom: ",round(nu, 2)))
           },
           stop("Wrong 'distribution.model'"))
    ## Plots
    layout(t(1:2))
    plot(resi, ylab = paste("Residuals of GARCH(1,1) with",
                            tolower(lab["distr"]), "innovations")) # Scatterplot
    mtext(lab["info"], side = 4, adj = 0, line = 0.5)
    qq_plot(resi, FUN = qF, xlab = paste(lab["distr"], "quantiles")) # Q-Q plot
    mtext(lab["info"], side = 4, adj = 0, line = 0.5)
    layout(1)
    invisible(fit)
}

## Pure GARCH(1,1)
GARCH11_check(X.d, armaOrder = c(0, 0), distribution.model = "norm") # standard normal innovations
GARCH11_check(X.d, armaOrder = c(0, 0), distribution.model = "std") # standardized t innovations
GARCH11_check(X.d, armaOrder = c(0, 0), distribution.model = "sstd") # standardized skewed t innovations

## By a small margin, the best seem to be the standardized skewed t distribution as
## innovation distribution.


## b)

## AR(1)-GARCH(1,1)
GARCH11_check(X.d, armaOrder = c(1, 0), distribution.model = "norm") # standard normal innovations
GARCH11_check(X.d, armaOrder = c(1, 0), distribution.model = "std") # standardized t innovations
GARCH11_check(X.d, armaOrder = c(1, 0), distribution.model = "sstd") # standardized skewed t innovations

## MA(1)-GARCH(1,1)
GARCH11_check(X.d, armaOrder = c(0, 1), distribution.model = "norm") # standard normal innovations
GARCH11_check(X.d, armaOrder = c(0, 1), distribution.model = "std") # standardized t innovations
GARCH11_check(X.d, armaOrder = c(0, 1), distribution.model = "sstd") # standardized skewed t innovations

## ARMA(1,1)-GARCH(1,1)
GARCH11_check(X.d, armaOrder = c(1, 1), distribution.model = "norm") # standard normal innovations
GARCH11_check(X.d, armaOrder = c(1, 1), distribution.model = "std") # standardized t innovations
GARCH11_check(X.d, armaOrder = c(1, 1), distribution.model = "sstd") # standardized skewed t innovations

## Overall, we see that the resulting fit of the models cannot be improved by adding
## a low-order ARMA model for the conditional mean.


## c)

## We work with the GARCH(1,1) model with skewed t innovations.
fit <- GARCH11_check(X.d, armaOrder = c(0, 0), distribution.model = "sstd") # get fitted model
vola <- sigma(fit) # = fit@fit$sigma = sqrt(fit@fit$var); in-sample (sigma_t) estimate
alpha <- 0.99 # confidence level
mu <- fitted(fit) # fitted conditional mean (here: constant, as armaOrder = c(0, 0))
VaR <- mu + vola * qdist("sstd", p = alpha, # in-sample VaR_alpha estimate
                         skew = coef(fit)[["skew"]], shape = coef(fit)[["shape"]]) # manually computed
VaR. <- quantile(fit, probs = alpha) # via rugarch()'s quantile method
stopifnot(all.equal(VaR, VaR., check.attributes = FALSE)) # comparison


## d)

## EWMA estimate of volatility (see also MFE (2015, Section 4.2.2))
uspec.EWMA <- ugarchspec(variance.model = list(model = "iGARCH", garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))
fit <- ugarchfit(spec = uspec.EWMA, data = X.d) # fit
vola.EWMA <- as.numeric(sigma(fit)) # extract vola; in-sample (sigma_t)

## Comparison with the GARCH volatility estimate
time <- index(vola)
plot(time, as.numeric(vola), type = "l",
     xlab = "Time t", ylab = expression("Estimated volatility"~(sigma[t])))
lines(time, vola.EWMA, col = "royalblue3")
legend("topleft", bty = "n", lty = rep(1, 2), col = c("black", "royalblue3"),
       legend = c("GARCH(1,1) volatility estimate", "EWMA estimate of volatility"))

## Relative difference
plot(time, abs(vola - vola.EWMA)/abs(vola), type = "l", xlab = "Time t",
     ylab = "Relative difference between the two volatility estimates")
## Note: |(y-x)/y| = eps => |y-x| = eps*|y| => (+/-)(y-x) = (+/-)*eps*y
##       => y-x = (+/-)*eps*y => x = (1 -/+ eps) * y, so x approximates y
##       up to eps * 100%.

## Summary statistics of the relative differencea
summary(as.numeric(abs(vola - vola.EWMA)/abs(vola)))
## => The relative difference between the two estimated volatilities
##    is typically around 3% but can reach up to 10%.


### Exercise 4.17 (Fitting GARCH models to weekly equity index return data) ####

## Data preparation
data(SP500)
S <- SP500['2000-01-01/2015-12-31'] # pick out the data we work with
X.w <- apply.weekly(-returns(S), FUN = colSums) # weekly
nrow(X.d)
nrow(X.w) # a bit less data


## a)

## Pure GARCH(1,1)
GARCH11_check(X.w, armaOrder = c(0, 0), distribution.model = "norm") # standard normal innovations
GARCH11_check(X.w, armaOrder = c(0, 0), distribution.model = "std") # standardized t innovations
GARCH11_check(X.w, armaOrder = c(0, 0), distribution.model = "sstd") # standardized skewed t innovations

## The best seem to be the standardized skewed t distribution as innovation
## distribution (as before). The fitted degrees of freedom are a bit larger
## (hinting at a Central Limit Theorem effect).


## b)

## AR(1)-GARCH(1,1)
GARCH11_check(X.w, armaOrder = c(1, 0), distribution.model = "norm") # standard normal innovations
GARCH11_check(X.w, armaOrder = c(1, 0), distribution.model = "std") # standardized t innovations
GARCH11_check(X.w, armaOrder = c(1, 0), distribution.model = "sstd") # standardized skewed t innovations

## MA(1)-GARCH(1,1)
GARCH11_check(X.w, armaOrder = c(0, 1), distribution.model = "norm") # standard normal innovations
GARCH11_check(X.w, armaOrder = c(0, 1), distribution.model = "std") # standardized t innovations
GARCH11_check(X.w, armaOrder = c(0, 1), distribution.model = "sstd") # standardized skewed t innovations

## ARMA(1,1)-GARCH(1,1)
GARCH11_check(X.w, armaOrder = c(1, 1), distribution.model = "norm") # standard normal innovations
GARCH11_check(X.w, armaOrder = c(1, 1), distribution.model = "std") # standardized t innovations
GARCH11_check(X.w, armaOrder = c(1, 1), distribution.model = "sstd") # standardized skewed t innovations

## Overall, we see that the resulting fit of the models cannot be improved by adding
## a low-order ARMA model for the conditional mean (as before).


## c)

## We work with the GARCH(1,1) model with skewed t innovations.
fit <- GARCH11_check(X.w, armaOrder = c(0, 0), distribution.model = "sstd") # get fitted model
vola <- sigma(fit) # = fit@fit$sigma = sqrt(fit@fit$var); in-sample (sigma_t) estimate
alpha <- 0.99 # confidence level
mu <- fitted(fit) # fitted conditional mean (here: constant, as armaOrder = c(0, 0))
VaR <- mu + vola * qdist("sstd", p = alpha, # in-sample VaR_alpha estimate
                         skew = coef(fit)[["skew"]], shape = coef(fit)[["shape"]])


## d)

## EWMA estimate of volatility (see also MFE (2015, Section 4.2.2))
fit <- ugarchfit(spec = uspec.EWMA, data = X.w) # fit
vola.EWMA <- as.numeric(sigma(fit)) # extract vola; in-sample (sigma_t)

## Comparison with the GARCH volatility estimate
time <- index(vola)
plot(time, as.numeric(vola), type = "l",
     xlab = "Time t", ylab = expression("Estimated volatility"~(sigma[t])))
lines(time, vola.EWMA, col = "royalblue3")
legend("topleft", bty = "n", lty = rep(1, 2), col = c("black", "royalblue3"),
       legend = c("GARCH(1,1) volatility estimate", "EWMA estimate of volatility"))

## Relative difference
plot(time, abs(vola - vola.EWMA)/abs(vola), type = "l", xlab = "Time t",
     ylab = "Relative difference between the two volatility estimates")

## Summary statistics of the relative differencea
summary(as.numeric(abs(vola - vola.EWMA)/abs(vola)))
## => The relative difference between the two estimated volatilities
##    is typically around 7% but can reach up to 36%.


### Exercise 4.18 (Fitting GARCH models to foreign-exchange return data) #######

## Data preparation
data(EUR_USD)
S <- EUR_USD['2012-01-01/2015-12-31'] # pick out the data we work with
X <- apply.weekly(returns(S), FUN = colSums)


## a)

## Pure GARCH(1,1)
GARCH11_check(X, armaOrder = c(0, 0), distribution.model = "norm") # standard normal innovations
GARCH11_check(X, armaOrder = c(0, 0), distribution.model = "std") # standardized t innovations
GARCH11_check(X, armaOrder = c(0, 0), distribution.model = "sstd") # standardized skewed t innovations

## There is barely any visible difference between the models and so we opt for
## the standard normal innovation distribution.


## b)

## AR(1)-GARCH(1,1)
GARCH11_check(X, armaOrder = c(1, 0), distribution.model = "norm") # standard normal innovations
GARCH11_check(X, armaOrder = c(1, 0), distribution.model = "std") # standardized t innovations
GARCH11_check(X, armaOrder = c(1, 0), distribution.model = "sstd") # standardized skewed t innovations

## MA(1)-GARCH(1,1)
GARCH11_check(X, armaOrder = c(0, 1), distribution.model = "norm") # standard normal innovations
GARCH11_check(X, armaOrder = c(0, 1), distribution.model = "std") # standardized t innovations
GARCH11_check(X, armaOrder = c(0, 1), distribution.model = "sstd") # standardized skewed t innovations

## ARMA(1,1)-GARCH(1,1)
GARCH11_check(X, armaOrder = c(1, 1), distribution.model = "norm") # standard normal innovations
GARCH11_check(X, armaOrder = c(1, 1), distribution.model = "std") # standardized t innovations
GARCH11_check(X, armaOrder = c(1, 1), distribution.model = "sstd") # standardized skewed t innovations

## Overall, we see that the resulting fit of the models cannot be improved by adding
## a low-order ARMA model for the conditional mean (as before).


## c)

## We work with the GARCH(1,1) model with standard normal innovations.
fit <- GARCH11_check(X, armaOrder = c(0, 0), distribution.model = "norm") # get fitted model
vola <- sigma(fit) # = fit@fit$sigma = sqrt(fit@fit$var); in-sample (sigma_t) estimate
alpha <- 0.99 # confidence level
mu <- fitted(fit) # fitted conditional mean (here: constant, as armaOrder = c(0, 0))
VaR <- mu + vola * qdist("norm", p = alpha) # in-sample VaR_alpha estimate


## d)

## EWMA estimate of volatility (see also MFE (2015, Section 4.2.2))
fit <- ugarchfit(spec = uspec.EWMA, data = X) # fit
vola.EWMA <- as.numeric(sigma(fit)) # extract vola; in-sample (sigma_t)

## Comparison with the GARCH volatility estimate
time <- index(vola)
plot(time, as.numeric(vola), type = "l",
     xlab = "Time t", ylab = expression("Estimated volatility"~(sigma[t])))
lines(time, vola.EWMA, col = "royalblue3")
legend("topleft", bty = "n", lty = rep(1, 2), col = c("black", "royalblue3"),
       legend = c("GARCH(1,1) volatility estimate", "EWMA estimate of volatility"))

## Relative difference
plot(time, abs(vola - vola.EWMA)/abs(vola), type = "l", xlab = "Time t",
     ylab = "Relative difference between the two volatility estimates")

## Summary statistics of the relative differencea
summary(as.numeric(abs(vola - vola.EWMA)/abs(vola)))
## => The relative difference between the two estimated volatilities
##    is typically around 2% but can reach up to 6%.
