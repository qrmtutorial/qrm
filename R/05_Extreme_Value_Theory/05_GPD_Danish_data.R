## By Alexander McNeil and Marius Hofert


### Setup ######################################################################

library(xts)
library(qrmtools)
data(fire)


### 1 Find the threshold(s) ####################################################

## Plot
plot.zoo(fire, ylab = "Danish fire insurance claim losses in 1M DKK")

## Sample mean excess plot (with two plausible thresholds)
mean_excess_plot(fire) # sample mean excess plot; omits largest three points ('omit = 3')
u10 <- 10 # threshold 1
u20 <- 20 # threshold 2
abline(v = c(u10, u20))

## Effect of changing the threshold on xi
GPD_shape_plot(fire)
abline(v = c(u10, u20))
abline(h = 0.5) # models above line have infinite variance (E(X^k) = Inf <=> xi >= 1/k; k = 2 here)


### 2 Fit GPDs to the excesses #################################################

## Fit GPD models to excesses via MLE
exceed.u10 <- fire[fire > u10] # exceedances
exceed.u20 <- fire[fire > u20] # exceedances
excess.u10 <- exceed.u10 - u10 # excesses
excess.u20 <- exceed.u20 - u20 # excesses
(fit.u10 <- fit_GPD_MLE(excess.u10)) # MLE
shape.u10 <- fit.u10$par[["shape"]]
scale.u10 <- fit.u10$par[["scale"]]
(fit.u20 <- fit_GPD_MLE(excess.u20)) # MLE
shape.u20 <- fit.u20$par[["shape"]]
scale.u20 <- fit.u20$par[["scale"]]


### 3 Visually compare the sample excess df and the fitted (GPD) excess df #####

df.u10 <- function(q) # define fitted GPD
    pGPD(q, shape = shape.u10, scale = scale.u10) # fitted GPD
df.u20 <- function(q) # define fitted GPD
    pGPD(q, shape = shape.u20, scale = scale.u20) # fitted GPD

## Plot empiricial excess df vs fitted GPD
## u10
res <- edf_plot(excess.u10, log = "x")
z <- tail(res$t, n = -1)
lines(z, pGPD(z, shape = shape.u10, scale = scale.u10)) # fitted GPD
## u20
res <- edf_plot(excess.u20, log = "x")
z <- tail(res$t, n = -1)
lines(z, pGPD(z, shape = shape.u20, scale = scale.u20)) # fitted GPD

## Plot empiricial exceedance df vs shifted fitted GPD
## u10
res <- edf_plot(exceed.u10, log = "x")
z <- tail(res$t, n = -1)
lines(z, pGPD(z-u10, shape = shape.u10, scale = scale.u10)) # shifted fitted GPD
## u20
res <- edf_plot(exceed.u20, log = "x")
z <- tail(res$t, n = -1)
lines(z, pGPD(z-u20, shape = shape.u20, scale = scale.u20)) # shifted fitted GPD

## Corresponding Q-Q plots (more meaningful)
qf.u10 <- function(p) # quantile function of df
    qGPD(p, shape = shape.u10, scale = scale.u10)
qf.u20 <- function(p) # quantile function of df
    qGPD(p, shape = shape.u20, scale = scale.u20)
qq_plot(excess.u10, FUN = qf.u10)
qq_plot(excess.u20, FUN = qf.u20)
qq_plot(exceed.u10, FUN = function(p) u10 + qf.u10(p))
qq_plot(exceed.u20, FUN = function(p) u20 + qf.u20(p))


### 4 Compute semi-parametric risk measure estimators ##########################

## VaR_alpha, ES_alpha for two alphas and both thresholds
alpha <- c(0.99, 0.995)
(VaR.u10 <- VaR_POT(alpha, threshold = u10, p.exceed = mean(fire > u10),
                    shape = shape.u10, scale = scale.u10))
(VaR.u20 <- VaR_POT(alpha, threshold = u20, p.exceed = mean(fire > u20),
                    shape = shape.u20, scale = scale.u20))
(ES.u10 <- ES_POT(alpha, threshold = u10, p.exceed = mean(fire > u10),
                  shape = shape.u10, scale = scale.u10))
(ES.u20 <- ES_POT(alpha, threshold = u20, p.exceed = mean(fire > u20),
                  shape = shape.u20, scale = scale.u20))


### 5 Semi-parametric Smith estimator including VaR_0.99 #######################

## Empirical tail probabilities with Smith estimator overlaid
## u = 10
tail_plot(fire, threshold = u10, shape = shape.u10, scale = scale.u10)
abline(h = 1-alpha[1], v = VaR.u10[1], lty = 2) # 0.99, VaR_0.99
abline(h = 1-alpha[1], v = ES.u10[1],  lty = 2) # 0.99, ES_0.99
## u = 20
tail_plot(fire, threshold = u20, shape = shape.u20, scale = scale.u20)
abline(h = 1-alpha[1], v = VaR.u20[1], lty = 2) # 0.99, VaR_0.99
abline(h = 1-alpha[1], v = ES.u20[1],  lty = 2) # 0.99, ES_0.99

## A version including confidence intervals for VaR_0.99 and ES_0.99
## u = 10
opar <- par(mar = c(5, 4, 4, 2) + c(0, 1, 0, 0))
fit.u10. <- QRM::fit.GPD(fire, u10)
QRM::showRM(fit.u10., alpha = 0.99, RM = "VaR", method = "BFGS") # with VaR estimate and CIs
QRM::showRM(fit.u10., alpha = 0.99, RM = "ES",  method = "BFGS") # with ES  estimate and CIs
par(opar)
## u = 20
opar <- par(mar = c(5, 4, 4, 2) + c(0, 1, 0, 0))
fit.u20. <- QRM::fit.GPD(fire, u20)
QRM::showRM(fit.u20., alpha = 0.99, RM = "VaR", method = "BFGS") # with VaR estimate and CIs
QRM::showRM(fit.u20., alpha = 0.99, RM = "ES",  method = "BFGS") # with ES  estimate and CIs
par(opar)
## => CIs quite a bit wider (as we have less data)
