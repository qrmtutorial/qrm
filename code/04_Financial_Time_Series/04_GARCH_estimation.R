## By Alexander McNeil and Marius Hofert


### Setup ######################################################################

library(rugarch)
library(zoo)
library(ADGofTest) # for ad.test()
library(moments) # for skewness(), kurtosis()
library(qrmdata)
library(qrmtools)


### 1 S&P 500 data #############################################################

## Load S&P 500 data
data("SP500")
plot.zoo(SP500, xlab = "Time")

## Extract the data we work with and build log-returns
SPdata <- SP500['2006-01-01/2009-12-31'] # 4 years of data
X <- returns(SPdata)
plot.zoo(X, xlab = "Time")


### 2 Fit an AR(1)--GARCH(1,1) model with normal innovations ###################

## Model specification (without fixed.pars, so without keeping any parameter
## fixed during the optimization)
uspec.N <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), # AR(1) part
                                        include.mean = TRUE), # with mean
                      distribution.model = "norm") # normal innovations
(fit.N <- ugarchfit(spec = uspec.N, data = X))
## Note: Fit to 'X', not 'SPdata'!

## The fit contains a lot of information
## - The parameter estimates are the "Optimal Parameters" (at the top)
## - plot(fit.N) will take us into a menu system (below are the most important ones)

## Series with conditional quantiles
plot(fit.N, which = 2)
layout(matrix(1:4, ncol = 2, byrow = TRUE)) # specify (2,2)-matrix of plots
plot(fit.N, which = 6) # ACF of absolute data |X_t| (shows serial correlation)
plot(fit.N, which = 9) # Q-Q plot of standardized residuals Z_t (shows leptokurtosis of Z_t; normal assumption not supported)
plot(fit.N, which = 10) # ACF of standardized residuals Z_t (shows AR dynamics do a reasonable job of explaining conditional mean)
plot(fit.N, which = 11) # ACF of squared standardized residuals Z_t^2 (shows GARCH dynamics do a reasonable job of explaining conditional sd)
layout(1) # restore layout


### 3 Fit an AR(1)--GARCH(1,1) model with Student t innovations ################

## Now consider t innovations
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                      distribution.model = "std") # Student t innovations
fit.t <- ugarchfit(spec = uspec.t, data = X)

## The pictures are similar, but the Q-Q plot looks "better"
layout(matrix(1:4, ncol = 2, byrow = TRUE)) # specify (2,2)-matrix of plots
plot(fit.t, which = 6) # ACF of |X_t|
plot(fit.t, which = 9) # Q-Q plot of Z_t against a normal
plot(fit.t, which = 10) # ACF of Z_t
plot(fit.t, which = 11) # ACF of Z_t^2
layout(1)


### 4 Likelihood-ratio test for comparing the two models #######################

## The normal model is "nested within" the t model (appears from the t as
## special (limiting) case). The likelihood-ratio test tests:
## H0: normal innovations are good enough
## HA: t innovations necessary
## Decision: We reject the null (in favour of the alternative) if the
##           likelihood-ratio test statistic exceeds the 0.95 quantile of a
##           chi-squared distribution with 1 degree of freedom (1 here as that's
##           the difference in the number of parameters for the two models)

## The likelihood ratio statistic is twice the difference between the log-likelihoods
LL.N <- fit.N@fit$LLH # log-likelihood of the model based on normal innovations
LL.t <- fit.t@fit$LLH # log-likelihood of the model based on t innovations
LRT <- 2*(fit.t@fit$LLH-fit.N@fit$LLH) # likelihood-ratio test statistic
LRT > qchisq(0.95, 1) # => H0 is rejected
1-pchisq(LRT, df = 1) # p-value (probability of such an extreme result if normal hypothesis were true)


### 5 Using Akaike's information criterion to compare the two models ###########

## Akaike's information criterion (AIC) can be used for model comparison.
## Favour the model with smallest AIC.
## Note: AIC = 2*<number of estimated parameters> - 2 * log-likelihood
(AIC.N <- 2*length(fit.N@fit$coef) - 2*LL.N)
(AIC.t <- 2*length(fit.t@fit$coef) - 2*LL.t)
## => t model is preferred


### 6 Exploring the structure of the fitted model object #######################

## Basic components
class(fit.t)
getClass("uGARCHfit")
getSlots("uGARCHfit")
fit.t.fit <- fit.t@fit
names(fit.t.fit)
fit.t.model <- fit.t@model
names(fit.t.model)
fit.t.model$modeldesc
fit.t.model$pars

## To find out about extraction methods, see ?ugarchfit
## Then follow a link to documentation for "uGARCHfit" object

## Some extraction methods
(param <- coef(fit.t)) # estimated coefficients
sig <- sigma(fit.t) # estimated volatility
VaR.99 <- quantile(fit.t, probs = 0.99) # estimated VaR at level 99%
Z <- residuals(fit.t, standardize = TRUE) # estimated standardized residuals Z_t

## Plots
plot.zoo(sig,    xlab = "Time", ylab = expression(hat(sigma)[t]))
plot.zoo(VaR.99, xlab = "Time", ylab = expression(widehat(VaR)[0.99]))
plot.zoo(Z,      xlab = "Time", ylab = expression(hat(Z)[t]))

## More residual checks
(mu.Z <- mean(Z)) # ok (should be ~= 0)
(sd.Z <- as.vector(sd(Z))) # ok (should be ~= 1)
skewness(Z) # should be 0 (is < 0 => left skewed)
hist(Z) # => left skewed
kurtosis(Z) # should be 6/(nu-4)
nu <- param[["shape"]] # estimated degrees-of-freedom
6/(nu-4) # => sample kurtosis larger than it should be
pt.hat <- function(q) pt((q-mu.Z)/sd.Z, df = nu) # estimated t distribution for Z
## Z ~ t_{hat(nu)}(hat(mu), hat(sig)^2) -> F_Z(z) = t_{hat(nu)}((z-hat(mu))/hat(sig))
ad.test(as.numeric(Z), distr.fun = pt.hat) # AD test
## => Anderson--Darling test rejects the estimated t distribution

## Possible things to try:
## 1) A GARCH model with asymmetric dynamics - GJR-GARCH
## 2) A GARCH model with asymmetric innovation distribution
## 3) Different ARMA specifications for the conditional mean


### 7 Forecasting/predicting from the fitted model #############################

(forc <- ugarchforecast(fit.t, n.ahead = 10))
plot(forc, which = 1) # prediction of X_t
plot(forc, which = 3) # prediction of sigma_t
