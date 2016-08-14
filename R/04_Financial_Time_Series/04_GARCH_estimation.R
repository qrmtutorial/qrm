## By Alexander McNeil

library(rugarch)
library(qrmdata)
library(zoo)

## Load some real data
data("SP500")
plot.zoo(SP500)
## Compute log returns
SP500.r <- diff(log(SP500))[-1]
## Take 4 years of data
SP500.r <- SP500.r['2006-01-01/2009-12-31']


## Fit an AR(1)-GARCH(1,1) model with normal innovations
## We first have to create a model specification
AR.GARCH.N.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                              mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                              distribution.model = "norm")

fit.a <- ugarchfit(spec = AR.GARCH.N.spec, data = SP500.r)
## The fit contains a lot of information
## The parameter estimates are the "Optimal Parameters" at the top
fit.a
## plot(fit.a) will take us into a menu system
## We select the most important plots
## First the series with conditional quantiles
plot(fit.a, which = 2)
opar <- par(mfrow = c(2,2))
## acf of absolute data - shows serial correlation
plot(fit.a, which = 6)
## QQplot of data - shows leptokurtosis of standardized rediduals - normal assumption not supported
plot(fit.a, which = 9)
## acf of standardized residuals - shows AR dynamics do a reasonable job of explaining conditional mean
plot(fit.a, which = 10)
## acf of squared standardized residuals - shows GARCH dynamics do a reasonable job of explaining conditional sd
plot(fit.a, which = 11)
par(opar)

## Fit an AR(1)-GARCH(1,1) model with student innovations
AR.GARCH.T.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                              mean.model = list(armaOrder = c(1,0), include.mean=TRUE),
                              distribution.model = "std")

fit.b <- ugarchfit(spec = AR.GARCH.T.spec,data = SP500.r)
fit.b
## The pictures are similar, but QQplot looks "better"
opar <- par(mfrow = c(2,2))
plot(fit.b, which = 6)
plot(fit.b, which = 9)
plot(fit.b, which = 10)
plot(fit.b, which = 11)
par(opar)

## Now compare log-likelihoods and Akaike numbers
fit.a@fit$LLH
fit.b@fit$LLH

## The normal model is "nested within" the t model
## normal is a special case
## The likelihood ratio statistic is twice the difference between the log-likelihoods
(LRT <- 2*(fit.b@fit$LLH-fit.a@fit$LLH))
## It tests:
## H0: normal innovations are good enough
## HA: t innovations necessary
## We reject the null in favour of the alternative if LRT exceeds the 0.95 quantile of a chi-squared distribution with 1 degree of freedom
## The df is the difference in number of parameters
qchisq(0.95, 1)
## The null is clearly rejected
## Alternatively we can give p-value (probability of such an extreme result if normal hypothesis were true)
1-pchisq(LRT, 1)

## Akaike numbers can be used for nested or non-nested comparisons, particularly the latter
## Favour the model with smallest AIC
## in this case: normal model -6.2654; t model -6.2997

(forc <- ugarchforecast(fit.b, n.ahead = 10))

## In the next few commands we will explore the structure of the fitted model object
class(fit.a)
getClass("uGARCHfit")
getSlots("uGARCHfit")
fit.a.fit <- fit.a@fit
names(fit.a.fit)
fit.a.model <- fit.a@model
names(fit.a.model)
fit.a.model$modeldesc
fit.a.model$pars

## To find out about extraction methods:
?ugarchfit
## Then follow a link to documentation for "uGARCHfit" object
## Here are some common extractions:
coef(fit.a)
plot.zoo(sigma(fit.a))
plot.zoo(quantile(fit.a, 0.99))
Z.hat <- residuals(fit.a, standardize = TRUE)


## We have a look in more detail at the residuals
plot.zoo(Z.hat)
## Note how they are left skewed
hist(Z.hat)
mean(Z.hat);var(Z.hat)
shapiro.test(as.numeric(Z.hat))

## For these commands you will need the package "moments"
library(moments)
skewness(Z.hat)
kurtosis(Z.hat)
jarque.test(as.numeric(Z.hat))

## Things to try:
## A GARCH model with asymmetric dynamics - GJR-GARCH
## A GARCH model with asymmetric innovation distribution
## Different ARMA specifications for the conditional mean

