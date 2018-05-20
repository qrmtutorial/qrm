## By Alexander McNeil

## Original idea by Alexios Ghalanos


### Setup ######################################################################

library(rugarch)
library(zoo)
library(qrmdata)
library(qrmtools)


### 1 S&P 500 data #############################################################

## Load S&P 500 data
data("SP500")
plot.zoo(SP500, xlab = "Time")

## Extract the data we work with and build log-returns
SPdata <- SP500['2006-01-01/2009-12-31'] # 4 years of data
X <- returns(SPdata)


### 2 Model specifications and fitting #########################################

## Specify models
GARCH.spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                         variance.model = list(model = "sGARCH",  garchOrder = c(1,1)),
                         distribution.model = "std")
EWMAfixed.spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                             variance.model = list(model = "iGARCH"),
                             fixed.pars = list(alpha1 = 1-0.94, omega = 0))
EWMAest.spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                           variance.model = list(model = "iGARCH"),
                           fixed.pars = list(omega = 0))

## Model fitting
mod1 <- ugarchfit(GARCH.spec, X)
mod2 <- ugarchfilter(EWMAfixed.spec, X) # nothing to estimate, apply filter instead
mod3 <- ugarchfit(EWMAest.spec, X)


### 3 Plotting #################################################################

## Plot volatility estimates
plot(sigma(mod1), main = "", auto.grid = FALSE, major.ticks = "auto", minor.ticks = FALSE)
lines(sigma(mod2), col = 2, lty = 2)
lines(sigma(mod3), col = 3, lty = 3)
legend("topleft", bty = "n", col = 1:3, lty = 1:3, cex = 0.9,
       legend = as.expression(c("GARCH", substitute("EWMA[fix."*lambda == x*"]", list(x = 0.94)),
                                substitute("EWMA[est."*lambda == x*"]", list(x = round(coef(mod3)[3],2))))))
