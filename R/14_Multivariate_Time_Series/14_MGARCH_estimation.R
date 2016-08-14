## By Alexander McNeil

library(xts)
library(qrmdata)
library(rmgarch)

## Load some real data
data("FTSE")
data("SMI")

INDEXES <- merge(FTSE, SMI, all = FALSE)
plot.zoo(INDEXES)

## Compute returns
FTSE.X <- diff(log(FTSE))[-1]
SMI.X <- diff(log(SMI))[-1]
INDEXES.X <- merge(FTSE.X, SMI.X, all = FALSE)
plot.zoo(INDEXES.X)

## Take 4 years of data
data <- INDEXES.X['2006-01-01/2009-12-31']
pairs(as.zoo(data))
dim(data)

## Specify univariate AR(1)-GARCH(1,1) for both marginal processes
uspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                    mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                    distribution.model = "std")

## Check the univariate specification for the two component series
fit.marg1 <- ugarchfit(spec = uspec, data = data[,1])
fit.marg2 <- ugarchfit(spec = uspec, data = data[,2])

## Combine univariate specs to obtain spec for marginal models
marginspec <- multispec(replicate(2, uspec))

## Create spec for DCC
mspec <- dccspec(marginspec, dccOrder = c(1,1), model = "DCC", distribution = "mvt")

mod <- dccfit(mspec,data)
mod

## Check marginal coefficients are same in joint model
coef(mod)
coef(fit.marg1)
coef(fit.marg2)
## Note: a two stage fit is undertaken
## First univariate GARCH is fitted to each margin
## The standardized residuals are extracted
## Then a model with dynamically changing conditional correlation matrix is fitted
## (Ideally, all parameters should be estimated in one step.)

## Some pictures of fit
plot(mod, which = 2)
plot(mod, which = 3)
plot(mod, which = 4)
plot(mod, which = 5)


## A model with a changing copula
copspec <- cgarchspec(uspec = marginspec,
                      distribution.model = list(copula = "mvt", method = "ML",
                                                time.varying = TRUE, transformation = "parametric"))
mod2 <- cgarchfit(copspec, data)
mod2

## The only difference here is that a meta-t model is fitted
## with marginal parameters given by estimated shape parameters
## in first step

likelihood(mod2)
likelihood(mod)
## Fit is slightly superior

cbind(coef(mod), coef(mod2))
## Compare coefficients with DCC model
## Only last 3 are different