## by Alexander McNeil

library(xts)
library(qrmdata)
library(rmgarch)

# Load some real data
data("FTSE")
data("SMI")

INDEXES <- merge(FTSE,SMI,all=FALSE)
plot.zoo(INDEXES)
# compute returns
FTSE.X <- diff(log(FTSE))[-1]
SMI.X <- diff(log(SMI))[-1]
INDEXES.X <- merge(FTSE.X,SMI.X,all=FALSE)
plot.zoo(INDEXES.X)

# take 4 years of data
data <- INDEXES.X['2006-01-01/2009-12-31']
pairs(as.zoo(data))
dim(data)

## Specify univariate AR(1)-GARCH(1,1) for both marginal processes
uspec.norm <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                       mean.model=list(armaOrder=c(1,0),include.mean=TRUE),
                    distribution.model="norm")
uspec.std <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                         mean.model=list(armaOrder=c(1,0),include.mean=TRUE),
                         distribution.model="std")

# Check the univariate specification for the two component series
fit.marg1 <- ugarchfit(spec=uspec.std,data=data[,1])
fit.marg2 <- ugarchfit(spec=uspec.std,data=data[,2])

# Combine univariate specs to obtain spec for marginal models
marginspec.norm <- multispec(replicate(2,uspec.norm))
marginspec.std <- multispec(replicate(2,uspec.std))

# Create spec for DCC
mspec <- dccspec(marginspec.std, dccOrder = c(1,1), model="DCC", distribution="mvt")

mod <- dccfit(mspec,data)
mod


# Some pictures of fit
plot(mod,which=2)
plot(mod,which=3)
plot(mod,which=4)
plot(mod,which=5)


# Check marginal coefficients are same in joint model
class(mod)
class(fit.marg1)
getSlots("DCCfit")
modfit <- mod@mfit
names(modfit)
modfit$coef
getSlots("uGARCHfit")
marg1fit <- fit.marg1@fit
marg1fit$coef


# A model with a changing copula

copspec = cgarchspec(uspec = marginspec.std, 
                     distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, transformation = "parametric"))
mod2 <- cgarchfit(copspec,data)
mod2

# Compare with DCC model

likelihood(mod2)
likelihood(mod)

cbind(coef(mod),coef(mod2))
