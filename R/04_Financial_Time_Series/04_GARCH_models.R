## by Alexander McNeil
require(rugarch)
require(QRM)
load("INDEXES-2000-2012.RData")


# Generate a GARCH model
GARCHspec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),mean.model=list(armaOrder=c(0,0),include.mean=FALSE),distribution.model="norm",fixed.pars=list(omega=0.02,alpha1=0.15,beta1=0.8))

# Generate a realization
path <- ugarchpath(GARCHspec,n.sim=2000,n.start=50,m.sim=1)
plot(path,which=1)
plot(path,which=2)
plot(path,which=3)
plot(path,which=4)

vol <- sigma(path)
plot(vol,type="h")
names(attributes(path))
names(path@path)
X <- path@path$seriesSim
acf(X)
acf(abs(X))
qqnorm(X)
qqline(X,col=2)
shapiro.test(X)

# In the next section we will analyse real data

X.INDEXES <- returns(INDEXES0012)
X.sp500 <- X.INDEXES[,2]

# Fit an AR(1)-GARCH(1,1) model with normal innovations
# We first have to create a model specification
AR.GARCH.N.spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),mean.model=list(armaOrder=c(1,0),include.mean=TRUE),distribution.model="norm")

fit.a <- ugarchfit(spec=AR.GARCH.N.spec,data=X.sp500)
# The fit contains a lot of information
# The parameter estimates are the "Optimal Parameters" at the top
fit.a
# plot(fit.a) will take us into a menu system
# We select the most important plots
# First the series with conditional quantiles
plot(fit.a,which=2)
par(mfrow=c(2,2))
# acf of absolute data - shows serial correlation
plot(fit.a,which=6)
# QQplot of data - shows leptokurtosis of standardized rediduals - normal assumption not supported
plot(fit.a,which=9)
# acf of standardized residuals - shows AR dynamics do a reasonable job of explaining conditional mean
plot(fit.a,which=10)
# acf of squared standardized residuals - shows GARCH dynamics do a reasonable job of explaining conditional sd
plot(fit.a,which=11)
par(mfrow=c(1,1))

# Fit an AR(1)-GARCH(1,1) model with student innovations
AR.GARCH.T.spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),mean.model=list(armaOrder=c(1,0),include.mean=TRUE),distribution.model="std")

fit.b <- ugarchfit(spec=AR.GARCH.T.spec,data=X.sp500)
fit.b
# The pictures are similar, but QQplot looks "better"
par(mfrow=c(2,2))
plot(fit.b,which=6)
plot(fit.b,which=9)
plot(fit.b,which=10)
plot(fit.b,which=11)
par(mfrow=c(1,1))

# Now compare log-likelihoods and Akaike numbers
fit.a@fit$LLH
fit.b@fit$LLH

# The normal model is "nested within" the t model
# normal is a special case
# The likelihood ratio statistic is twice the difference between the log-likelihoods
(LRT = 2*(fit.b@fit$LLH-fit.a@fit$LLH))
# It tests:
# H0: normal innovations are good enough
# HA: t innovations necessary
# We reject the null in favour of the alternative if LRT exceeds the 0.95 quantile of a chi-squared distribution with 1 degree of freedom
# The df is the difference in number of parameters
qchisq(0.95,1)
# The null is clearly rejected
# Alternatively we can give p-value (probability of such an extreme result if normal hypothesis were true)
1-pchisq(LRT,1)

# Akaike numbers can be used for nested or non-nested comparisons, particularly the latter
# Favour the model with smallest AIC
# in this case: normal model -6.2654; t model -6.2997

(forc <- ugarchforecast(fit.b, n.ahead=10))


