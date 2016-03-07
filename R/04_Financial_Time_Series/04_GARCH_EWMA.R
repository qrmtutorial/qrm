## by Alexander McNeil
## original idea by Alexios Ghalanos
library(rugarch)
library(qrmdata)



# Load some real data
data("SP500")
plot(SP500)
# Compute log returns
SP500.r <- diff(log(SP500))[-1]
# take 4 years of data
SP500.r <- SP500.r['2006-01-01/2009-12-31']

# Specify models
GARCH.spec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                             variance.model=list(model="sGARCH",  garchOrder=c(1,1)),
                             distribution.model="std")
EWMAfixed.spec = ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                             variance.model=list(model="iGARCH"), 
                            fixed.pars=list(alpha1=1-0.94, omega=0))
EWMAest.spec = ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                           variance.model=list(model="iGARCH"), 
                          fixed.pars=list(omega=0))

# Fit models
# Note second model has nothing to estimate so filter applied instead
mod1 = ugarchfit(GARCH.spec, SP500.r)
mod2 = ugarchfilter(EWMAfixed.spec, SP500.r)
mod3 = ugarchfit(EWMAest.spec, SP500.r)

# plot volatility estimates
plot(sigma(mod1), main="", auto.grid = FALSE, major.ticks = "auto",
     minor.ticks = FALSE)
lines(sigma(mod2), col=2, lty=2)
lines(sigma(mod3), col=3, lty=3)

# add a legend
l1 = as.expression("GARCH")
l2 = as.expression(substitute(paste("EWMA[fix.",lambda,"=",x,"]"), list(x=0.94)))
cf1=round(coef(mod3)[3],2)
l3 = as.expression(substitute(paste("EWMA[est.",lambda,"=",x,"]"), list(x=cf1)))
legend("topleft", c(l1, l2, l3), col=1:3, lty=1:3, bty="n", cex=0.9)



