### by Alexander J. McNeil

require(xts)
require(qrmdata)
require(termstrc)
data(ZCB_CA)

# NELSON-SIEGEL FACTOR MODEL OF THE YIELD CURVE

# We will start by taking a 100 day excerpt
head(ZCB_CA)
ZCBexcerpt <- ZCB_CA['2011-08-08/2011-12-30']
yields <- 100*as.matrix(ZCBexcerpt)
dates <- time(ZCBexcerpt)
maturities = (1:120)*0.25

(datazeroyields = zeroyields(maturities, yields, dates))

NSfit <- estim_nss(datazeroyields,"ns")
summary(NSfit)

# Now plot the Nelson-Siegel factors, after first renaming factors
factors <- NSfit$optparam
factors[,4] <- 1/factors[,4]
dimnames(factors)[[2]] <- c("Z1","Z2","Z3","eta")
params <- xts(factors,dates)
plot.zoo(params)

# Let's have a look at the quality of the fit on a particular day
day <- 1
fitted <- NSfit$yhat[day,]
thedate <- dates[day]

plot(maturities,yields[day,],main=as.character(thedate),ylim=range(yields[day],fitted),xlab="Time to maturity",ylab="Yield",pch=20)
lines(maturities,fitted)


k2 <- function(s,eta){
 (1-exp(-eta*s))/(eta*s)
}

k3 <- function(s,eta){
  k2(s,eta)-exp(-eta*s)
}
eta <- factors[,4]
k2.curve <- k2(maturities,eta[day])
k3.curve <- k3(maturities,eta[day])

par(oma=c(0,0,0,4))
plot(maturities,k2.curve,type="l",xlab="Time to Maturity",ylab="k2")
par(new=TRUE)
plot(maturities,k3.curve,type="l", ylab="",xlab="Time to Maturity",axes=FALSE)
axis(4)
mtext("k3",side=4,line=3)
par()

# Create a "zero-yields" object with 10 years of data now
ZCB10yr <- ZCB_CA['2002-01-02/2011-12-30']

yields <- 100*as.matrix(ZCB10yr)
dates <- time(ZCB10yr)
maturities = (1:120)*0.25
datazeroyields = zeroyields(maturities, yields, dates)
datazeroyields
# This time we hold eta fixed throughout the estimation (much quicker)
# This is the so-called Diebold-Li method (dl)
NSfit <- estim_nss(datazeroyields,method="dl",optimtype="global")


# Plot the Nelson-Siegel surface for different times and maturities
plot(NSfit)
# Now plot the Nelson-Siegel factors, after first renaming factors
factors <- NSfit$optparam
dimnames(factors)[[2]] <- paste("Z",1:3,sep="")
params <- xts(factors,dates)
plot.zoo(params)

# The next piece of code shows how we can mark a significant date on all 3 plots
plotfunc <- function(x, ...)
	{	
	lines(x,...)
	abline(v=as.Date("2008-09-15"),lty=2)
	}                                 
plot.zoo(params,panel=plotfunc,yax.flip=TRUE)

# These data are now the risk-factor changes
X <- diff(params)[-1,]
plot.zoo(X)
cor(X)




