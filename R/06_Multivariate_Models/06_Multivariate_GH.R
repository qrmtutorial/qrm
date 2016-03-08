## By Alexander McNeil

## Libraries required
library(xts)
library(qrmdata)
library(ghyp)

## Dow Jones Data
data("DJ_const")
Sdata <- DJ_const['2000-01-01/',1:10]
Xdata <- diff(log(Sdata))[-1,]
Xdata.w <- apply.weekly(Xdata,FUN=colSums) # weekly data
plot.zoo(Xdata.w)

## Fit elliptical models
mod.NIG.sym <- fit.NIGmv(Xdata.w,symmetric=TRUE, nit = 10000)
summary(mod.NIG.sym)
mod.T.sym <- fit.tmv(Xdata.w,symmetric=TRUE, nit = 10000)
summary(mod.T.sym)
mod.GHYP.sym <- fit.ghypmv(Xdata.w,symmetric=TRUE, nit = 10000)
summary(mod.GHYP.sym)

## Fit skewed models
mod.NIG.skew <- fit.NIGmv(Xdata.w,symmetric=FALSE, nit = 10000)
summary(mod.NIG.skew)
mod.T.skew <- fit.tmv(Xdata.w,symmetric=FALSE, nit = 10000)
summary(mod.T.skew)
mod.GHYP.skew <- fit.ghypmv(Xdata.w,symmetric=FALSE, nit = 10000)
summary(mod.GHYP.skew)

mod.NIG.sym@llh
mod.NIG.skew@llh
mod.T.sym@llh
mod.T.skew@llh
mod.GHYP.sym@llh
mod.GHYP.skew@llh

## Hypothesis: elliptical T is good enough
LRstat <- 2*(mod.T.skew@llh-mod.T.sym@llh)
1-pchisq(LRstat,10) # => H0 not rejected

## Hypothesis: elliptical ghyp is good enough
LRstat <- 2*(mod.GHYP.skew@llh-mod.GHYP.sym@llh)
1-pchisq(LRstat,10) # => H0 not rejected

## Hypothesis: skewed t no worse than skewed ghyp
LRstat <- 2*(mod.GHYP.skew@llh-mod.T.skew@llh)
1-pchisq(LRstat,1) # => H0 rejected at 5% level

## Hypothesis: elliptical t no worse than elliptical ghyp
LRstat <- 2*(mod.GHYP.sym@llh-mod.T.sym@llh)
1-pchisq(LRstat,1) # => H0 rejected at 5% level

## Find skewness parameters
getSlots("mle.ghyp")
mod.GHYP.skew@gamma

## Make picture of two returns
Xdata <-  Xdata.w[,c(5,8)]
plot(as.numeric(Xdata[,1]),as.numeric(Xdata[,2]),xlab=names(Xdata)[1],ylab=names(Xdata)[2])
mod <- fit.ghypmv(Xdata,symmetric=FALSE, nit = 10000)
summary(mod)
xvals <- seq(from=min(Xdata[,1]),to=max(Xdata[,2]),length=50)
yvals <- seq(from=min(Xdata[,2]),to=max(Xdata[,2]),length=50)
outer.func <- function(x,y){dghyp(cbind(x,y),mod)}
zvals <- outer(xvals,yvals,outer.func)
zvals.un <- unique(as.numeric(zvals))
contour(xvals,yvals,zvals,add=TRUE,col=2,levels=seq(from=quantile(zvals.un,0.6),to=quantile(zvals.un,0.9),length=10))

