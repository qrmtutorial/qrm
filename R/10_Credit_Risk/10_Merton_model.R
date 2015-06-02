### by Alexander McNeil

## required packages
require(sde)
require(QRM)

### function
## to be moved to QRMtools

BlackScholes <- function(t, S, r, sigma, K, T, type=c("call","put"))
{
    d1 <- (log(S/K) + (r+sigma^2/2)*(T-t))/(sigma*sqrt(T-t))
    d2 <- d1-sigma*sqrt(T-t)
    type <- match.arg(type)
    switch(type,
           "call" = {
               S*pnorm(d1)-K*exp(-r*(T-t))*pnorm(d2)
           },
           "put" =  {
               S*(pnorm(d1)-1)+K*exp(-r*(T-t))*(1-pnorm(d2))
           },
           stop("Wrong type"))
}

## Parameters for Merton model
V0=1
muV=0.03
sigmaV = 0.25
r=0.02
B=0.85
T=1
N=364

## Simulated asset value trajectories for Merton model
npaths <- 50
paths <- matrix(NA,nrow=N+1,ncol=npaths)
for (i in 1:npaths)
{paths[,i] <- GBM(x=V0,r=muV,sigma=sigmaV,T=T,N=N)}
 head(paths)

## make timeSeries for pretty plotting
## Remove QRM dependence; make nicer plot
times <- timeSequence(from="2015-01-01",length=N+1,by="day")
paths <- timeSeries(paths,times)
plot(paths,plot.type="single")
abline(v=timeDate("2015-12-31"))
abline(h=B)

## Overlay a lognormal density
lnV.mu = log(V0) + (muV-0.5*sigmaV^2)*T
lnV.sigma = sigmaV*sqrt(T)
rvals = range(paths)
vals = seq(from=rvals[1], to=rvals[2],length=100)
dens = dlnorm(vals,meanlog=lnV.mu,sdlog=lnV.sigma)
par(new=TRUE)
plot(max(dens)-dens,vals,type="l",axes=FALSE,ann=FALSE)

## Example of a default path also showing value of debt
set.seed(63)
Vt = GBM(x=V0,r=muV,sigma=sigmaV,T=T,N=N)
times = seq(from=0,to=1,length=N+1)
plot(times,Vt,type="l",ylim=range(0.6,max(Vt)),xlab="t",ylab=expression(V[t]))
abline(v=1)
# Add default free debt value
lines(times,B*exp(-r*(T-times)),col=3)
# Add defaultable debt
Bt = B*exp(-r*(T-times)) - BlackScholes(times,Vt,r,sigmaV,B,T,type="put")
lines(times,Bt,col=2)

## Example of a non-default path also showing value of debt
set.seed(77)
Vt = GBM(x=V0,r=muV,sigma=sigmaV,T=T,N=N)
times = seq(from=0,to=1,length=N+1)
plot(times,Vt,type="l",ylim=range(0.6,max(Vt)),xlab="t",ylab=expression(V[t]))
abline(v=1)
# Add default free debt value
lines(times,B*exp(-r*(T-times)),col=3)
# Add defaultable debt
Bt = B*exp(-r*(T-times)) - BlackScholes(times,Vt,r,sigmaV,B,T,type="put")
lines(times,Bt,col=2)

