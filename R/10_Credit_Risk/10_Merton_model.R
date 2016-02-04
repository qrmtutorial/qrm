### by Alexander McNeil

## required packages
library(qrmtools)
library(sde)


## Parameters for Merton model
V0 <- 1
muV <- 0.03
sigmaV <- 0.25
r <- 0.02
B <- 0.85
T <- 1
N <- 364

## Simulated asset value trajectories for Merton model
npaths <- 50
paths <- matrix(NA,nrow=N+1,ncol=npaths)
for (i in 1:npaths)
{paths[,i] <- GBM(x=V0,r=muV,sigma=sigmaV,T=T,N=N)}
 head(paths)

## plot paths

def.par <- par(no.readonly = TRUE)

times <- (1:(N+1))/(N+1)
par(mar=c(3,3,2,1),mgp=c(2,1,0))
plot(times,paths[,1],type="l",xlab="t",ylab=expression(V[t]),ylim=range(paths))
for (i in 1:npaths){
	lines(times,paths[,i],col=sample(1:20,1))
}
abline(v=1)
abline(h=B)

## Overlay a lognormal density
lnV.mu <- log(V0) + (muV-0.5*sigmaV^2)*T
lnV.sigma <- sigmaV*sqrt(T)
rvals <- range(paths)
vals <- seq(from=rvals[1], to=rvals[2],length=100)
dens <- dlnorm(vals,meanlog=lnV.mu,sdlog=lnV.sigma)
dvals <- (max(dens)-dens)/(max(dens)-min(dens))
lines(dvals,vals)

## Example of a default path also showing value of debt
set.seed(63)
Vt <- GBM(x=V0,r=muV,sigma=sigmaV,T=T,N=N)
times <- seq(from=0,to=1,length=N+1)
plot(times,Vt,type="l",ylim=range(0.6,max(Vt)),xlab="t",ylab=expression(V[t]))
abline(v=1)
# Add default free debt value
lines(times,B*exp(-r*(T-times)),col=3)
# Add defaultable debt
Bt <- B*exp(-r*(T-times)) - Black_Scholes(times,Vt,r,sigmaV,B,T,type="put")
lines(times,Bt,col=2)

## Example of a non-default path also showing value of debt
set.seed(77)
Vt <- GBM(x=V0,r=muV,sigma=sigmaV,T=T,N=N)
times <- seq(from=0,to=1,length=N+1)
plot(times,Vt,type="l",ylim=range(0.6,max(Vt)),xlab="t",ylab=expression(V[t]))
abline(v=1)
# Add default free debt value
lines(times,B*exp(-r*(T-times)),col=3)
# Add defaultable debt
Bt <- B*exp(-r*(T-times)) - Black_Scholes(times,Vt,r,sigmaV,B,T,type="put")
lines(times,Bt,col=2)

par(def.par)


