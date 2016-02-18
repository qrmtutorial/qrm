delta.t = 0.25
N = 20
r = 0.02
xstar = 0.0042
delta = 0.6

premium.leg = function(gamma,xstar,N,delta.t,r)
{
	k = 1:N
	xstar*sum(exp(-r*k*delta.t)*delta.t*exp(-gamma*k*delta.t))
}

gamma = seq(from = 0, to = 0.02, length = 20)
pleg=rep(NA,length(gamma))
for (i in 1:length(gamma))
 	pleg[i] = premium.leg(gamma[i], xstar, N, delta.t, r)


default.leg = function(gamma,delta,N,delta.t,r)
{
	delta*gamma*(1-exp(-(r+gamma)*N*delta.t))/(r+gamma)
}

dleg=rep(NA,length(gamma))
for (i in 1:length(gamma))
 	dleg[i] = default.leg(gamma[i], delta, N, delta.t, r)

plot(gamma,pleg,ylim=range(dleg,pleg),type="n")
lines(gamma,dleg,col=2)
lines(gamma,pleg,col=3)


rootfunc = function(gamma,xstar,N,delta.t,r,delta)
{
	premium.leg(gamma,xstar,N,delta.t,r)-default.leg(gamma,delta,N,delta.t,r)
}
(tmp=uniroot(rootfunc,interval=c(0,0.02),xstar=xstar,N=N,delta.t=delta.t,r=r,delta=delta))
abline(v=tmp$root)

# approximate value
abline(v=xstar/delta,lty=2,col=4)

# simpler solution possible
(tmp2=uniroot(rootfunc,interval=c(0,0.02),xstar=xstar,N=1,delta.t=delta.t,r=r,delta=delta))
