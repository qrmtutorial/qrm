delta = 0.5
gamma = 0.04
T = 5
epsilon = 0.01
r = 0.02

RTspread = function(t,T,gamma,delta)
{
	-log(1- delta*(1-exp(-gamma*(T-t))))/(T-t)
}

t = seq(from=0,to=(T-epsilon),length=50)
c1 = RTspread(t,T,gamma,delta)

RFspread = function(t,T,gamma,delta,r)
{
	ratio1 = (1-delta)*gamma/(gamma+r)
	ratio2 = (r+delta*gamma)/(gamma+r)
	-r-log(ratio1+ratio2*exp(-(gamma+r)*(T-t)))/(T-t)
}

c2 = RFspread(t,T,gamma,delta,r)

plot(t,c1,ylim=range(0,c1,c2),type="n",ylab="spread")
abline(h=delta*gamma)
lines(t,c1,col=2)
lines(t,c2,col=3)
