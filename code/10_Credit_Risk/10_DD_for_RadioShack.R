## By Alexander McNeil

library(xts)
library(qrmdata)
library(qrmtools)
data(RSHCQ)


# Use last couple of years of data

RadioShack = RSHCQ['2010-01-01/2012-03-31']
# Number of shares in millions (approximately)
N.shares <- 100
Svalues <- RadioShack*N.shares
# Value of equity in millions of dollars (approx)
plot(Svalues)

Vvalues <- xts(rep(NA,length(Svalues)),time(Svalues))
# Value of one-year debt in millions approximately
B <- 1042


# Root finding equation
rooteqn <- function(V,S,t,r,sigmaV,B,T)
{
	S - Black_Scholes(t,V,r,sigmaV,B,T,"call")
}


# Initial estimate of volatility
Svalues.X <- diff(log(Svalues))[-1]
sigmaV <- as.numeric(sd(Svalues.X))*sqrt(250)


# First iteration
for (i in 1:length(Svalues)){
	tmp <- uniroot(rooteqn, interval =c(Svalues[i],10*Svalues[i]),S=Svalues[i],t=0,r=0.03,sigmaV=sigmaV,B=B,T=1)
	Vvalues[i] <- tmp$root
}
sigmaV.old <- sigmaV
Vvalues.X <- diff(log(Vvalues))[-1]
sigmaV <- as.numeric(sd(Vvalues.X)*sqrt(250))
plot(Vvalues)


# Further iterations
it <- 1
while (abs(sigmaV-sigmaV.old)/sigmaV.old > 0.000001)
{
	it <- it + 1
for (i in 1:length(Svalues)){
	tmp <- uniroot(rooteqn, interval =c(Svalues[i],10*Svalues[i]),S=Svalues[i],t=0,r=0.03,sigmaV=sigmaV,B=B,T=1)
	Vvalues[i] <- tmp$root
}
sigmaV.old <- sigmaV
Vvalues.X <- diff(log(Vvalues))[-1]
sigmaV <- as.numeric(sd(Vvalues.X)*sqrt(250))
}

sigmaV
it
tail(Vvalues)
plot(Vvalues,ylim=range(Vvalues,Svalues))
lines(Svalues,col=2)

# Merton Distance to default
DD <- (log(Vvalues[length(Vvalues)])-log(B))/sigmaV
DD

# Value given by Moody's
(log(1834)-log(1042))/0.24


