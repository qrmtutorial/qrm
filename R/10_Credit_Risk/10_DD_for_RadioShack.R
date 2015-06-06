### By Alexander McNeil

require(timeSeries)
load("RadioShack.RData")
source("Black_Scholes.R")

# Use last couple of years of data

RadioShack = window(RadioShack, start="2010-01-01", end="2012-03-31")
# Number of shares in millions (approximately)
N.shares <- 100
Svalues <- RadioShack*N.shares
# Value of equity in millions of dollars (approx)
plot(Svalues)

Vvalues <- timeSeries(rep(NA,length(Svalues)),time(Svalues))
# Value of one-year debt in millions approximately
B <- 1042


# Root finding equation
rooteqn <- function(V,S,t,r,sigmaV,B,T)
{
	S - BlackScholes(t,V,r,sigmaV,B,T,"call")
}


# Initial estimate of volatility
sigmaV <- as.numeric(sd(returns(Svalues)))*sqrt(250)


# First iteration
for (i in 1:length(Svalues)){
	tmp <- uniroot(rooteqn, interval =c(Svalues[i],10*Svalues[i]),S=Svalues[i],t=0,r=0.03,sigmaV=sigmaV,B=B,T=1)
	Vvalues[i] <- tmp$root
}
sigmaV.old <- sigmaV
sigmaV <- as.numeric(sd(returns(Vvalues)))*sqrt(250)


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
sigmaV <- as.numeric(sd(returns(Vvalues)))*sqrt(250)
}

sigmaV
it
tail(Vvalues)
plot(Vvalues,ylim=range(Vvalues,Svalues))
lines(Svalues,col=2)

# Merton Distance to default
DD <- (log(Vvalues[length(Vvalues)])-log(B))/sigmaV

# Value given by Moody's
(log(1834)-log(1042))/0.24


