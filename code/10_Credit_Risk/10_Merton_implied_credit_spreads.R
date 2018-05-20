### by Alexander McNeil

## functions

MertonSpread <- function(leverage,tau,sigmaV)
{
dt1 <- (-log(leverage) + 0.5*tau*sigmaV^2)/(sigmaV*sqrt(tau))
dt2 <-	dt1-sigmaV*sqrt(tau)
-log(pnorm(dt2)+pnorm(-dt1)/leverage)/tau
}


## Spread versus asset volatility; leverage and time-to-maturity fixed
leverage = 0.6
tau  = 2
sigmaV = seq(from=0.01, to = 0.5, length=50)
cta = MertonSpread(leverage,tau,sigmaV)
plot(sigmaV,cta,type="l")

## Spread versus time-to-maturity; leverage and volatility fixed
leverage = 0.6
sigmaV = 0.25
tau = seq(from = 0.01, to = 5, length = 50)
ctb = MertonSpread(leverage,tau,sigmaV)
plot(tau,ctb,type="l")

## Spread versus leverage; volatility and time-to-maturity fixed
leverage = seq(from=0.01, to = 0.99, length = 50)
# low asset volatility
sigmaV = 0.25
tau = 2
ctc = MertonSpread(leverage,tau,sigmaV)
plot(leverage,ctc,type="l")
# high asset volatility
sigmaV = 0.5
tau = 2
ctc = MertonSpread(leverage,tau,sigmaV)
plot(leverage,ctc,type="l")
# high time-to-maturity
sigmaV = 0.5
tau = 4
ctc = MertonSpread(leverage,tau,sigmaV)
plot(leverage,ctc,type="l")


