## By Alexander J. McNeil

library(xts)
library(qrmdata)
library(rugarch)
library(qrmdata)


data(SP500)
data(VIX)

SP500.r <- diff(log(SP500))['1990-01-03/']
VIX.r <- diff(log(VIX))['1990-01-03/']
length(SP500.r)
length(VIX.r)

X.d <- merge(SP500.r, VIX.r)
X.w <- apply.weekly(X.d, FUN = colSums)

## BACKTESTING SECTION


## Recall key risk-factor changes for option price
plot(as.matrix(X.d))
plot(as.matrix(X.w))
## Consider these as history of risk-factor changes up to present time.

## A bank is short a put option
## Put option characteristics
K <- 100
# Option currently "in of the money"
S <- 80
r <- 0.03
sigma <- 0.25
T <- 1


(Vt <- -Black_Scholes(0, S, r, sigma, K, T, "put"))
Vtplus1.d <- -Black_Scholes(1/250, exp(log(S)+X.d[,1]), r, exp(log(sigma)+X.d[,2]), K, T, "put")
Vtplus1.w <- -Black_Scholes(1/52, exp(log(S)+X.w[,1]), r, exp(log(sigma)+X.w[,2]), K, T, "put")
gains.d <- Vtplus1.d-Vt
gains.w <- Vtplus1.w-Vt
## The above time series are historically simulated "gains"
plot(gains.d)
qqnorm(gains.d)
plot(gains.w)
qqnorm(gains.w)

## These series bear all the signs of stochastic volatility
acf(gains.d)
acf(abs(gains.d))

spec.n <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
                     distribution.model = "norm")
spec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(1,1), include.mean = TRUE), distribution.model = "std")
spec.tskew <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                         mean.model = list(armaOrder = c(1,1), include.mean = TRUE), distribution.model = "sstd")


## We now carry out a series of rolling out-of-sample VaR predictions
## We refit the model every 10 days

roll.n <- ugarchroll(spec.n, gains.d, n.start = 1000, window.size = 1000, refit.every = 10,
                     refit.window = "moving", solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = c(0.05, 0.01))
show(roll.n)
plot(roll.n,which = 4)
report(roll.n, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)

## A little more detail on how VaR numbers are extracted
VaR.table.n <- as.data.frame(roll.n, which = "VaR")
names(VaR.table.n)
realized <- VaR.table.n[,"realized"]
VaR05.n = VaR.table.n[,1]
VaR01.n = VaR.table.n[,2]
VaRTest(alpha = 0.01, actual = realized, VaR = VaR01.n)
par.n <- as.data.frame(roll.n, which = "density")
names(par.n)
ES05.n = par.n$Mu - par.n$Sigma*ES_t(0.95, df = Inf)
ES01.n = par.n$Mu - par.n$Sigma*ES_t(0.99, df = Inf)
ESTest(alpha = 0.05, actual = realized, ES = ES05.n, VaR = VaR05.n)
ESTest(alpha = 0.01, actual = realized, ES = ES01.n, VaR = VaR01.n)


## Now try fitting the Student t Model

roll.t <- ugarchroll(spec.t, gains.d, n.start = 1000, window.size = 1000, refit.every = 10,
                     refit.window = "moving", solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = c(0.05, 0.01))
show(roll.t)
plot(roll.t,which = 4)
report(roll.t, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)


VaR.table.t <- as.data.frame(roll.t, which = "VaR")
names(VaR.table.t)
realized <- VaR.table.t[,"realized"]
VaR05.t = VaR.table.t[,1]
VaR01.t = VaR.table.t[,2]

par.t <- as.data.frame(roll.t, which = "density")
names(par.t)
ES05.t = par.t$Mu - par.t$Sigma*ESst(0.95, df = par.t$Shape, scale = TRUE)
ES01.t = par.t$Mu - par.t$Sigma*ESst(0.99, df = par.t$Shape, scale = TRUE)
ESTest(alpha = 0.05, actual = realized, ES = ES05.t, VaR = VaR05.t)
ESTest(alpha = 0.01, actual = realized, ES = ES01.t, VaR = VaR01.t)

## Now try skew-t innovations. SLOW

roll.tskew <- ugarchroll(spec.tskew, gains.d, n.start = 1000, window.size = 1000, refit.every = 10,
                         refit.window = "moving", solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = c(0.05, 0.01))
show(roll.tskew)
plot(roll.tskew, which = 4)
report(roll.tskew, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)


