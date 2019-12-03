## By Marius Hofert

## See MFE (2015, Example 2.16)


### Setup ######################################################################

library(qrmtools)


### 1 Auxiliary computations ###################################################

## Computing VaR and ES for a normal and t distribution for X_{t+1} (or the
## linearized loss derived from X_{t+1})
V <- 10000 # value of the portfolio today
sig <- 0.2/sqrt(250) # daily volatility (annualized volatility of 20%)
nu <- 4 # degrees of freedom for the t distribution
alpha <- 1-2^seq(log(1-0.001, base = 2), -10, length.out = 256) # confidence levels; concentrated near 1
VaR.n <- VaR_t(alpha, scale = V*sig, df = Inf) # VaR_alpha under normal
VaR.t <- VaR_t(alpha, scale = V*sig*sqrt((nu-2)/nu), df = nu) # VaR_alpha under t
ES.n  <- ES_t (alpha, scale = V*sig, df = Inf) # ES_alpha under normal
ES.t  <- ES_t (alpha, scale = V*sig*sqrt((nu-2)/nu), df = nu) # ES_alpha under t


### 2 Plot #####################################################################

plot(1-alpha, ES.t, type = "l", ylim = range(VaR.n, VaR.t, ES.n, ES.t), log = "x",
     col = "maroon3", xlab = expression(1-alpha), ylab = "")
lines(1-alpha, VaR.t, col = "darkorange2")
lines(1-alpha, ES.n,  col = "royalblue3")
lines(1-alpha, VaR.n, col = "black")
legend("topright", bty = "n", lty = rep(1,4), col = c("maroon3", "darkorange2",
                                                      "royalblue3", "black"),
       legend = c(substitute(ES[alpha]~~"for "*italic(t[nu.])*" model", list(nu. = nu)),
                substitute(VaR[alpha]~~"for "*italic(t[nu.])*" model", list(nu. = nu)),
                expression(ES[alpha]~~"for normal model"),
                expression(VaR[alpha]~~"for normal model")))
## Results:
## This shows that VaR_alpha (or ES_alpha) is not always 'riskier' for the
## t distribution than it is for the normal distribution (only for
## sufficiently large alpha)
