## By Alexander McNeil and Marius Hofert

## Comparison of (some) mean excess plots

library(QRM)

set.seed(271)
n <- 5000

MEplot(abs(rnorm(n))) # mean excess plot of |N(0,1)| data
MEplot(abs(rt(n, df = 4))) # mean excess plot of |t_4| data
MEplot(rexp(n)) # mean excess plot of Exp(1) data; should theoretically be a constant!
MEplot(runif(n)) # mean excess plot of U(0,1) data
