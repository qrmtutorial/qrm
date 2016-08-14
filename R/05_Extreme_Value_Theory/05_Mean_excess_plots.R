## By Alexander J. McNeil and Marius Hofert

## Comparison of (some) mean excess plots


set.seed(271)
n <- 5000

dat.abs.t4 <- abs(rt(n, df = 4)) # |t_4| data
MEplot(dat.abs.t4)

dat.exp <- rexp(n)
MEplot(dat.exp) # should theoretically be a constant!

dat.abs.n <- abs(rnorm(n))
MEplot(dat.abs.n)

dat.unif <- runif(n)
MEplot(dat.unif)
