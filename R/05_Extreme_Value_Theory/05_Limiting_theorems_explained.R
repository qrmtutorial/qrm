## By Marius Hofert

## Visualization of major limiting theorems in Probability (SLLN, CLT) and
## Extreme Value Theory (Gnedenko's Theorem, Pickands--Balkema--de Haan)
## based on the *same* (simulated) data.


### Setup ######################################################################

library(qrmtools)

## Data from a Par(theta) distribution
n <- 2e5 # sample size = number of iid random variables
th <- 3 # parameter theta
set.seed(271) # set seed for reproducibility
X <- rPar(n, theta = th) # generate data


### 1 Strong Law of Large Numbers (SLLN) #######################################

## Building cumulative averages (X_1/1, (X_1+X_2)/2, (X_1+X_2+X_3)/3,...)
stopifnot(th > 1) # for mean to exist
Xn <- cumsum(X)/(1:n)

## Plot (this one path of the stochastic process (bar{X}_n)_{n=1}^{infty})
plot(1:n, Xn, type = "l", log = "x", ylab = "",
     xlab = expression("Number n of iid random variables"~(X[i])[i == 1]^n),
     main = substitute(bold("Strong Law of Large Numbers for Par("*th.*") data"), list(th. = th)))
mu <- 1/(th-1) # true mean
abline(h = mu, col = "royalblue3")
legend("bottomright", lty = c(1,1), col = c("black", "royalblue3"),
       bty = "n", y.intersp = 1.2,
       legend = c(expression((bar(X)[n])[n]),
                  substitute("true mean"~mu == mu., list(mu. = mu))))


### 2 Central Limit Theorem (CLT) ##############################################

## Build blocks of data
m <- 500 # number of blocks (each of size n/m = 400)
X. <- split(X, f = rep(1:m, each = floor(n/m))) # split data into blocks

## Location-scale transform blocked sums via sqrt(n) * (bar{X}_n - mu) / sigma
## = (S_n - n * mu) / (sqrt(n) * sigma)
stopifnot(th > 2) # for variance to exist
mu <- 1/(th-1) # true mean
sig2 <- 2/((th-2)*(th-1)^2) # true variance
Z <- sapply(X., function(x) (sum(x) - length(x) * mu) / (sqrt(length(x) * sig2))) # standardize by mean(<sum>) and sd(<sum>)

## Histogram with overlaid densities
dens <- density(Z)
hist(Z, probability = TRUE, ylim = c(0, max(dnorm(0), dens$y)), breaks = 20,
     main = substitute(bold("Central Limit Theorem for Par("*th.*") data"),
                       list(th. = th)), xlab = expression("Realizations of"~sqrt(n)*(bar(X)[n]-mu)/sigma))
lines(dens, col = "royalblue3") # overlaid density estimate
curve(dnorm, from = min(Z), to = max(Z), add = TRUE, col = "darkorange2") # overlaid N(0,1) density
box()
legend("topright", lty = c(1,1), col = c("royalblue3", "darkorange2"), bty = "n",
       legend = c("Density estimate", "N(0,1) density"))

## Q-Q plot
qq_plot(Z, FUN = qnorm, main = "Central Limit Theorem")


### 3 Gnedenko's Theorem #######################################################

## bar{F}(x) = x^{-theta} L(x) for L(x) = (1+1/x)^{-theta} (slowly varying at
## infinity). By Gnedenko's Theorem (1943), such F is in MDA(H_{1/theta}) for
## all theta>0 and as normalizing sequences one can take c_n = F^-(1-1/n) and
## d_n = 0 for all n so that (M_n-d_n)/c_n is approximately distributed as
## H_{xi = 1/theta, mu, sig} for some xi, mu and sig. One could take estimates
## of xi, mu and sig, but one can show that besides xi = 1/theta, one has
## mu = 1 and sig = 1/theta here. Let's check that.

## Location-scale transform blocked maxima with c_n = F^-(1-1/n) and d_n = 0
M <- sapply(X., function(x) (max(x) - 0) / qPar(1-1/length(x), theta = th))

## Histogram with overlaid densities
dens <- density(M, adjust = 2) # smoothed density estimate
x <- seq(0, max(M), length.out = 257)
true.dens <- dGEV(x, xi = 1/th, mu = 1, sigma = 1/th) # true density
hist(M, probability = TRUE, ylim = c(0, max(dens$y, true.dens)), breaks = 60,
     main = substitute(bold("Gnedenko's Theorem for Par("*th.*") data"),
                       list(th. = th)), xlab = expression("Realizations of"~(M[n]-d[n])/c[n]~~
                                                          "(for"~c[n] == {F^{-1}}(1-1/n)*";"~~d[n] == 0*")"))
lines(dens, col = "royalblue3") # density estimate
lines(x, true.dens, col = "darkorange2")
box()
legend("topright", lty = c(1,1), col = c("royalblue3", "darkorange2"), bty = "n",
       legend = c("Density estimate", expression("Limiting GEV density"~h[list(xi,mu,sigma)])))

## Q-Q plot
qq_plot(M, FUN = function(p) qGEV(p, xi = 1/th, mu = 1, sigma = 1/th),
        main = substitute(bold("Gnedenko's Theorem for Par("*th.*") data"),
                          list(th. = th)))
## => For smaller block sizes (e.g. take n = 2e4 => n/m = 40), there is
##    significant departure visible (and even here there is departure visible).
##    This indicates that one typically needs quite a large (block and thus)
##    sample size to get a sufficiently good approximation to the limiting GEV.
##    The following approach is less 'wasteful' with the data and already
##    works well/better for smaller n.


### 4 Pickands--Balkema--de Haan (1974/1975) ###################################

## For sufficiently large thresholds u, excesses over u follow a GPD(xi, beta)
## distribution if and only if F is in MDA(H_xi). If F is GPD(xi, beta)
## one can show that the excess distribution F_u over u is GPD(xi, beta+xi*u).
## Since that's the case for Pareto distributions (with xi = 1/theta,
## beta = 1/theta), excesses over u should follow a GPD(1/theta, (1/theta)*(1+u))
## distribution. Let's check that.

## Determine excesses
u <- quantile(X, 0.9) # use the 90% quantile (rule of thumb)
Y <- X[X>u] - u # excesses over u

## Histogram (note: true density peaks near 0)
dens <- density(Y) # density estimate
x <- seq(0, max(Y), length.out = 257)
true.dens <- dGPD(x, xi = 1/th, beta = (1/th)*(1+u))
hist(Y, probability = TRUE, ylim = c(0, max(dens$y, true.dens)), breaks = 60,
     main = substitute(bold("Pickands-Balkema-de Haan Theorem for Par("*th.*") data"),
            list(th. = th)), xlab = "Realizations of excesses Y over the threshold u (90% quantile)")
lines(dens, col = "royalblue3") # density estimate
lines(x, true.dens, col = "darkorange2")
box()
legend("topright", lty = c(1,1), col = c("royalblue3", "darkorange2"), bty = "n",
       legend = c("Density estimate", expression("Limiting GPD density"~g[list(xi,beta(u))])))

## Just the density estimates on the log-scale
x <- 10^seq(-2, 2, length.out = 65)
true.dens <- dGPD(x, xi = 1/th, beta = (1/th)*(1+u))
ii <- dens$x > 0 # only use those values with density > 0 (otherwise log-scale fails)
plot(dens$x[ii], dens$y[ii], type = "l", log = "x", col = "royalblue3",
     ylim = c(0, max(dens$y[ii], true.dens)), xlab = "x", ylab = "Density")
lines(x, true.dens, col = "darkorange2")
legend("topright", lty = c(1,1), col = c("royalblue3", "darkorange2"), bty = "n",
       legend = c("Density estimate", expression("Limiting GPD density"~g[list(xi,beta(u))])))

## Q-Q plot
qq_plot(Y, FUN = function(p) qGPD(p, xi = 1/th, beta = (1/th)*(1+u)),
        main = substitute(bold("Pickands-Balkema-de Haan Theorem for Par("*th.*") data"),
                          list(th. = th)))
## => In particular, this looks better than the Q-Q plot before (this
##    'peaks-over-threshold' method is less wasteful with the data).