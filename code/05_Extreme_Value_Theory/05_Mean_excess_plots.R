## By Alexander McNeil and Marius Hofert

## Comparison of some mean excess plots

## Note: Often helpful for computing the mean excess function e(u) of a df F with
##       finite first moment is the formula
##
##           e(u) = (1/bar(F)(u)) int_u^{x_F} bar(F)(x) dx.                  (*)


### Setup ######################################################################

library(qrmtools)

n <- 50000 # sample size
## Note: n has to be large if you want to compare with the theoretical ME functions


### 1 F in the Frechet domain of attraction ####################################

### 1.1 GPD(xi, beta) ##########################################################

## We know that for F = GPD(xi, beta), e(u) = (beta + xi * u) / (1 - xi)
## for all u such that beta + xi * u > 0. Let's check!
xi <- 0.4
beta <- 3
stopifnot(xi < 1, beta > 0) # for E|X| < Inf
set.seed(271)
x <- rGPD(n, shape = xi, scale = beta)
mean_excess_plot(x)
abline(a = beta/(1-xi), b = xi/(1-xi), col = "royalblue3") # intercept, slope
## => Okay, but seems to look worse for smaller n (even n = 5000) and larger xi

## Note in particular that for xi > 0, the slope is positive.


### 1.2 Par(th, 1) #############################################################

## Since Par(th, 1) = GPD(1/th, 1/th), we should have e(u) = (1 + u) / (th - 1)
th <- 3
stopifnot(th > 1) # for E|X| < Inf
set.seed(271)
x <- rPar(n, shape = 3)
mean_excess_plot(x)
abline(a = 1/(th-1), b = 1/(th-1), col = "royalblue3") # intercept, slope
## => Okay


### 1.3 t_nu on the positive real line #########################################

## t_nu is also in the MDA of a Frechet, with xi = 1/nu (which can be shown
## based on the t_nu density and Karamata's Theorem; see MFE (2015, Theorem A.7))
nu <- 3
stopifnot(nu > 1) # for E|X| < Inf
set.seed(271)
x <- rt(n, df = nu)
mean_excess_plot(x) # => should maybe also omit the *smallest* three points (here)
## Formula (*) and the substitution v = t_nu(x) provide us with e(u)
## (with numerical integration on a compact interval)
f <- function(v) (1-v) / dt(qt(v, df = nu), df = nu)
mef <- function(u) integrate(f, lower = pt(u, df = nu), upper = 1)$value / pt(u, df = nu, lower.tail = FALSE)
sx <- head(sort(x), n = -3) # omits the largest three unique values
u <- seq(sx[4], tail(sx, n = 1), length.out = 129) # omit the first 3, too
y <- sapply(u, mef)
lines(u, y, col = "royalblue3")
## => Okay
## Note: One would typically only consider large (positive) losses x (for which
##       Pickands--Balkema--de Haan is valid) for EVT modeling purposes.


### 2 F in the Gumbel domain of attraction #####################################

### 2.1 Exp(lam) ###############################################################

## In the Gumbel domain, we have xi = 0, so the limiting GPD(0, beta) has
## e(u) = beta / (1 - xi). Since Exp(lam) = GPD(1/lam), we should have e(u) = 1/lam.
lam <- 3
stopifnot(lam > 0)
set.seed(271)
x <- rexp(n, rate = lam)
mean_excess_plot(x)
abline(a = 1/lam, b = 0, col = "royalblue3")
## => Okay

## Note in particular that for xi = 0, the slope is 0.


### 2.2 N(0, 1) ################################################################

## N(0, 1) is not a special case of a GPD but in the MDA of a Gumbel.
## Hence xi is 0 and e(u) = (beta(u) + 0 * u) / (1 - 0) = beta(u) eventually.
set.seed(271)
x <- rnorm(n)
mean_excess_plot(x)
## Formula (*) and the substitution v = Phi(x) provide us with e(u)
## (with numerical integration on a compact interval)
f <- function(v) (1-v) / dnorm(qnorm(v))
mef <- function(u) integrate(f, lower = pnorm(u), upper = 1)$value / pnorm(u, lower.tail = FALSE)
sx <- head(sort(x), n = -3) # omits the largest three unique values
u <- seq(sx[4], tail(sx, n = 1), length.out = 129) # omit the first 3, too
y <- sapply(u, mef)
lines(u, y, col = "royalblue3")
## Note: One would typically only consider large (positive) losses x (for which
##       Pickands--Balkema--de Haan is valid) for EVT modeling purposes.


### 2.3 LN(0, 1) ###############################################################

## Same here (but more heavy-tailed than before).
set.seed(271)
x. <- rlnorm(n)
stopifnot(all.equal(log(x.), x)) # sanity check
mean_excess_plot(x.)
## Formula (*) and substitutions y = log(x) and v = Phi(y) provide us with e(u)
## (with numerical integration on a compact interval)
f <- function(v) {
    q.v <- qnorm(v)
    (1-v) * exp(q.v) / dnorm(q.v)
}
mef <- function(u) integrate(f, lower = pnorm(log(u)), upper = 1)$value / pnorm(log(u), lower.tail = FALSE)
sx <- head(sort(x.), n = -3) # omits the largest three unique values
u <- seq(0, tail(sx, n = 1), length.out = 129)
y <- sapply(u, mef)
lines(u, y, col = "royalblue3")
## => A bit more difficult here as more heavy tailed
## => Virtually indistinguishable from xi > 0 case!


### 3 F in the Weibull domain of attraction ####################################

## U(0,1) is in the MDA of a Weibull, so xi < 0. We obtain e(u) from (*)
## as e(u) = (1/2 - u * (1-u/2)) / (1-u)
set.seed(271)
x <- runif(n)
mean_excess_plot(x, col = adjustcolor("black", alpha.f = 0.002))
sx <- head(sort(x), n = -3) # omits the largest three unique values
u <- seq(0, tail(sx, n = 1), length.out = 129)
mef <- function(u) (0.5-u*(1-u/2)) / (1-u)
y <- sapply(u, mef)
lines(u, y, col = "royalblue3", lwd = 2)

## Note in particular that for xi < 0, the slope is negative.
