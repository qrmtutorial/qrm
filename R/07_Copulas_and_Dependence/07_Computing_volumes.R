## By Marius Hofert

## Computing C-volumes


### Setup ######################################################################

library(copula)

## Volume we consider
a <- c(1/4, 1/2) # lower left end point
b <- c(1/3, 1) # upper right end point
stopifnot(0 <= a, a <= 1, 0 <= b, b <= 1, a <= b)


### 1 Independence copula ######################################################

## Define the independence copula
ic <- indepCopula()

## Compute the Pi_{(a,b]} volume
vol <- prob(ic, l = a, u = b)

## Let's manually compute this volume. Under independence, the probability of falling
## in the two-dimensional interval (a, b] is (b[1] - a[1]) * (b[2] - a[2])
p <- (b[1] - a[1]) * (b[2] - a[2])

## Check whether they are (numerically) the same
stopifnot(all.equal(vol, p))

## d = 3
vol. <- prob(indepCopula(dim = 3), l = c(a, 1/8), u = c(b, 3/8))
p. <- p * (3/8 - 1/8)
stopifnot(all.equal(vol., p.))


### 2 Frechet--Hoeffding lower bound W #########################################

## Define W
W <- function(u1, u2) max(u1 + u2 - 1, 0)

## Compute the W_{(a,b]} volume
vol <- W(b[1], b[2]) - W(b[1], a[2]) - W(a[1], b[2]) + W(a[1], a[2])

## Let's manually compute this volume. Under W, the probability of falling
## in the two-dimensional interval (a, b] is given by the length of the line
## segment cut out of the secondary diagonal by (a, b]. Given how a and b
## are located here, this is:
p <- sqrt(2*(b[1] - a[1])^2) / sqrt(2)

## Check whether they are (numerically) the same
stopifnot(all.equal(vol, p))

## Compute the volume of a rectangle below or above the secondary diagonal
## Rectangle below the secondary diagonal
a. <- a
b. <- c(1/3, 1/2)
stopifnot(0 <= a., a. <= 1, 0 <= b., b. <= 1, a.[1]+a.[2] <= 1, b. >= a.)
W(b.[1], b.[2]) - W(b.[1], a.[2]) - W(a.[1], b.[2]) + W(a.[1], a.[2])
## Rectangle above the secondary diagonal
a. <- c(2/3, 1/2)
b. <- c(3/4, 3/4)
stopifnot(0 <= a., a. <= 1, 0 <= b., b. <= 1, a.[1]+a.[2] >= 1, b. >= a.)
W(b.[1], b.[2]) - W(b.[1], a.[2]) - W(a.[1], b.[2]) + W(a.[1], a.[2])


### 3 t copula #################################################################

## Define the t copula
nu <- 3 # note: pCopula() is currently only available for integer degrees of freedom
tc <- tCopula(iTau(tCopula(df = nu), tau = 0.5), df = nu)

## Compute the C_{(a,b]} volume
vol <- prob(tc, l = a, u = b)

## Let's manually compute this volume and check by Monte-Carlo
p <- pCopula(b, copula = tc) - pCopula(c(b[1], a[2]), copula = tc) -
    pCopula(c(a[1], b[2]), copula = tc) + pCopula(a, copula = tc)

## Check whether they are (numerically) the same
stopifnot(all.equal(vol, p))

## Let's also compute it by MC by approximating the volume by the relative
## frequency of samples falling in the two-dimensional interval (a, b]
N <- 1e6 # (large) sample size
set.seed(271) # set a seed (for reproducibility)
U <- rCopula(N, copula = tc) # sample the t copula
volMC <- mean(a[1] < U[,1] & U[,1] <= b[1] & a[2] < U[,2] & U[,2] <= b[2])

## Check whether we get (numerically) the same
stopifnot(all.equal(volMC, p, tol = 1e-2))
