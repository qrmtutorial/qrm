## By Marius Hofert

## R script for Chapter 11 of The QRM Exercise Book


### Setup ######################################################################

library(copula)
doPDF <- require(crop)
options(digits = 10)


### Exercise 11.6 (Calculating moments of portfolio credit loss distributions) #

m  <- 1000
delta <- 0.4
p <- 0.01

## a)
delta * m * p # E(L)
sqrt(delta^2 * m * p * (1-p)) # sd(L)

## b)
delta * m * p # E(L)
rho <- 0.005
sqrt(delta^2 * m * p * (1-p) * (1 + (m-1) * rho)) # sd(L)

## c)
delta * m * p # E(L)
rho <- 0.005
a <- 9.2 # a ...
b <- 13.8 # ... and b ...
mean.Delta.1 <- a/(a+b) # ... are chosen so that E(\Delta_1) = 0.4 = delta
stopifnot(all.equal(mean.Delta.1, delta))
var.Delta.1 <- a*b/((a+b)^2 * (a+b+1))
sqrt(m * p * var.Delta.1 + m * p * (1-p) * mean.Delta.1^2 * (1 + (m-1) * rho)) # sd(L)


### Exercise 11.7 (Correlation bounds for default indicators) ##################

## a) Reproducing code

##' @title Correlation bounds for B(1, p_.) margins (e.g., default indicators)
##' @param p 2-column matrix or 2-vector of probabilities
##' @param method character string indicating the minimal or maximal correlation
##'        bound to be plotted
##' @return correlation bound
##' @author Marius Hofert
cor_bound_B <- function(p, method = c("max", "min"))
{
    ## p = (p_1, p_2)
    if(!is.matrix(p)) p <- rbind(p, deparse.level = 0)
    method <- match.arg(method)
    C.bounds <- if(method == "min") pmax(p[,1]+p[,2]-1, 0) else pmin(p[,1], p[,2])
    (C.bounds - p[,1]*p[,2]) / sqrt(p[,1]*(1-p[,1]) * p[,2]*(1-p[,2]))
}

## B(1, p.) margins
n.grid <- 26 # number of grid points in each dimension
p <- seq(0.01, 0.99, length.out = n.grid) # subdivision points in each dimension
grid <- expand.grid("p[1]" = p, "p[2]" = p) # build a grid
val.min <- cbind(grid, "rho[min](p[1],p[2])" =
                 cor_bound_B(grid, method = "min"))
val.max <- cbind(grid, "rho[max](p[1],p[2])" =
                 cor_bound_B(grid))
if(doPDF) pdf(file = (file <- "fig_11_cor_bounds_B_low.pdf"))
wireframe2(val.min) # minimal correlation bound
if(doPDF) dev.off.crop(file)
if(doPDF) pdf(file = (file <- "fig_11_cor_bounds_B_up.pdf"))
wireframe2(val.max) # maximal correlation bound
if(doPDF) dev.off.crop(file)


### Exercise 11.12 (Exchangeable beta mixture model) ###########################

## Input parameters
m <- 1000 # number m of loans
pi. <- 0.01 # individual default probability pi = pi_1 (of each loan)
rho. <- 0.005 # default correlation rho_Y = Cor(Y_{j_1}, Y_{j_2}) (Y_. = default indicator)


## b)

## Now consider the conditional independence case with mixing variable Q ~ B(a,b).
## M then follows a so-called 'beta-binomial' distribution.
(a <- pi. * (1/rho. - 1))
(b <- a*(1/pi. - 1))
stopifnot(a > 0, b > 0)


## e)

## Define the beta-binomial pmf of M
dbetaBinom <- function(x, size, a, b) # x = k, size = m
    exp(lchoose(size, x) + lbeta(a + k, b + size - k) - lbeta(a, b)) # see (11.16)

## Independent case; here M ~ B(m, pi.)
k <- 0:60 # number k of defaulted loans we consider (we plot P(M = k) for those k)
p.tilde <- dbinom(k, size = m, prob = pi.) # pmf

## Compute the pmf P(M = k) for the given k
p <- dbetaBinom(k, size = m, a = a, b = b)

## Plot P(M = k) for the given k including the independence case
if(doPDF) pdf(file = (file <- "fig_11_exchangeable_beta_mixture_pmf_comparison.pdf"))
ylim <- range(p.tilde, p)
plot(k, p, type = "l", ylim = ylim, col = "maroon3",
     xlab = "Number of losses", ylab = "")
lines(k, p.tilde) # add pmf under independence
legend("topright", bty = "n", lty = 1, col = c("maroon3", "black"),
       legend = c("Probability under beta-binomial", "Probability under independence"))
if(doPDF) dev.off.crop(file)

## Plot their quotient
if(doPDF) pdf(file = (file <- "fig_11_exchangeable_beta_mixture_pmf_quotient.pdf"))
plot(k, p/p.tilde, type = "l", log = "y", xlab = "Number of losses", ylab = "")
abline(h = 1, lty = 2)
legend("topleft", bty = "n", lty = c(1, 2),
       legend = c(expression(frac("Probability under beta-binomial",
                                  "Probability under independence")),
                  "Reference line y = 1"))
if(doPDF) dev.off.crop(file)


### Exercise 11.16 (Large portfolio asymptotics application) ###################

## b)
lbar <- function(psi)
    (0.9/3) * ( pnorm((qnorm( 0.01) + sqrt(0.8) * psi)/sqrt(0.2)) +
                pnorm((qnorm(0.001) + sqrt(0.2) * psi)/sqrt(0.8)) )
(VaR.999 <- 30000 * lbar(qnorm(0.999))) # approximate VaR_0.999


### Exercise 11.17 (Large portfolio asymptotics with stochastic LGD) ###########

## b)
lbar <- function(psi)
    (1/6) * (1.2 * pnorm(qnorm(0.01) / 0.8 + (3/4) * psi) +
             2.4 * pnorm(qnorm(0.05) / 0.6 + (4/3) * psi))
(VaR.99 <- 30000 * lbar(qnorm(0.99))) # approximate VaR_0.99

## c)
(scale <- pnorm(0.5 + qnorm(0.99)) / 0.6)
(VaR.99. <- scale * VaR.99)