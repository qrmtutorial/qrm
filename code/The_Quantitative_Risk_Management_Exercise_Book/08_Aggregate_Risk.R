## By Marius Hofert

## R script for Chapter 8 of The QRM Exercise Book


### Setup ######################################################################

library(sn)
library(copula)
library(qrmtools)
doPDF <- require(crop)
options(digits = 10)


### Exercise 8.23 (A basic version of the rearrangement algorithm) #############

## Reproducing code

## b) i) Rearrange A until all columns are oppositely ordered to the sum of all others
(A <- matrix(c(2, 3, 8,
               4, 1, 1,
               1, 5, 2,
               3, 7, 4), ncol = 3, byrow = TRUE))
rearrange(A, tol = NULL, sample = FALSE, trace = TRUE) # rearrange


## b) ii) Plot of rearranged sample

## Define parameters of the three margins
th <- 2.5 # Pareto parameter
m <- 10 # mean of the log-normal
v <- 20 # variance of the log-normal
s <- 4 # shape of the gamma underlying the log-gamma
r <- 5 # rate of the gamma underlying the log-gamma
## Define list of marginal dfs
qF <- list(qPar = function(p) (1 - p)^(-1/th) - 1,
           qLN  = function(p) qlnorm(p, meanlog = log(m)-log(1+v/m^2)/2,
                                          sdlog = sqrt(log(1+v/m^2))),
           qLG  = function(p) exp(qgamma(p, shape = s, rate = r)))
## Apply ARA()
set.seed(271) # for reproducibility
alpha <- 0.99 # confidence level
wVaR <- ARA(alpha, qF = qF) # compute worst VaR (bounds)
X <- wVaR[["X.rearranged"]]$up # extract rearranged matrix (upper bound)
U <- pobs(X) # compute pseudo-observations
colnames(U) <- c("U[1]", "U[2]", "U[3]")
## Plot
if(doPDF) pdf(file = (file <- "fig_08_RA_worst_VaR_copula_sample.pdf"))
cloud2(U, screen = list(z = 15, x = -60)) # default 'screen': z = 40, x = -60
if(doPDF) dev.off.crop(file)


## b) iii) Visual intuition (see also R code Chapter 8)

## Parameters
a <- 0.95 # VaR confidence level
omega <- 10 # shape parameter of the skew t distribution
alpha <- 4 # slant parameter of the skew t distribution
nu <- 3 # degree of freedom

## Auxiliary functions for exactly locating VaR and ES
get_VaR <- function(a, omega, alpha, nu)
    qst(a, omega = omega, alpha = alpha, nu = nu)
get_ES <- function(a, omega, alpha, nu)
    integrate(qst, omega = omega, alpha = alpha, nu = nu, lower = a, upper = 1)$value / (1-a)

## Skewed t_3 density (assumed to be the density of L = L_1 + ... + L_d)
x <- seq(-8, 80, length.out = 201) # x values
y <- dst(x, omega = omega, alpha = alpha, nu = nu) # skewed t density
VaR <- get_VaR(a = a, omega = omega, alpha = alpha, nu = nu) # VaR
ES  <- get_ES (a = a, omega = omega, alpha = alpha, nu = nu) # ES
d <- ES-VaR # difference
x. <- seq(ES-d, ES+d, length.out = 201) # x values for evaluating the 'tail densities'
z1 <- (1-a) * dnorm(x., mean = ES, sd = 4) # wider 'tail density'
z2 <- (1-a) * dnorm(x., mean = ES, sd = 0.8) # narrower 'tail density'

## Plot
if(doPDF) pdf(file = (file <- "fig_08_RA_intuition.pdf"), width = 8, height = 6)
plot(x, y, type = "l", xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "")
text(15, 0.8*max(y), labels = expression(f[L](x)))
## x axis
axis(1, at = VaR, label = expression(VaR[alpha](L)), hadj = 0.25) # ES
axis(1, at = ES,  label = expression(ES[alpha](L)),  hadj = 0.25) # ES
## Shaded alpha region
sq <- seq(VaR, max(x), length.out = 201)
poly.x <- c(VaR, sq, max(x))
poly.y <- c(0, dst(sq, omega = omega, alpha = alpha, nu = nu), 0)
polygon(poly.x, poly.y, col = "gray80", border = NA)
text(0.44*diff(range(x)), 0.07*max(y), labels = expression(1-alpha), col = "gray50")
## Overlay peaked densities
lines(x., z1, lty = 2)
lines(x., z2, lty = 3)
if(doPDF) dev.off.crop(file)
