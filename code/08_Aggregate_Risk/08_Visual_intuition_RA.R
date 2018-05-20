## By Marius Hofert

## Visual intuition behind the Rearrangement Algorithm (RA)

library(sn) # skew normal and skew t distributions

## Parameters
a <- 0.95 # VaR confidence level
omega <- 10 # shape parameter of the skew t distribution
alpha <- 4 # slant parameter of the skew t distribution
nu <- 3 # degree of freedom

get_VaR <- function(a, omega, alpha, nu)
    qst(a, omega = omega, alpha = alpha, nu = nu)
get_ES <- function(a, omega, alpha, nu)
    integrate(qst, omega = omega, alpha = alpha, nu = nu, lower = a, upper = 1)$value / (1-a)

## Skewed t_3 density (assumed to be the density of L = L_1 + ... + L_d)
x <- seq(-8, 80, length.out = 201)
y <- dst(x, omega = omega, alpha = alpha, nu = nu)
plot(x, y, type = "l", xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "")
text(15, 0.8*max(y), labels = expression(f[L](x)))

## x axis
VaR <- get_VaR(a = a, omega = omega, alpha = alpha, nu = nu)
ES  <- get_ES (a = a, omega = omega, alpha = alpha, nu = nu)
axis(1, at = VaR, label = expression(VaR[alpha](L)), hadj = 0.25) # ES
axis(1, at = ES,  label = expression(ES[alpha](L)),  hadj = 0.25) # ES

## Shaded alpha region
sq <- seq(VaR, max(x), length.out = 201)
poly.x <- c(VaR, sq, max(x))
poly.y <- c(0, dst(sq, omega = omega, alpha = alpha, nu = nu), 0)
polygon(poly.x, poly.y, col = "gray80", border = NA)
text(0.44*diff(range(x)), 0.07*max(y), labels = expression(1-alpha), col = "gray50")

## Overlay peaked densities
d <- ES-VaR
x <- seq(ES-d, ES+d, length.out = 201)
z1 <- (1-a) * dnorm(x, mean = ES, sd = 4)
z2 <- (1-a) * dnorm(x, mean = ES, sd = 0.8)
lines(x, z1, lty = 2)
lines(x, z2, lty = 3)
