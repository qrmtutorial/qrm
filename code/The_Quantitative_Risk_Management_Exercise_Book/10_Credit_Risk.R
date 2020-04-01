## By Marius Hofert

## R script for Chapter 10 of The QRM Exercise Book


### Setup ######################################################################

options(digits = 10)


### Exercise 10.3 (Rating transition matrices) #################################

## Second part of the question
M <- matrix(c(0.8, 0.15, 0.05,
              0.1, 0.7,  0.2,
              0,   0,    1), ncol = 3, byrow = TRUE)
M %*% M


### Exercise 10.4 (Credit risk in the Merton model) ############################

## b)
pnorm((log(1/2) - (0.1-0.04/2)) / 0.2)


### Exercise 10.8 (Risk-neutral valuation) #####################################

## c)
A <- matrix(c(0.4, 1.025,
              1,   1.025), ncol = 2, byrow = TRUE)
b <- c(1, 0)
solve(A, b = b)


### Exercise 10.14 (Pricing and calibration of a CDS) ##########################

## b)
obj <- function(gamma, x.star = 0.0042, r = 0.02, delta = 0.6) # note: use -expm1() instead of 1-exp()
    (x.star/4) * exp(-(r + gamma)/4) + delta * (gamma / (r + gamma)) * expm1(-(r + gamma)/4)
uniroot(obj, interval = c(0, 1))


### Exercise 10.17 (Simplified CDS pricing) ####################################

## Reproducing code

## c)

## Q(tau_R > t) = E(exp(-\int_t^T \Psi_s ds)) = exp(alpha(T-t) + beta(T-t) * \Psi_t);
## see MFE (2015, p. 419)
Q <- function(t, Psi.t, kappa, th, sig, rho.0 = 0, rho.1 = 1) {
    tau <- t - 0
    gamma <- sqrt(kappa^2 + 2 * sig^2 * rho.1)
    denom <- gamma - kappa + exp(gamma * tau) * (gamma + kappa)
    alpha <- -rho.0 * tau + 2 * (kappa * th / sig^2) * log((2 * gamma * exp(tau * (gamma + kappa) / 2)) / denom)
    beta <- (-2 * rho.1 * expm1(gamma * tau)) / denom
    exp(alpha + beta * Psi.t)
}

## p_0(0, t)
p0 <- function(t, r) exp(-r * t)

## Premium leg for x = 1 (see MFE (2015, (10.39) for a deterministic hazard rate, t = 0)
V0_prem_leg <- function(t, r, Psi.t, kappa, th, sig, ...)
    sum(p0(t, r) * diff(c(0, t)) * Q(t, Psi.t = Psi.t, kappa = kappa, th = th, sig = sig, ...))

## Simplified default leg
V0_def_leg <- function(t, r, Psi.t, kappa, th, sig, ...) {
    p0.vec <- p0(t, r) # p_0(0, t_i)'s
    delta * (p0(t[1], r) - p0(t[N], r) * Q(t[N], Psi.t = Psi.t, kappa = kappa, th = th, sig = sig) +
             sum(Q(t[1:(N-1)], Psi.t = Psi.0, kappa = kappa, th = th, sig = sig) * diff(p0(t, r))))
}

## Setup
N <- 20 # number of discretization points
T <- 5 # maturity
t <- ((1:N)/N) * T # discretization points t_i's (\Delta_t = 1/4)
Psi.0 <- 0.03 # CIR specification at 0
kappa <- 1 # CIR parameter
th <- 0.04 # CIR parameter
sig <- 0.2 # CIR parameter
r <- 0.01 # interest rate
stopifnot(kappa * th >= sig^2 / 2) # see MFE (2015, p. 419)
delta <- 0.6 # LGD

## Corresponding fair spread
V0_def_leg(t, r = r, Psi.t = Psi.0, kappa = kappa, th = th, sig = sig) /
    V0_prem_leg(t, r = r, Psi.t = Psi.0, kappa = kappa, th = th, sig = sig) # 0.09994508 / 4.436493 = 0.02252794

## For larger sig
sig. <- 0.25
V0_def_leg(t, r = r, Psi.t = Psi.0, kappa = kappa, th = th, sig = sig.) /
    V0_prem_leg(t, r = r, Psi.t = Psi.0, kappa = kappa, th = th, sig = sig.) # 0.02236618

## For smaller sig
sig.. <- 0.15
V0_def_leg(t, r = r, Psi.t = Psi.0, kappa = kappa, th = th, sig = sig..) /
    V0_prem_leg(t, r = r, Psi.t = Psi.0, kappa = kappa, th = th, sig = sig..) # 0.02265696
