## By Marius Hofert

## R script for Chapter 2 of The QRM Exercise Book


### Setup ######################################################################

library(sfsmisc) # for eaxis()
library(qrmtools)
doPDF <- require(crop)
options(digits = 10)


### Exercise 2.13 (VaR and ES for a distribution function with jumps) ##########

## a) Reproducing code

## Compute values of F
x1 <- seq(1, 3, length.out = 401) # x values after first jump
y1 <- 1 - 1 / (1 + x1) # values of F at x1
x2 <- seq(3, 5, length.out = 401) # x values after second jump
y2 <- 1 - 1 / x2^2 # values of F at x2

## Plot
if(doPDF) pdf(file = (file <- "fig_02_df_with_jumps.pdf"))
opar <- par(pty = "s")
plot(NA, xlim = range(0, x1, x2), ylim = 0:1, xlab = "x", ylab = "F(x)")
lines(c(0, 1), c(0, 0))
lines(x1, y1)
points(min(x1), min(y1))
lines(x2, y2)
points(min(x2), min(y2))
par(opar)
if(doPDF) dev.off.crop(file)


### Exercise 2.17 (VaR and ES for bivariate normal risks) ######################

## Covariance matrices of X
sig.d.1 <- 0.2 /sqrt(250)
sig.d.2 <- 0.25/sqrt(250)
rho  <- 0.4
rho. <- 0.6
(Sigma  <- matrix(c(sig.d.1^2, rho *sig.d.1*sig.d.2, rho *sig.d.1*sig.d.2, sig.d.2^2), ncol = 2))
(Sigma. <- matrix(c(sig.d.1^2, rho.*sig.d.1*sig.d.2, rho.*sig.d.1*sig.d.2, sig.d.2^2), ncol = 2))

## Setup
w <- c(0.7, 0.3)
V.t <- 10^6
alpha <- 0.99
q <- qnorm(alpha)

## VaR and ES for rho = 0.4
(VaR.99 <- V.t * sqrt(drop(t(w) %*% Sigma  %*% w)) * q)
stopifnot(all.equal(VaR.99, VaR_t(0.99, loc = 0, scale = V.t * sqrt(drop(t(w) %*% Sigma %*% w)), df = Inf))) # sanity check
(ES.99  <- V.t * sqrt(drop(t(w) %*% Sigma  %*% w)) * dnorm(q)/(1-alpha))
stopifnot(all.equal(ES.99, ES_t(0.99, loc = 0, scale = V.t * sqrt(drop(t(w) %*% Sigma %*% w)), df = Inf))) # sanity check

## VaR and ES for rho = 0.6
(VaR.99 <- V.t * sqrt(drop(t(w) %*% Sigma. %*% w)) * q)
stopifnot(all.equal(VaR.99, VaR_t(0.99, loc = 0, scale = V.t * sqrt(drop(t(w) %*% Sigma. %*% w)), df = Inf))) # sanity check
(ES.99  <- V.t * sqrt(drop(t(w) %*% Sigma. %*% w)) * dnorm(q)/(1-alpha))
stopifnot(all.equal(ES.99, ES_t(0.99, loc = 0, scale = V.t * sqrt(drop(t(w) %*% Sigma. %*% w)), df = Inf))) # sanity check


### Exercise 2.18 (Basic historical simulation) ################################

## Monthly log-returns
X.tp1.1 <- c(-16.1, 5.1, -0.4, -2.5, -4, 10.5, 5.2, -2.9, 19.1,  0.4) / 100 # X_{t+1,1}
X.tp1.2 <- c( -8.2, 3.1,  0.4, -1.5, -3,  4.5, 2,   -3.7, 10.9, -0.4) / 100 # X_{t+1,2}

## a)
(L <- -2 * 1e3 * X.tp1.1) # -lambda_1 S_{t,1} X_{t+1,1}
(VaR.90.a <- quantile(L, probs = 0.9, names = FALSE, type = 1)) # VaR_{0.9}
## Note: Changing sort(-X.tp1.1)[9] changes VaR.90.a


## b)
(L <- -1 * 1e3 * X.tp1.1 - 10 * 100 * X.tp1.2) # -lambda_1 S_{t,1} X_{t+1,1} - lambda_2 S_{t,2} X_{t+1,2}
(VaR.90.b <- quantile(L, probs = 0.9, names = FALSE, type = 1)) # VaR_{0.9}
## Note: Changing the second-largest loss of S_1 changes VaR.90.b:
##       1) X.tp1.2[5] = -5 < -4 = X.tp1.1[5] => VaR.90.b = 90 > 80 = VaR.90.a
##       2) X.tp1.2[5] = -4      = X.tp1.1[5] => VaR.90.b = 80      = VaR.90.a
##       3) X.tp1.2[5] = -3 > -4 = X.tp1.1[5] => VaR.90.b = 70 < 80 = VaR.90.a


### Exercise 2.22 (Subadditivity of VaR for bivariate normal random variables) #

## c)

## Setup
n <- 1024 # number of alphas
alpha <- seq(1e-2, 1-1e-8, length.out = n) # range of confidence levels

## Compute VaR_alpha(L_1 + L_2) and VaR_alpha(L_1) + VaR_alpha(L_2)
## Specification
one <- as.matrix(rep(1, 2)) # column vector of 1s
mu <- rep(0, 2) # means
sd <- rep(1, 2) # standard deviations
rho <- c(-1, -0.5, 0, 0.5, 1) # correlations; note: for rho = 1, we obtain additivity (comonotonicity)
Sig <- lapply(rho, function(rho.) {
    res <- outer(sd, sd, function(s1, s2) rho. * s1 * s2)
    diag(res) <- sd^2
    res
})
mean <- as.numeric(t(one) %*% as.matrix(mu)) # mean of the sum
sig <- sapply(seq_along(rho), function(k) as.numeric((sqrt(t(one) %*% Sig[[k]] %*% one)))) # sd of the sum
## VaR_alpha(L_1) + VaR_alpha(L_2)
VaR <- qnorm(alpha, mean = mu[1], sd = sd[1]) + qnorm(alpha, mean = mu[2], sd = sd[2]) # independent of rho
## VaR_alpha(L_1 + L_2)
VaR.sum <- sapply(seq_along(rho), function(k) qnorm(alpha, mean = mean, sd = sig[k]))

## Plot
if(doPDF) pdf(file = (file <- "fig_02_superadd_joint_norm.pdf"))
pal <- colorRampPalette(c("black", "royalblue3", "darkorange2", "maroon3"), space = "Lab")
cols <- pal(length(rho)) # colors
opar <- par(pty = "s")
plot(alpha, VaR, ylim = range(VaR, VaR.sum), type = "l", lty = 2, lwd = 2,
     xlab = expression(alpha), ylab = expression(VaR[alpha]))
for(k in seq_along(rho))
    lines(alpha, VaR.sum[,k], col = cols[k])
abline(v = 1/2, lty = 2)
legend("topleft", bty = "n", lty = c(rep(1, length(rho)), 2), lwd = c(rep(1, length(rho)), 2),
       col = c(cols, "black"),
       legend = c(as.expression(lapply(seq_along(rho), function(k)
           substitute(VaR[alpha](L[1]+L[2])*", "*rho == r, list(r = rho[k])))),
           expression(VaR[alpha](L[1])+VaR[alpha](L[2]))))
par(opar)
opar <- par(usr = c(0,1,0,1))
text(0.3, 0.1, labels = "Superadditivity\nregion")
text(0.7, 0.1, labels = "Subadditivity\nregion")
par(opar)
mtext(expression("Jointly normal"~group("(",list(L[1], L[2]),")")~"with N(0,1) margins and correlation"~rho),
      side = 4, line = 1, adj = 0)
if(doPDF) dev.off.crop(file)


### Exercise 2.24 (Matching VaR and ES) ########################################

## a) Exact confidence level alpha for which ES_alpha = VaR_{0.99} for N(0,1)
h <- function(alpha) ES_t(alpha, df = Inf) - VaR_t(0.99, df = Inf)
eps <- 1e-8 # optimize over alpha in (eps, 1-eps)
uniroot(h, interval = c(eps, 1 - eps))$root


## b) Exact confidence level alpha for which ES_alpha = VaR_{0.99} for t_3.5
h <- function(alpha) ES_t(alpha, df = 3.5) - VaR_t(0.99, df = 3.5)
uniroot(h, interval = c(eps, 1 - eps))$root


### Exercise 2.26 (Superadditivity of VaR for two iid exponential random variables)

## c) Determine unique root in (0,1) numerically
h <- function(a) (1-a)*(1-2*log(1-a))-1
eps <- 1e-8 # optimize over alpha in (eps, 1-eps)
uniroot(h, interval = c(eps, 1 - eps))$root
