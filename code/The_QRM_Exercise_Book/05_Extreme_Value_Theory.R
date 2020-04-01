## By Marius Hofert

## R script for Chapter 5 of The QRM Exercise Book


### Setup ######################################################################

library(xts)
library(qrmdata)
library(qrmtools)
doPDF <- require(crop)
options(digits = 10)


### Exercise 5.7 (Threshold choice in the peaks-over-threshold method) #########

## Reproducing code

## Simulated data
n <- 1e4 # need large sample size here
shape <- 1/2
scale <- 1/10
stopifnot(0 < shape, shape < 1, scale > 0)
set.seed(271)
L.sim <- rGPD(n, shape = shape, scale = scale)
shape/(1-shape) # slope of the mean excess plot

## Mean excess plot of simulated data
if(doPDF) pdf(file = (file <- "fig_05_POT_threshold_choice_sim.pdf"))
opar <- par(pty = "s")
mean_excess_plot(L.sim[L.sim > 0])
## abline(a = scale/(1-shape), b = shape/(1-shape)) # true mean excess function
par(opar)
if(doPDF) dev.off.crop(file)

## Ethereum data
data(crypto)
ETH <- crypto['2015-08-06/2017-12-31', "ETH"]
X <- returns(ETH) # log-returns
L.ETH <- -X # losses; here: negative log-returns
length(L.ETH)
length(L.ETH[L.ETH > 0])
u <- 0.05
fit <- fit_GPD_MLE(L.ETH[L.ETH > u]-u) # fit a GPD
fit$par["shape"] # shape ~= 0.25 = 1/4
fit$par["shape"]/(1-fit$par["shape"]) # slope of the mean excess plot xi/(1-xi) ~ 1/3

## Mean excess plot of ETH data
if(doPDF) pdf(file = (file <- "fig_05_POT_threshold_choice_ETH.pdf"))
opar <- par(pty = "s")
mean_excess_plot(L.ETH[L.ETH > 0])
par(opar)
if(doPDF) dev.off.crop(file)


### Exercise 5.21 (Block maxima method applied to S&P 500 return data) #########

## a)

## Data preparation
data(SP500)
S <- SP500
X <- returns(S) # log-returns
L <- -X # losses; here: negative log-returns
plot.zoo(L, xlab = "Time", ylab = "S&P 500 losses (-log-returns)")

## Find losses related to the financial crisis of 2007--2008
u <- 0.075
abline(h = u, col = "maroon3")
ii <- which((L > 0.075) & (index(L) >= "2000-01-01")) # after Black Monday event
index(L)[ii] # dates of large losses due to financial crisis
L[ii] # corresponding losses
## Note: The loss on 2008-09-29 was due to Congress rejecting a 700B USD plan to
##       rescue the banking industry; see http://money.cnn.com/2013/06/20/investing/stocks-markets/index.html
(end <- index(L)[min(ii)-1]) # take this (Fri, 2008-09-26) as 'last available'
crash <- index(L)[min(ii)] # date of crash (Mon, 2008-09-29)
L[crash] # loss

## Drops in S&P 500
## Recall: L_t = -log(S_t/S_{t-1}) = -log(1 + beta), where beta denotes the 'gain'
##         => drop = -beta = -(exp(-L_t)-1) = -expm1(-L_t)
-expm1(-sum(L['2008-09-22/2008-09-26'])) # ~= 3.33% (drop since Mon)
-expm1(-sum(L['2008-08-26/2008-09-26'])) # ~= 4.23% (drop past month)
-expm1(-sum(L['2008-05-26/2008-09-26'])) # ~= 11.82% (drop past quarter)

## Data we work with
time <- c("1960-01-01", as.character(end)) # from 1960 until 'end'
L. <- -returns(S[paste0(time, collapse = "/")]) # grab out the data we work with

## Plot the S&P 500 -log-returns
plot.zoo(L., main = "S&P 500 losses (-log-returns)",
         xlab = "Time t", ylab = expression(L[t] == -log(S[t]/S[t-1])))


## Block Maxima Method (BMM)
M <- period.apply(L., INDEX = endpoints(L., "years"), FUN = max) # yearly maxima
fit <- fit_GEV_MLE(M) # GEV maximum likelihood estimator
stopifnot(fit$convergence == 0) # => converged
(xi  <- fit$par[["shape"]]) # ~= 0.5278 => Frechet domain
(mu  <- fit$par[["loc"]])
(sig <- fit$par[["scale"]])
fit$SE # standard errors
ceiling(1/xi) # => infinite 2nd moment


## b)

## Q: What is the probability that next year's maximal risk-factor change
##    exceeds the maximal loss over the past 20 years?
##    Note: 1-pGEV(max(head(M, n = -1)), shape = xi, loc = mu, scale = sig) # ~= 0.83% (because of Black Monday!)
1-pGEV(max(tail(head(M, n = -1), n = 20)), shape = xi, loc = mu, scale = sig) # exceedance prob. ~= 7.95%


## c)

## Q: What is the 10-year and 50-year return level? ... so the loss we expect
##    to be exceeded once every 10 (or 50) years.
##    Recall: k n-block return level = r_{n,k} = H^-(1-1/k) = level which is
##            expected to be exceeded in one out of every k n-blocks.
qGEV(1-1/10, shape = xi, loc = mu, scale = sig) # r_{n = 260, k = 10} ~=  6.32%; n ~ 1y
qGEV(1-1/50, shape = xi, loc = mu, scale = sig) # r_{n = 260, k = 50} ~= 14.50%


## d)

## Q: What is the return period of a risk-factor change at least as large as
##    on 'crash'? ... so the number of n-blocks for which we expect to see
##    at least one of them exceeding a loss as large as on 'crash'.
##    Recall: k_{n,u} = 1/\bar{H}(u) = period (= number of n-blocks) in which we
##            expect to see a single n-block exceeding u (= risk-factor change
##            as on 'crash')
1/(1-pGEV(as.numeric(L[crash]),
          shape = xi, loc = mu, scale = sig)) # ~= 20.78 years
## => 2008 - 21 = 1987 (Black Monday!)


### Exercise 5.22 (Peaks-over-threshold analysis of Bitcoin data) ##############

## a)

## Data preparation
## Note: The analysis differs depending on whether the huge loss shortly after
##       2014 is included.
data(crypto)
time <- c("2014-01-01", "2017-12-31")
time. <- paste0(time, collapse = "/")
BTC <- crypto[time., "BTC"] # price of 1 BTC in USD
X <- returns(BTC) # log-returns
L <- -X # losses; here: negative log-returns

## Plots
plot.zoo(BTC, xlab = "Time", ylab = "Price of 1 BTC in USD") # BTC
plot.zoo(L, main = "", xlab = "Time", ylab = "BTC -log-returns") # -log-returns


## b)

## Sample mean excess plot
mean_excess_plot(L[L > 0], xaxt = "n", yaxt = "n")
axis(1, at = seq(0, 0.25, by = 0.02))
axis(2, at = seq(0, 0.35, by = 0.02))
u <- 0.064 # another option: 0.075 (but less excesses then)
abline(v = u, lty = 3, lwd = 1.6)
legend("bottomright", bty = "n", lty = c(NA, 3), lwd = c(1, 1.6), pch = c(1, NA),
       legend = c("Mean excess", "Threshold choice"))

## Effect of changing the threshold on xi
GPD_shape_plot(L) #, thresholds = seq(quantile(L, 0.5), quantile(L, 0.99), length.out = 65))
abline(v = u, lty = 4, lwd = 1.6, col = "royalblue3") # threshold
abline(h = 0.5, lty = 3, lwd = 1.6, col = "darkorange2") # models above line have infinite variance (E(X^k) = Inf <=> xi >= 1/k; k = 2 here)
abline(h =   1, lty = 3, lwd = 1.6, col = "maroon3") # models above line have infinite mean
legend("topleft", bty = "n", lty = c(1, 2, 4, 3, 3), lwd = c(1, 1, 1.4, 1.6, 1.6),
       col = c(rep("black", 2), "royalblue3", "maroon3", "darkorange2"),
       legend = c("GPD shape", "95% CIs", "Threshold choice", "Infinite mean", "Infinite variance"))
## We see that our choice of u leads to a (roughly) stabilizing fitted GPD shape
## parameter (before confidence intervals get wider). Our choice of u leads to
## an infinite variance (but finite mean) model. Especially infinite variance
## models lie well between the pointwise asymptotic confidence intervals for
## thresholds such as u or larger.


## c)

## Compute GPD MLE
exceed <- L[L > u] # exceedances
excess <- exceed - u # excesses
(fit <- fit_GPD_MLE(excess)) # MLE
(shape <- fit$par[["shape"]]) # => infinite variance model; see ceiling(1/shape)
(scale <- fit$par[["scale"]])

## Q-Q plot
## (alternatively: qq_plot(exceed, FUN = function(p) qGPD(p, shape = shape, scale = scale) + u))
qq_plot(excess, FUN = function(p) qGPD(p, shape = shape, scale = scale))

## Plot empirical exceedance loss distribution function overlaid with the
## shifted fitted GPD
## Note: Replacing 'x-u' by 'x' would give the empirical excess distribution
##       function F[u](x) (= GPD(x))
res <- edf_plot(exceed, do.points = FALSE, ylab = "Exceedance loss distribution function")
x <- tail(res$t, n = -1)
lines(x, pGPD(x-u, shape = shape, scale = scale), col = "royalblue3") # shifted fitted GPD
legend("bottomright", bty = "n", lty = c(1, 1), col = c("black", "royalblue3"),
       legend = c("empirical", expression(F[u](x-u)~"for"~x>=u~"and"~F[u]~"being the fitted GPD")))


## d)

## Semi-parametric VaR_alpha and ES_alpha estimates
(p.exceed <- mean(L > u))
alpha <- c(0.975, 0.99)
(VaR.u <- VaR_GPDtail(alpha[2], threshold = u, p.exceed = p.exceed,
                      shape = shape, scale = scale))
(ES.u <- ES_GPDtail(alpha[1], threshold = u, p.exceed = p.exceed,
                    shape = shape, scale = scale))


## e)

## Empirical tail probabilities with Smith estimator overlaid
res <- tail_plot(L, threshold = u, shape = shape, scale = scale,
                 ylab = "Tail probability 1-F(x)")
axis(4, at = 1 - alpha, # 4th axis
     labels = as.expression(lapply(1:2, function(k) substitute(1-a, list(a = alpha[k])))))
abline(v = VaR.u, lty = 2, lwd = 1.6, col = "royalblue3") # VaR_alpha
abline(v = ES.u,  lty = 2, lwd = 1.6, col = "maroon3") # ES_alpha
legend("topright", bty = "n", inset = 0.02,
       pch = c(1, rep(NA, 3)), lty = c(NA, 1, 2, 2), lwd = c(NA, 1, 1.6, 1.6),
       col = c(rep("black", 2), "royalblue3", "maroon3"),
       legend = as.expression(c("Empirical tail probability", "Tail estimator",
                                substitute(VaR[a], list(a = alpha[2])),
                                substitute(ES[a], list(a = alpha[1])))))

## Sanity check
exceed <- res$np[,"exceed"]
Fn.bar.exceed <- res$np[,"Fn.bar.exceed"]
f <- function(x) (sum(L > x) + 1/2) / length(L) # or take 'exceed' instead of 'L' (equivalent here since x >= u)
Fn.bar.manual <- sapply(exceed, f)
stopifnot(all.equal(Fn.bar.manual, Fn.bar.exceed))
