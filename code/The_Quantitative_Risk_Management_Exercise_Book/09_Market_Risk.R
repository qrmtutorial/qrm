## By Marius Hofert

## R script for Chapter 9 of The QRM Exercise Book


### Setup ######################################################################

library(RColorBrewer)
library(rugarch)
library(sfsmisc) # for eaxis()
library(nvmix) # for rNorm(), fitStudent() and rStudent()
library(zoo) # for na.fill()
library(xts) # for xts()
library(ADGofTest) # for ad.test()
library(qrmdata)
library(qrmtools)
doPDF <- require(crop)
options(digits = 10)


### Exercise 9.11 (Variance-covariance method for multivariate t risk-factor changes)

## c)

## Covariance matrix of X
sig.d.1 <- 0.2 /sqrt(250)
sig.d.2 <- 0.25/sqrt(250)
rho  <- 0.4
(SigmaX  <- matrix(c(sig.d.1^2, rho *sig.d.1*sig.d.2, rho *sig.d.1*sig.d.2, sig.d.2^2), ncol = 2))

## Setup
w <- c(0.7, 0.3)
V.t <- 10^6
alpha <- 0.99

## VaR and ES for nu = 5
nu <- 5
Sigma <- ((nu-2)/nu) * SigmaX
q <- qt(alpha, df = nu)
(VaR.99 <- V.t * sqrt(drop(t(w) %*% Sigma %*% w)) * q)
stopifnot(all.equal(VaR.99,
                    VaR_t(0.99, loc = 0, scale = V.t * sqrt(drop(t(w) %*% Sigma %*% w)), df = nu))) # sanity check
(ES.99  <- (1/(1-alpha)) * V.t * sqrt(drop(t(w) %*% Sigma %*% w)) * (nu/(nu-1)) * dt(q, df = nu) * (1+(q^2)/nu))
stopifnot(all.equal(ES.99,
                    ES_t(0.99, loc = 0, scale = V.t * sqrt(drop(t(w) %*% Sigma %*% w)), df = nu))) # sanity check


## d)

## VaR and ES for nu = Inf
nu <- Inf
Sigma <- (1 - 2/nu) * SigmaX # rewritten in order to correctly deal with nu = Inf
q <- qt(alpha, df = nu)
(VaR.99 <- V.t * sqrt(drop(t(w) %*% Sigma %*% w)) * q)
stopifnot(all.equal(VaR.99,
                    VaR_t(0.99, loc = 0, scale = V.t * sqrt(drop(t(w) %*% Sigma %*% w)), df = nu))) # sanity check
(ES.99  <- (1/(1-alpha)) * V.t * sqrt(drop(t(w) %*% Sigma %*% w)) * (1/(1-1/nu)) * dt(q, df = nu) * (1+(q^2)/nu)) # rewritten in order to correctly deal with nu = Inf
stopifnot(all.equal(ES.99,
                    ES_t(0.99, loc = 0, scale = V.t * sqrt(drop(t(w) %*% Sigma %*% w)), df = nu))) # sanity check


### Exercise 9.15 (Tests based on VaR violations) ##############################

## b)

binom.test(14, n = 500, p = 0.01)[["p.value"]] # 14 violations among 500 losses, 'violation' probability 0.01 (under H0)

## Note (under H0):
## p = E(I_{\{L_{t+1} > VaR_{alpha}^t\}}) = P(L_{t+1} > VaR_{alpha}^t)
##   = 1 - P(L_{t+1} <= VaR_{alpha}^t) = 1 - alpha = 1 - 0.99 = 0.01


## c)

## Data preparation
T <- c(7, 61, 70, 75, 90, 130, 200, 245, 367, 371, 385, 403, 406, 487) # violation days
S <- diff(c(0, T)) # spacings

## We are checking the iid property of the VaR_alpha violations (Bernoulli trials
## process). If the violations are iid, the spacings between VaR_alpha violations
## should be iid Geo(p) (on IN) with p = P(X_1 > VaR_alpha) = 1 - alpha = 0.01;
## see Exercise 3.15. Let us check this with a Kolmogorov--Smirnov test.
ks.test(S, y = function(q) pgeom(q, prob = 0.01)) # p-value ~= 0.01653
## => Strong indication against the iid property of the VaR_alpha violations.
## Note that by Exercise 3.15, we could also test against an Exp(p) distribution.

## Also the Anderson--Darling test rejects at 5% (but keep in mind the problem
## of multiple testing; still, here, both tests also lead to rejection at
## Bonferroni-corrected level 0.05/2).
ad.test(S, function(q) pgeom(q, prob = 0.01)) # p-value ~= 0.002711

## Even though the sample size is small, we could consider a Q-Q plot
qq_plot(S, FUN = function(p) qgeom(p, prob = 0.01)) # => far off


### Exercise 9.16 (Elicitability and comparison of VaR estimates) ##############

## a)

##' @title Scoring Function
##' @param y point forecast
##' @param l loss
##' @param alpha confidence level
##' @return scores (scoring function values)
##' @author Marius Hofert
S <- function(y, l, alpha)
{
    n <- length(l)
    if(length(y) == 1) y <- rep(y, n)
    stopifnot(length(y) == n)
    lley <- l <= y
    df <- abs(l - y)
    res <- numeric(n)
    res[lley] <- (1-alpha) * df[lley]
    res[!lley] <- alpha * df[!lley]
    res
}

## Evaluation and plot
y <- 1
l <- seq(-10, 10, length.out = 257)
alpha <- c(0.5, 0.9)
S. <- lapply(alpha, function(a) S(y, l = l, alpha = a))
if(doPDF) pdf(file = (file <- "fig_09_elicitability_score_function.pdf"))
opar <- par(mar = c(5, 4, 4, 2) + 0.1 + c(0, 1, 0, 0)) # # set margins (bottom, left, top, right)
plot(l, S.[[1]], type = "l", ylim = range(S.), xlab = "Loss level l",
     ylab = substitute("Score function"~{S[alpha]^q}(y,l)~"at l for y"==y., list(y. = y)))
lines(l, S.[[2]], col = "royalblue3")
legend("bottomright", bty = "n", lty = c(1, 1), col = c("black", "royalblue3"),
       legend = as.expression(lapply(1:2, function(k) substitute(alpha == a, list(a = alpha[k])))))
par(opar)
if(doPDF) dev.off.crop(file)


## b) Reproducing code

## Model specification (for simulation)
nu <- 3.5 # d.o.f. of the innovation distribution (of Z_t)
fixed.p <- list(mu = 0, # mu (intercept)
                omega = 1, # alpha_0 (intercept)
                alpha1 = 0.4, # alpha_1 (GARCH parameter of X_{t-1}^2)
                beta1 = 0.2, # beta_1 (GARCH parameter of sigma_{t-1}^2)
                shape = nu) # d.o.f. nu of the standardized t_nu innovation distribution
armaOrder  <- c(0, 0) # ARMA order
garchOrder <- c(1, 1) # GARCH order
varModel <- list(model = "sGARCH", garchOrder = garchOrder)
spec.true <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                        fixed.pars = fixed.p, distribution.model = "std") # t standardized residuals

## Simulate (X_t)
n <- 12 # sample size (= length of simulated paths)
path <- ugarchpath(spec.true, n.sim = n, m.sim = 1, rseed = 271) # n.sim length of simulated path; m.sim = number of paths
## Note the difference:
## - ugarchpath(): simulate from a specified model
## - ugarchsim():  simulate from a fitted object
dig <- 2 # digits for rounding
L <- round(fitted(path), digits = dig) # simulated path => actually losses L here

## Fit a GARCH(1, 1) with N(0, 1) innovations
spec.N <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "norm") # without fixed parameters here
fit.N <- ugarchfit(spec.N, data = L) # fit

## Fit a GARCH(1, 1) with standardized t innovations (d.o.f. estimated)
spec.t <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "std") # without fixed parameters here
fit.t <- ugarchfit(spec.t, data = L) # fit

## Compute the corresponding implied VaR estimates
alpha <- 0.99
VaR.N <- round(as.numeric(quantile(fit.N, probs = alpha)), digits = dig)
VaR.t <- round(as.numeric(quantile(fit.t, probs = alpha)), digits = dig)


## b) Solution

## Compute estimates of the forecasting errors; see MFE (2015, Section 9.3.3)
mean(S(VaR.N, l = L, alpha = alpha))
mean(S(VaR.t, l = L, alpha = alpha)) # smaller

## Note: This is a somewhat constructed exercise, for other small choices
##       of 'n' one could equally well find that the score for the wrong model
##       is smaller. Clearly, we would need a much larger 'n' here.


### Exercise 9.18 (Standard methods for market risk) ###########################

## a), b)

## Data preparation
data(EURSTX_const)
time <- c("2000-01-03", "2009-12-31")
S <- EURSTX_const[paste0(time, collapse = "/"), c("BMW.DE", "SIE.DE")] # pick out data (data.frame)
colnames(S) <- c("BMW", "SIE")
S <- na.fill(S, fill = "extend") # fill NAs (see also NA_plot(S))

## Plot the time series to check if anything is 'suspicious'
## BMW
plot.zoo(S[,"BMW"], xlab = "Time t", ylab = expression(BMW~stock~price~S[t]))
## SIE
plot.zoo(S[,"SIE"], xlab = "Time t", ylab = expression(SIE~stock~price~S[t]))

## Compute and plot the risk-factor changes
X <- returns(S) # compute log-returns
## BMW
plot.zoo(X[,"BMW"], xlab = "Time t", ylab = expression(BMW~log-returns~X[t]))
## SIE
plot.zoo(X[,"SIE"], xlab = "Time t", ylab = expression(SIE~log-returns~X[t]))
## Jointly
X <- as.matrix(X)
plot(X, cex = 0.4)


## c)

##' @title Compute the Loss Operator
##' @param x matrix of risk-factor changes
##' @param weights weights tilde{w}_{t,j} = lambda_j * S_{t,j}
##' @return losses -sum_{j=1}^d(tilde{w}_{t,j}(exp(X_{t+1,j})-1))
##' @author Marius Hofert
loss_operator <- function(x, weights)
    -rowSums(expm1(x) * matrix(weights, nrow = nrow(x), ncol = length(weights), byrow = TRUE))

##' @title Estimate VaR and ES
##' @param S stock data, an (n, d)-matrix
##' @param lambda number of shares of each stock
##' @param alpha confidence level for VaR and ES
##' @param method character string specifying the estimator
##' @param ... additional arguments passed to the various methods
##' @return list containing the estimated risk measures VaR and ES, and
##'         possibly other results (depending on the estimator)
##' @author Marius Hofert
risk_measure <- function(S, lambda, alpha,
                         method = c("Var.Cov", "hist.sim", "MC.N", "MC.t", "POT"),
                         ...)
{
    ## Input checks and conversions
    if(!is.matrix(S)) S <- rbind(S, deparse.level = 0L) # guarantees that ncol() works
    stopifnot(0 < alpha, alpha < 1, # check whether alpha is in (0,1)
              length(lambda) == ncol(S), lambda > 0) # check length and sign of lambda
    method <- match.arg(method) # match correct method if not fully provided

    ## Ingredients required for *all* methods
    X <- returns(S) # compute risk-factor changes
    if(!length(X)) stop("'S' should have more than just one line") # check
    S. <- as.numeric(tail(S, n = 1)) # pick out last available stock prices ("today")
    w. <- lambda * S. # weights tilde{w} today

    ## Method switch (now consider the various methods)
    switch(method,
           "Var.Cov" = { # variance-covariance method
               ## Estimate a multivariate normal distribution
               mu <- colMeans(X) # estimate the mean vector mu
               Sigma <- var(X) # estimate the covariance matrix Sigma
               ## Compute the implied parameters of the normal distribution of the linearized loss
               L.Delta.mean <- -sum(w. * mu) # mean of the approx. normal df of the linearized loss L^{Delta}
               L.Delta.sd <- sqrt(t(w.) %*% Sigma %*% w.) # corresponding standard deviation
               ## Compute VaR and ES and return
               qa <- qnorm(alpha)
               list(VaR = L.Delta.mean + L.Delta.sd * qa,
                    ES  = L.Delta.mean + L.Delta.sd * dnorm(qa) / (1-alpha),
               ## => We could just return a vector here, but for other methods,
               ##    we might want to return additional auxiliary results,
               ##    and we should *always* return similar objects (here: lists)
               ## Additional quantities returned here
                    mu    = mu, # fitted mean vector
                    Sigma = Sigma) # fitted covariance matrix
           },
           "hist.sim" = { # historical simulation method
               ## Using nonparametrically estimated risk measures
               L <- loss_operator(X, weights = w.) # compute historical losses
               ## Nonparametrically estimate VaR and ES and return
               list(VaR = VaR_np(L, alpha),
                    ES  =  ES_np(L, alpha))
           },
           "MC.N" = { # Monte Carlo based on a fitted multivariate normal
               stopifnot(hasArg(N)) # check if the number 'N' of MC replications has been provided (via '...')
               N <- list(...)$N # pick out N from '...'
               mu <- colMeans(X) # estimate the mean vector mu
               Sigma  <- var(X) # estimate the covariance matrix Sigma
               X. <- rNorm(N, loc = mu, scale = Sigma) # simulate risk-factor changes
               L <- loss_operator(X., weights = w.) # compute corresponding (simulated) losses
               ## Compute VaR and ES and return
               list(VaR = VaR_np(L, alpha), # nonparametrically estimate VaR
                    ES  =  ES_np(L, alpha), # nonparametrically estimate ES
                    ## Additional quantities returned here
                    mu    = mu, # fitted mean vector
                    Sigma = Sigma) # fitted covariance matrix
           },
           "MC.t" = { # Monte Carlo based on a fitted multivariate t
               stopifnot(hasArg(N)) # check if the number 'N' of MC replications has been provided (via '...')
               N <- list(...)$N # pick out N from '...'
               fit <- fitStudent(X) # fit a multivariate t distribution
               X. <- rStudent(N, df = fit$df, loc = fit$loc, scale = fit$scale) # simulate risk-factor changes
               L <- loss_operator(X., weights = w.) # compute corresponding (simulated) losses
               ## Compute VaR and ES and return
               list(VaR = VaR_np(L, alpha), # nonparametrically estimate VaR
                    ES  =  ES_np(L, alpha), # nonparametrically estimate ES
                    ## Additional quantities returned here
                    mu    = fit$loc, # fitted location vector
                    sigma = fit$scale, # fitted dispersion matrix
                    df    = fit$df) # fitted degrees of freedom
           },
           "POT" = { # simulate losses from a fitted Generalized Pareto distribution (GPD); this is underlying the peaks-over-threshold (POT) method
               stopifnot(hasArg(q)) # check if the quantile-threshold 'q' has been provided
               L. <- loss_operator(X, weights = w.) # historical losses
               u <- quantile(L., probs = list(...)$q, names = FALSE) # determine the threshold as the q-quantile of the historical losses
               excess <- L.[L. > u] - u
               fit <- fit_GPD_MLE(excess) # fit a GPD to the excesses
               xi <- fit$par[["shape"]] # fitted xi
               beta <- fit$par[["scale"]] # fitted beta
               if(xi <= 0) stop("Risk measures only implemented for xi > 0.")
               ## Now compute semi-parametric VaR and ES estimates
               ## G_{xi,beta}(x) = 1-(1+xi*x/beta)^{-1/xi} if xi != 0
               Fbu <- length(excess) / length(L.) # number of excesses / number of losses = N_u / n
               VaR <- u + (beta/xi)*(((1-alpha)/Fbu)^(-xi)-1) # see MFE (2015, Section 5.2.3)
               ES <- (VaR + beta-xi*u) / (1-xi) # see MFE (2015, Section 5.2.3)
               if(xi >= 1) ES <- Inf # adjust to be Inf if xi >= 1 (i.e., ES < 0); see Coles (2001, p. 79)
               ## Return
               list(VaR = VaR, # parametrically estimate VaR
                    ES  = ES, # parametrically estimate ES
                    ## Additional quantities returned here
                    xi  = xi, # fitted xi
                    beta = beta, # fitted beta
                    converged = fit$converged, # did the fitting algorithm converge?
                    u = u, # threshold
                    excess = excess) # excesses over u
           },
           stop("Wrong 'method'"))
}


## d)

lambda <- c(1, 10) # (example) number of shares of the two stocks
alpha <- 0.99 # confidence levels for computing the risk measures
N <- 1e4 # Monte Carlo sample size

## Estimate VaR and ES with the various methods
set.seed(271) # set a seed so that all simulation results are reproducible; see ?set.seed
var.cov  <- risk_measure(S, lambda = lambda, alpha = alpha, method = "Var.Cov")
hist.sim <- risk_measure(S, lambda = lambda, alpha = alpha, method = "hist.sim")
MC.N     <- risk_measure(S, lambda = lambda, alpha = alpha, method = "MC.N", N = N)
POT      <- risk_measure(S, lambda = lambda, alpha = alpha, method = "POT",  N = N, q = 0.9)
MC.t     <- risk_measure(S, lambda = lambda, alpha = alpha, method = "MC.t", N = N)

## Pick out VaR and ES for all methods
rm <- rbind("MC (normal)"    = unlist(MC.N[c("VaR", "ES")]),
            "Var.-cov."      = unlist(var.cov[c("VaR", "ES")]),
            "Hist. sim."     = unlist(hist.sim[c("VaR", "ES")]),
            "POT"            = unlist(POT [c("VaR", "ES")]),
            "MC (Student t)" = unlist(MC.t[c("VaR", "ES")]))
round(rm, digits = 4)

## Visual model assessment with a Q-Q plot
excess <- POT$excess # excesses over the threshold
xi <- POT$xi # estimated xi
beta <- POT$beta # estimated beta
qq_plot(excess, FUN = function(p) qGPD(p, shape = xi, scale = beta)) # fine
mtext(substitute("Q-Q plot for the fitted GPD("*xi.*", "*beta.*") distribution",
                 list(xi. = round(xi, 2), beta. = round(beta, 2))),
      side = 4, line = 1, adj = 0)


## e)

## Compute historical losses (for creating a histogram); see risk.measure()
S. <- as.numeric(tail(S, n = 1)) # pick out last available stock prices ("today")
w. <- lambda * S. # weights tilde{w}
L <- loss_operator(X, weights = w.) # historical losses
summary(L) # get important statistics about the losses

## Plot a histogram of the losses including the VaR and ES estimates
hist(L, breaks = "Scott", probability = TRUE, xlim = c(0, max(L, rm)), main = "",
     xlab = substitute("Losses L > 0 from"~sd~"to"~ed~"with"~
                                         widehat(VaR)[a] <= widehat(ES)[a],
                       list(sd = time[1], ed = time[2], a = alpha)), col = "gray90"); box() # histogram
lty <- c(3, 2, 1, 4, 5)
lwd <- c(1.6, 2, 1, 1.2, 1.2)
for(k in seq_len(nrow(rm)))
    abline(v = rm[k,], lty = lty[k], lwd = lwd[k]) # colored vertical lines indicating VaR and ES
legend("topright", bty = "n", inset = 0.02, lty = lty, lwd = lwd, legend = rownames(rm),
       title = as.expression(substitute(widehat(VaR)[a] <= widehat(ES)[a], list(a = alpha)))) # legend

## From this figure, we can make the following observations:
## 1) The variance-covariance method and the Monte Carlo method based on a
##    fitted multivariate normal distribution lead to similar results (they
##    both assume multivariate normal risk-factor changes but differ in the
##    computation of the loss distribution (analytical vs empirical)).
## 2) Both underestimate VaR_0.99 and ES_0.99 as estimated by historical simulation.
## 3) The historical simulation method implies that the loss distribution
##    is more heavy-tailed. This is captured quite well by the POT method and
##    (possibly overestimated by) the Monte Carlo method based on a fitted
##    multivariate t distribution (the fitted degrees of freedom are 2.3880).
## 4) The EVT-based approach leads results comparable to the historical
##    simulation method.
## 5) It is overall reassuring that several methods (historical simulation, POT
##    method and, with slight departure for ES_0.99, MC for a Student t)
##    lead to similar results. Within this range, one would then determine an
##    adequate risk capital estimate.


### Exercise 9.19 (Standard methods for an option position) ####################

## a)

## Option characteristics
tau.0 <- 0 # tau_t for t = 0 (corresponds to 2009-06-04 when option was sold)
S.0 <- 942.46 # S&P 500 value S_t at t = 0
sig.0 <- 30.18 / 100 # volatility sigma_t at t = 0 (was quoted in % in the exercise)
r <- 0.01 # interest rate per annum
K <- 1050 # strike price
T <- 5 # maturity in years (so 2014-06-04)

## Solve 0 = lambda * S_0 - P^{BS}(tau_0, S_0; r, sigma_0, K, T) w.r.t. lambda
(lambda <- Black_Scholes(tau.0, S = S.0, r = r, sigma = sig.0, K = K, T = T,
                         type = "put") / S.0) # lambda


## c)

## Timeline:
## - Around 2009-06-04: Option was sold
## - 2006-06-14--2010-06-04: Available data window => the corresponding risk factor changes
##                           are used to estimate one-day VaR and ES.
## - End of 2010-06-04 ('today') = 't': Inputs for Black--Scholes formula:
##   + tau_t   = 1 (= 250/250)
##   + S_t     = 1064.88
##   + sigma_t = 35.48
## - Maturity: T = 5 (in years)

## Data preparation
data(SP500) # stock index (S_t)
data(VIX) # volatility (sigma_t) in %
time <- c("2006-06-14", "2010-06-04") # ~= 4 years of risk factor changes
horizon <- paste0(time, collapse = "/")
dat <- merge(SP500[horizon], VIX[horizon]/100)
colnames(dat) <- c("S", "sigma")

## Checks
stopifnot(all(!is.na(dat))) # no missing data
stopifnot(round(dat["2009-06-04",], 4) == c(S.0, sig.0)) # ... as claimed above in a)
stopifnot(round(dat["2010-06-04",], 4) == c(1064.88, 0.3548)) # ... as claimed in c)
stopifnot(nrow(dat) == 1001) # ... as claimed in c); 1001 => 1000 = 4 * 250 risk factor changes


### Compute and plot the risk factor changes ###################################

X <- returns(dat, method = c("logarithmic", "diff")) # risk factor changes
colnames(X) <- c("S", "sigma")
stopifnot(nrow(X) == 1000)
plot.zoo(X[,"S"], xlab = "Time t", ylab = expression("S&P 500 log-returns"~X[list(t,1)]))
plot.zoo(X[,"sigma"], xlab = "Time t", ylab = expression("Volatility differences"~X[list(t,2)]))
X <- as.matrix(X)
plot(X, cex = 0.4)


### Time point t (= 'today' = 2010-06-04) and implied ingredients ##############

## Time point t
tau.t <- 1 # for later
t <- 250 # corresponds to tau.t = 1, so 2010-06-04

## SP500 and volatility at t
S.t   <- as.numeric(dat["2010-06-04", "S"]) # S_t
sig.t <- as.numeric(dat["2010-06-04", "sigma"]) # sigma_t

## Compute Greeks at t
Greeks <- Black_Scholes_Greeks(tau.t, S = S.t, r = r, sigma = sig.t, K = K, T = T,
                               type = "put")

## Compute ingredients of approximations (see solution of b))
Delta.t <- 1/250 # time step in years (one-day horizon)
g.tau <- Greeks[["theta"]] * Delta.t # see solution b)
delta <- c((lambda - Greeks[["delta"]]) * S.t, -Greeks[["vega"]])
Gamma <- matrix(c(-S.t^2 * Greeks[["gamma"]] + (lambda - Greeks[["delta"]]) * S.t,
                  -Greeks[["vanna"]] * S.t, -Greeks[["vanna"]] * S.t,
                  -Greeks[["vomma"]]),
                ncol = 2, byrow = TRUE)


### i) Variance-covariance method ##############################################

## Estimate a multivariate normal distribution
mu <- colMeans(X) # estimate the mean vector mu
Sigma <- var(X) # estimate the covariance matrix Sigma

## Compute the implied parameters of the normal distribution of the linearized loss
## (see the linearized loss operator in the solution of b))
L.Delta.mean <- -(g.tau * Delta.t + sum(delta * mu)) # mean of the approximate normal distribution
L.Delta.sd <- sqrt(t(delta) %*% Sigma %*% delta) # corresponding standard deviation

## Compute implied risk measures
alpha <- 0.99 # confidence level
qa <- qnorm(alpha)
(VaR.VC <- L.Delta.mean + L.Delta.sd * qa) # VaR_alpha
(ES.VC  <- L.Delta.mean + L.Delta.sd * dnorm(qa) / (1-alpha)) # ES_alpha


### ii) Historical simulation with full revaluation ############################

## V_t (t = 'today', so tau.t = 1)
P.BS.t <- Black_Scholes(tau.t, S = S.t, r = r, sigma = sig.t, K = K, T = T, type = "put")
V.t <- lambda * S.t - P.BS.t # V_t

## V_s (s = t+1,...,T based on the 4 years of data in X)
S.s <- exp(log(S.t) + X[,"S"]) # stock prices implied by X
## Derivation: Z_t = log(S_t), X_{t+1} = Z_{t+1} - Z_t
##             => S_{t+1} = e^{Z_{t+1}} = e^{Z_t + X_{t+1}} = e^{log(S_t) + X_{t+1}}
sig.s <- sig.t + X[,"sigma"] # volatilities implied by X
## Derivation: Z_t = sigma_t, X_{t+1} = Z_{t+1} - Z_t
##             => sigma_{t+1} = Z_{t+1} = Z_t + X_{t+1} = sigma_t + X_{t+1}
P.BS.s <- Black_Scholes(Delta.t, S = S.s, r = r, sigma = sig.s, K = K, T = T, type = "put")
V.s <- lambda * S.s - P.BS.s # V_s for s = t+1,...,T

## Losses and corresponding VaR and ES
L.HS.full <- -(V.s - V.t) # losses (implicit loss operator)
(VaR.HS.full <- VaR_np(L.HS.full, level = alpha)) # VaR_alpha
(ES.HS.full  <- ES_np (L.HS.full, level = alpha)) # ES_alpha


### iii) Historical simulation with delta-gamma approximation ##################

## Losses and corresponding VaR and ES
L.HS.DG <- -(g.tau * Delta.t + X %*% delta + 0.5 * rowSums(X %*% Gamma * X)) # (see the code of mahalanobis())
(VaR.HS.DG <- VaR_np(L.HS.DG, level = alpha)) # VaR_alpha
(ES.HS.DG  <- ES_np (L.HS.DG, level = alpha)) # ES_alpha


### iv) Univariate dynamic historical simulation as in MFE (2015, Section 9.2.4)

## We implement the following steps which are listed under (2) in MFE (2015,
## Section 9.2.4); others are possible, too.
## 1) Take the losses L_s, s = 1,...,t from the loss operator and assume they
##    are of the form L_s = mu_s + sigma_s Z_s for SWN(0,1) Z_s
## 2) Fit one ARMA-GARCH model to obtain estimates hat(mu)_s, hat(sigma)_s
##    of mu_s, sigma_s
## 3) Predict mu_{t+1} and sigma_{t+1}
## 4) Compute VaR and ES based on the fitted innovation distribution from the
##    extracted Z_s => hat(VaR), hat(ES)
## 5) The dynamic VaR and ES estimates are thus mu_{t+1} + sigma_{t+1} * hat(VaR)
##    and mu_{t+1} + sigma_{t+1} * hat(ES)

## 1) Historical simulation losses
L <- L.HS.full

## 2) Fit an ARMA(1,1)-GARCH(1,1) model with standardized t innovations to L
##    Specification
mean.model <- list(armaOrder = c(1, 1))
var.model <- list(model = "sGARCH", garchOrder = c(1, 1))
spec <- ugarchspec(var.model, mean.model = mean.model, distribution.model = "std")
##    Fit & extract
fit <- ugarchfit(spec, data = L) # fit
mu <- fitted(fit) # fitted hat(mu)_t
sig <- sigma(fit) # fitted hat(sigma)_t
Z <- residuals(fit, standardize = TRUE) # standardized residuals
##    Check
nu <- fit@fit$coef[["shape"]] # extract (fitted) d.o.f. nu
qq_plot(Z, FUN = function(p) sqrt((nu-2)/nu) * qt(p, df = nu),
        main = substitute("Q-Q plot of ("*Z[t]*") against a standardized"~italic(t)[nu.],
                          list(nu. = round(nu, 2))))
##    => The fit is not good, but we keep this 'broad brush approach' model here

## 3) Predict mu_{t+1}, sigma_{t+1}
pred <- ugarchforecast(fit, n.ahead = 1)
mu.tp1 <- as.numeric(fitted(pred))
sig.tp1 <- as.numeric(sigma(pred))

## 4) Compute VaR and ES based on the fitted innovation distribution
VaR.t <- VaR_t(alpha, loc = 0, scale = sqrt((nu-2)/nu), df = nu)
ES.t  <- ES_t (alpha, loc = 0, scale = sqrt((nu-2)/nu), df = nu)

## 5) Univariate dynamic VaR and ES estimates
(VaR.DHS <- mu.tp1 + sig.tp1 * VaR.t)
(ES.DHS  <- mu.tp1 + sig.tp1 * ES.t)
stopifnot(all.equal(VaR.DHS, as.numeric(quantile(pred, probs = alpha)))) # sanity check with rugarch


### Comparison of VaR and ES for the various methods ###########################

rm <- rbind("Var.-cov." = c(VaR.VC, ES.VC),
            "Hist. sim. full reval" = c(VaR.HS.full, ES.HS.full),
            "Hist. sim. Delta Gamma" = c(VaR.HS.DG, ES.HS.DG),
            "Univ. dynamic hist. sim." = c(VaR.DHS, ES.DHS))
colnames(rm) <- paste(c("VaR", "ES"), alpha, sep = "_")
rm
## Note:
## It is not clear what method is best. We would have to repeat the exercise over
## a period of time and do some backtesting; probably the var-cov method
## underestimates and the dynamic historical simulation method overestimates.


### Exericse 9.21 (PCA factor model for US zero-coupon bond yields) ############

## a)

## Load the data, pick out the components/times we work with and fill gaps
data(ZCB_USD)
yields <- ZCB_USD['2002-01-02/2011-12-30'] # in %
X <- returns(yields, method = "diff") # risk-factor changes


## b)

## PCA
PCA <- prcomp(X) # see also summary(PCA)
Gamma <- PCA$rotation # principal axes (jth column is orthonormal eigenvector of cov(X) corresponding to jth largest eigenvalue) or 'loadings'
mu <- PCA$center # estimated centers (all close to 0)
Y <- PCA$x # estimated principal components of X or 'scores'; (2012, 10)-matrix
var <- PCA$sdev^2 # explained variances per principal component
## Note: var equals the sorted eigenvalues of Cov(X) since diag(<sorted sigma^2>)
##       = Cov(Y) = Cov(Gamma^T (X - mu)) = Gamma^T Cov(X) Gamma = diag(<sorted lambda>)

## Proportion of variability explained by the first three principal components
prop <- cumsum(var)/sum(var)
npr <- 3
prop[npr] # => ~= 98.70%


## c)

## Pick out the first npr-many principal axes/loadings and plot each
## of them over time
Gamma. <- Gamma[,seq_len(npr)] # first npr-many principal axes/loadings
maturities <- seq_len(ncol(X))
plot(NA, xlim = range(maturities), ylim = range(Gamma.),
     xlab = "Time to maturity (years)",
     ylab = "Principal axes value (loading)")
abline(h = 0)
for(j in seq_len(npr))
    lines(maturities, Gamma.[,j], lty = j+1)
legend("top", bty = "n", lty = 1:(npr+1),
       legend = c("Reference line", paste("PC", seq_len(npr))))
## After mirroring the curves corresponding to the first two principal components
## at y = 0, the curves roughly share their shape with MFE (2015, Figure 9.3).


## d)

## Build time series of the first npr-many principal components of X and plot them
Y. <- xts(Y[,seq_len(npr)], time(X)) # factor series
plot.zoo(Y., type = "h", xlab = "Time", ylab = paste("Component", seq_len(npr)),
         main = "Principal components of X")
