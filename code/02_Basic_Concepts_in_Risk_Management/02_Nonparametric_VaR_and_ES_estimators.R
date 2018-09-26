## By Marius Hofert

## We want to...
## - non-parametrically estimate VaR_alpha and ES_alpha;
## - compute bootstrapped estimators;
## - compute bootstrapped confidence intervals;
## - estimate Var(hat{VaR}_alpha) and Var(hat{ES}_alpha);
## ... as functions in alpha

## Note:
## 1) This is an actual problem from the realm of Quantitative Risk Management.
##    Bigger banks and insurance companies are required to compute (on a daily
##    basis) risk measures such as VaR and ES based on loss distributions
##    (according to the Basel II/Solvency II guidelines).
## 2) Value-at-risk is defined by
##
##       VaR_alpha(X) = F^-(alpha) = inf{x in IR: F(x) >= alpha}
##
##    (that is, the alpha-quantile of the df F of X). An estimator for VaR_alpha(X)
##    can be obtained by replacing F by it's empirical distribution function.
##    The corresponding empirical quantile is then simply
##
##       hat{VaR}_alpha(X) = X_{(ceiling(n*alpha))};
##
##    here X_{(k)} denotes the k-th smallest value among X_1,..,X_n and
##    ceiling() denotes the ceiling function.
## 3) Expected shortfall is defined by
##
##       ES_alpha(X) = (1/(1-alpha)) * int_alpha^1 F^-(u) du
##
##    (that is, the u-quantile integrated over all u in [alpha, 1] divided by
##    the length of the integration interval 1-alpha). Under continuity,
##
##       ES_alpha(X) = E(X | X > VaR_alpha(X))
##                   = (1/(1-alpha)) * E(X * I_{X > VaR_alpha(X)})
##
##    where I{} denotes the indicator function of the given event. This can be
##    estimated via
##
##       hat{ES}_alpha(X) = (1/(1-alpha)) * (1/n) * sum(X_i * I_{X_i > hat{VaR}_alpha(X)})
##                        = (1 / E(N)) * sum(X_i * I_{X_i > hat{VaR}_alpha(X)})
##                        = mean(Y_i)
##
##    where E(N) = (1-alpha) * n = P(X > VaR_alpha(X)) * n is the expected number
##    of X_i's exceeding VaR_alpha(X) and Y_i are those X_i which exceed
##    hat{VaR}_alpha(X).


### Setup ######################################################################

library(sfsmisc) # for eaxis()
library(qrmtools)

n <- 2500 # sample size (~= 10y of daily data)
B <- 1000 # number of bootstrap replications (= number or realizations of VaR, ES)


### 1 Auxiliary functions ######################################################

##' @title Empirically estimate confidence intervals (default: 95%)
##' @param x values
##' @param alpha significance level
##' @return estimated (lower, upper) beta confidence interval
##' @author Marius Hofert
CI <- function(x, alpha = 0.05, na.rm = FALSE)
    quantile(x, probs = c(alpha/2, 1-alpha/2), na.rm = na.rm, names = FALSE)

##' @title Bootstrap the nonparametric VaR or ES estimator for all alpha
##' @param x vector of losses
##' @param B number of bootstrap replications
##' @param level confidence level
##' @param method risk measure used
##' @return (length(alpha), B)-matrix where the bth column contains the estimated
##'         risk measure at each alpha based on the bth bootstrap sample of the
##'         losses
##' @author Marius Hofert
##' @note vectorized in x and level
bootstrap <- function(x, B, level, method = c("VaR", "ES"))
{
    stopifnot(is.vector(x), (n <- length(x)) >= 1, B >= 1) # sanity checks
    ## Define the risk measure (as a function of x, level)
    method <- match.arg(method) # check and match 'method'
    rm <- if(method == "VaR") {
        VaR_np # see qrmtools; essentially quantile(, type = 1)
    } else {
        function(x, level) ES_np(x, level = level, verbose = TRUE) # see qrmtools; uses '>' and 'verbose'
    }
    ## Construct the bootstrap samples (by drawing with replacement)
    ## from the underlying empirical distribution function
    x.boot <- matrix(sample(x, size = n * B, replace = TRUE), ncol = B) # (n, B)-matrix
    ## For each bootstrap sample, estimate the risk measure
    apply(x.boot, 2, rm, level = level) # (length(level), B)-matrix
}


### 2 Simulate losses and nonparametrically estimate VaR and ES ################

## Simulate losses (as we don't have real ones and we want to investigate
## the performance of the estimators)
th <- 2 # Pareto parameter (true underlying distribution; *just* infinite Var)
set.seed(271) # set a seed (for reproducibility)
L <- rPar(n, shape = th) # simulate losses with the 'inversion method'
plot(L)


### 2.1 Nonparametric estimates of VaR_alpha and ES_alpha for a fixed alpha ####

alpha <- 0.99
(VaR. <- VaR_np(L, level = alpha))
(ES.  <-  ES_np(L, level = alpha, verbose = TRUE))
abline(h = c(VaR., ES.), lty = 2)
## ... but single numbers don't tell us much
## More interesting: The behavior in alpha


### 2.2 As functions in alpha vs true values ###################################

## Compute the nonparametric VaR_alpha and ES_alpha estimators as functions of alpha
alpha <- 1-10^-seq(0.5, 5, by = 0.05) # alphas we investigate (concentrated near 1)
stopifnot(0 < alpha, alpha < 1)
VaR. <- VaR_np(L, level = alpha) # estimate VaR_alpha for all alpha
ES.  <-  ES_np(L, level = alpha, verbose = TRUE) # estimate ES_alpha for all alpha
if(FALSE)
    warnings()

## True values (known here)
VaR.Par. <- VaR_Par(alpha, shape = th) # theoretical VaR_alpha values
ES.Par.  <-  ES_Par(alpha, shape = th) # theoretical ES_alpha values

## Plot estimates with true VaR_alpha and ES_alpha values
ran <- range(ES.Par., ES., VaR.Par., VaR., na.rm = TRUE)
plot(alpha, ES.Par., type = "l", ylim = ran,
     xlab = expression(alpha),
     ylab = expression("True"~VaR[alpha]*","~ES[alpha]~"and estimates")) # true ES_alpha
lines(alpha, ES., type = "l", col = "maroon3") # ES_alpha estimate
lines(alpha, VaR.Par., type = "l") # true VaR_alpha
lines(alpha, VaR., type = "l", col = "royalblue3") # VaR_alpha estimate
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 3),
       col = c("black", "maroon3", "royalblue3"),
       legend = c(expression(ES[alpha]~"and"~VaR[alpha]),
                  expression(widehat(ES)[alpha]),
                  expression(widehat(VaR)[alpha])))
## => We don't see much

## With logarithmic y-axis
ran <- range(ES.Par., ES., VaR.Par., VaR., na.rm = TRUE)
plot(alpha, ES.Par., type = "l", ylim = ran, log = "y",
     xlab = expression(alpha),
     ylab = expression("True"~VaR[alpha]*","~ES[alpha]~"and estimates")) # true ES_alpha
lines(alpha, ES., type = "l", col = "maroon3") # ES_alpha estimate
lines(alpha, VaR.Par., type = "l") # true VaR_alpha
lines(alpha, VaR., type = "l", col = "royalblue3") # VaR_alpha estimate
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 3),
       col = c("black", "maroon3", "royalblue3"),
       legend = c(expression(ES[alpha]~"and"~VaR[alpha]),
                  expression(widehat(ES)[alpha]),
                  expression(widehat(VaR)[alpha])))
## => We already see that ES_alpha always underestimates its true value
##    (ES_alpha is more difficult to estimate than VaR_alpha)

## With logarithmic x- and y-axis and in 1-alpha
## Note: This is why we chose an exponential sequence (powers of 10) in alpha
##       concentrated near 1, so that the x-axis points are equidistant when
##       plotted in log-scale)
plot(1-alpha, ES.Par., type = "l", ylim = ran, log = "xy",
     xlab = expression(1-alpha),
     ylab = expression("True"~VaR[alpha]*","~ES[alpha]~"and estimates")) # true ES_alpha
lines(1-alpha, ES., type = "l", col = "maroon3") # ES_alpha estimate
lines(1-alpha, VaR.Par., type = "l") # true VaR_alpha
lines(1-alpha, VaR., type = "l", col = "royalblue3") # VaR_alpha estimate
legend("topright", bty = "n", y.intersp = 1.2, lty = rep(1, 3),
       col = c("black", "maroon3", "royalblue3"),
       legend = c(expression(ES[alpha]~"and"~VaR[alpha]),
                  expression(widehat(ES)[alpha]),
                  expression(widehat(VaR)[alpha])))
## => - Critical for ES_0.99; severely critical for alpha ~ 0.999.
##    - By using nonparametric estimators we underestimate the risk
##      capital and this although we have comparably large 'n' here!
##    - Note that VaR_alpha is simply the largest max(L) for alpha
##      sufficiently large.

## Q: Why do the true VaR_alpha and ES_alpha seem 'linear' in alpha for
##    large alpha in log-log scale?
## A: - Linearity in log-log scale means that y is a power function in x:
##
##         y = x^beta   =>   log(y) = beta * log(x),
##
##      so (log(x), log(y)) is a line with slope beta if y = x^beta.
##    - For a Pareto distribution,
##
##         VaR_alpha(L) = (1-alpha)^{-1/theta} - 1,
##          ES_alpha(L) = (theta/(theta-1)) * (1-alpha)^{-1/theta} - 1.
##
##      Omitting the '-1' for large alpha, we see that
##
##         VaR_alpha(L) ~= (1-alpha)^{-1/theta},
##          ES_alpha(L) ~= (theta/(theta-1)) * (1-alpha)^{-1/theta}
##
##      from which it follows that
##
##         (log(1-alpha), log(VaR_alpha(L))) ~= a line with slope -1/theta,
##         (log(1-alpha), log(ES_alpha(L)))  ~= a line with slope -1/theta and
##                                              intercept log(theta/(theta-1))


### 3 Bootstrapping confidence intervals for VaR_alpha, ES_alpha ###############

### 3.1 VaR_alpha ##############################################################

## Bootstrap the VaR estimator
set.seed(271)
VaR.boot <- bootstrap(L, B = B, level = alpha) # (length(alpha), B)-matrix containing the bootstrapped VaR estimators
stopifnot(all(!is.na(VaR.boot))) # no NA or NaN

## Compute statistics
VaR.boot.    <- rowMeans(VaR.boot) # bootstrapped mean of the VaR estimator; length(alpha)-vector
VaR.boot.var <- apply(VaR.boot, 1, var) # bootstrapped variance of the VaR estimator; length(alpha)-vector
VaR.boot.CI  <- apply(VaR.boot, 1, CI) # bootstrapped 95% CIs; (2, length(alpha))-matrix


### 3.2 ES_alpha ###############################################################

## Bootstrap the ES estimator
set.seed(271)
system.time(ES.boot <- bootstrap(L, B = B, level = alpha, method = "ES")) # (length(alpha), B)-matrix containing the bootstrapped ES estimators
if(FALSE)
    warnings()

## Investigate appearing NaNs (due to too few losses exceeding hat(VaR)_alpha)
isNaN <- is.nan(ES.boot) # => contains some NaN
image(x = alpha, y = seq_len(B), z = isNaN, col = c("white", "black"),
      xlab = expression(alpha), ylab = "Bootstrap replication")
mtext("Black: NaN; White: OK", side = 4, line = 1, adj = 0)

## Plot % of NaN
percNaN <- 100 * apply(isNaN, 1, mean) # % of NaNs for all alpha
plot(1-alpha, percNaN, type = "l", log = "x",
     xlab = expression(1-alpha), ylab = expression("% of NaN when estimating"~ES[alpha]))
## => Becomes a serious issue for alpha >= 0.999 (1-alpha <= 10^{-3})

## Compute statistics with NaNs removed
na.rm <- TRUE # only partially helps (if no data available, the result is still NaN)
ES.boot.    <- rowMeans(ES.boot, na.rm = na.rm) # bootstrapped mean of the ES estimator; length(alpha)-vector
ES.boot.var <- apply(ES.boot, 1, var, na.rm = na.rm) # bootstrapped variance of the ES estimator; length(alpha)-vector
ES.boot.CI  <- apply(ES.boot, 1, CI, na.rm = na.rm) # bootstrapped 95% CIs; (2, length(alpha))-matrix
## => Compute it over the remaining non-NaN values (contains NAs then)


### 4 Plots ####################################################################

## Plot as a function in alpha:
## 1) True VaR_alpha and ES_alpha
## 2) The nonparametric estimates hat{VaR}_alpha and hat{ES}_alpha
## 3) The bootstrapped mean of hat{VaR}_alpha and hat{ES}_alpha
## 4) The bootstrapped variance of hat{VaR}_alpha and hat{ES}_alpha
## 5) The bootstrapped 95% confidence intervals for VaR_alpha and ES_alpha
ran <- range(VaR.Par., # true VaR
             VaR., # nonparametric estimate
             VaR.boot., # bootstrapped estimate (variance Var(VaR.)/B)
             VaR.boot.var, # bootstrapped variance
             VaR.boot.CI, # bootstrapped confidence intervals
             ES.Par., ES., ES.boot., ES.boot.var, ES.boot.CI, na.rm = TRUE) # same for ES_alpha
## VaR_alpha
plot(1-alpha, VaR.Par., type = "l", log = "xy", xaxt = "n", yaxt = "n",
     ylim = ran, xlab = expression(1-alpha), ylab = "") # true VaR_alpha
lines(1-alpha, VaR., lty = "dashed", lwd = 2, col = "royalblue3") # nonparametrically estimated VaR_alpha
lines(1-alpha, VaR.boot., lty = "solid", col = "royalblue3") # bootstrapped nonparametric estimate of VaR_alpha
lines(1-alpha, VaR.boot.var, lty = "dotdash", lwd = 1.4, col = "royalblue3") # bootstrapped Var(hat(VaR)_alpha)
lines(1-alpha, VaR.boot.CI[1,], lty = "dotted", col = "royalblue3") # bootstrapped 95% CI
lines(1-alpha, VaR.boot.CI[2,], lty = "dotted", col = "royalblue3")
## ES_alpha
lines(1-alpha, ES.Par.) # true ES_alpha
lines(1-alpha, ES., lty = "dashed", lwd = 2, col = "maroon3") # nonparametrically estimated ES_alpha
lines(1-alpha, ES.boot., lty = "solid", col = "maroon3") # bootstrapped nonparametric estimate of ES_alpha
lines(1-alpha, ES.boot.var, lty = "dotdash", lwd = 1.4, col = "maroon3") # bootstrapped Var(hat(ES)_alpha)
lines(1-alpha, ES.boot.CI[1,], lty = "dotted", col = "maroon3") # bootstrapped 95% CI
lines(1-alpha, ES.boot.CI[2,], lty = "dotted", col = "maroon3")
## Also adding peaks-over-threshold-based estimators (as motivation for Chapter 5)
u <- quantile(L, probs = 0.95, names = FALSE)
L. <- L[L > u] - u # compute the excesses over u
fit <- fit_GPD_MLE(L.) # fit a GPD to the excesses
xi <- fit$par[["shape"]] # fitted shape xi
beta <- fit$par[["scale"]] # fitted scale beta
if(xi <= 0) stop("Risk measures only implemented for xi > 0.")
Fbu <- length(L.) / length(L) # number of excesses / number of losses = N_u / n
VaR.POT <- u + (beta/xi)*(((1-alpha)/Fbu)^(-xi)-1) # see McNeil, Frey, Embrechts (2015, Section 5.2.3)
ES.POT <- (VaR.POT + beta-xi*u) / (1-xi) # see McNeil, Frey, Embrechts (2015, Section 5.2.3)
lines(1-alpha, VaR.POT, lty = "dashed")
lines(1-alpha, ES.POT,  lty = "dashed")
## Misc
eaxis(1) # a nicer exponential y-axis
eaxis(2) # a nicer exponential y-axis
mtext(substitute(B == B.~~"replications of size"~~n == n.~~"from Par("*th.*")",
                 list(B. = B, n. = n, th. = th)), side = 4, line = 1, adj = 0) # secondary y-axis label
legend("bottomleft", bty = "n", lwd = c(1, rep(c(2, 1, 1.4, 1), 2), 1),
       lty = c("solid", rep(c("dashed", "solid", "dotdash", "dotted"), times = 2), "dashed"),
       col = c("black", rep(c("royalblue3", "maroon3"), each = 4), "black"),
       legend = c(## VaR_alpha
                  expression("True"~VaR[alpha]~"and"~ES[alpha]),
                  expression(widehat(VaR)[alpha]),
                  expression("Bootstrapped"~~widehat(VaR)[alpha]),
                  expression("Bootstrapped"~~Var(widehat(VaR)[alpha])),
                  "Bootstrapped 95% CIs",
                  ## ES_alpha
                  expression(widehat(ES)[alpha]),
                  expression("Bootstrapped"~~widehat(ES)[alpha]),
                  expression("Bootstrapped"~~Var(widehat(ES)[alpha])),
                  "Bootstrapped 95% CIs",
                  ## POT
                  expression("POT-based"~VaR[alpha]~"and"~ES[alpha]~"estimates")))

## Results:
## - hat(VaR)_alpha < hat(ES)_alpha (clear)
## - hat(ES)_alpha (and bootstrapped quantities) are not available for large alpha
##   (there is no observation beyond hat(VaR)_alpha anymore)
## - hat(VaR)_alpha and hat(ES)_alpha underestimate VaR_alpha and ES_alpha for
##   large alpha; hat(ES)_alpha even a bit more
##   (a non-parametric estimator cannot adequately capture the tail)
## - Var(hat(VaR)_alpha) and Var(hat(ES)_alpha) (mostly) increase in alpha
##   (less data to the right of VaR_alpha => higher variance)
## - Var(hat(ES)_alpha) > Var(hat(VaR)_alpha)
##   (hat(ES)_alpha requires information from the *whole* tail)
## - For not too extreme alpha, the estimated CI for ES_alpha are larger than
##   those for VaR_alpha
## - The POT-based method (see later) allows more accurate estimates for large alpha
