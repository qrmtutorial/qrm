## By Marius Hofert

## We ...
## - non-parametrically estimate value-at-risk (VaR_alpha) and expected shortfall (ES_alpha);
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

library(qrmtools)

n <- 2500 # sample size (~= 10y of daily data)
B <- 1000 # number of bootstrap replications (= number or realizations of VaR, ES)
th <- 2 # Pareto parameter (true underlying distribution)


### 1 Auxiliary functions ######################################################

##' @title Empirically estimate confidence intervals (default: 95%)
##' @param x The values
##' @param beta The significance level
##' @return Estimated (lower, upper) beta confidence interval
##' @author Marius Hofert
CI <- function(x, beta = 0.05, na.rm = FALSE)
    quantile(x, probs = c(beta/2, 1-beta/2), na.rm = na.rm, names = FALSE)

##' @title Bootstrap the nonparametric VaR or ES estimator for all alpha
##' @param x The vector of losses
##' @param B The number of bootstrap replications
##' @param alpha The confidence level
##' @param method The risk measure used
##' @return (length(alpha), B)-matrix where the bth column contains the estimated
##'         risk measure at each alpha based on the bth bootstrap sample of the
##'         losses
##' @author Marius Hofert
##' @note Vectorized in x and alpha
bootstrap <- function(x, B, alpha, method = c("VaR", "ES"))
{
    stopifnot(is.vector(x), (n <- length(x)) >= 1, B >= 1) # sanity checks
    method <- match.arg(method) # check and match 'method'
    rm <- if(method == "VaR") { # risk measure
        VaR_np
    } else {
        function(x, alpha) ES_np(x, alpha = alpha, verbose = TRUE)
    }
    x.boot <- matrix(sample(x, size = n*B, replace = TRUE), ncol = B) # (n, B)-matrix of bootstrap samples
    apply(x.boot, 2, rm, alpha = alpha) # (length(alpha), B)-matrix
}


### 2 Simulate losses and nonparametrically estimate VaR and ES ################

## Simulate losses (as we don't have real ones and we want to investigate
## the performance of the estimators)
set.seed(271) # set a seed (for reproducibility)
L <- rPar(n, theta = th) # simulate losses with the 'inversion method'


### 2.1 Nonparametric estimates of VaR_alpha and ES_alpha for a fixed alpha ####

alpha <- 0.99
(VaR. <- VaR_np(L, alpha = alpha))
(ES.  <-  ES_np(L, alpha = alpha, verbose = TRUE))
## ... but single numbers don't tell us much
## More interesting: The behavior in alpha


### 2.2 As functions in alpha vs true values ###################################

## Compute the nonparametric VaR_alpha and ES_alpha estimators as functions of alpha
alpha <- 1-1/10^seq(0.5, 5, by = 0.05) # alphas we investigate (concentrated near 1)
stopifnot(0 < alpha, alpha < 1)
VaR. <- VaR_np(L, alpha = alpha) # estimate VaR_alpha for all alpha
ES.  <-  ES_np(L, alpha = alpha, verbose = TRUE) # estimate  ES_alpha for all alpha
warnings()

## True values (known here)
VaR.Par. <- VaR_Par(alpha, theta = th) # theoretical VaR_alpha values
ES.Par.  <-  ES_Par(alpha, theta = th) # theoretical ES_alpha values

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
## => We already see that ES_alpha is more difficult to estimate than VaR_alpha

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
VaR.boot <- bootstrap(L, B = B, alpha = alpha) # (length(alpha), B)-matrix containing the bootstrapped VaR estimators
stopifnot(all(!is.na(VaR.boot))) # no NA or NaN

## Compute statistics
VaR.boot.    <- rowMeans(VaR.boot) # bootstrapped mean of the VaR estimator; length(alpha)-vector
VaR.boot.var <- apply(VaR.boot, 1, var) # bootstrapped variance of the VaR estimator; length(alpha)-vector
VaR.boot.CI  <- apply(VaR.boot, 1, CI) # bootstrapped 95% CIs; (2, length(alpha))-matrix


### 3.2 ES_alpha ###############################################################

## Bootstrap the ES estimator
system.time(ES.boot <- bootstrap(L, B = B, alpha = alpha, method = "ES")) # (length(alpha), B)-matrix containing the bootstrapped ES estimators
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
library(sfsmisc) # for eaxis()
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
## Misc
eaxis(1) # a nicer exponential y-axis
eaxis(2) # a nicer exponential y-axis
mtext(substitute(B == B.~~"replications of size"~~n == n.~~"from Par("*th.*")",
                 list(B. = B, n. = n, th. = th)), side = 4, line = 1, adj = 0) # secondary y-axis label
legend("bottomleft", bty = "n", lwd = c(1, rep(c(2, 1, 1.4, 1), 2)),
       lty = c("solid", rep(c("dashed", "solid", "dotdash", "dotted"), times = 2)),
       col = c("black", rep(c("royalblue3", "maroon3"), each = 4)),
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
                  "Bootstrapped 95% CIs"))

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