## By Marius Hofert

## Non-parametrically estimating value-at-risk (VaR_alpha) and expected shortfall
## (ES_alpha) and providing bootstrapped estimators and confidence intervals (as
## functions in alpha); also estimate Var(hat{VaR}_alpha) and Var(hat{ES}_alpha)

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
##       ES_alpha(X) = int_alpha^1 F^-(u) du / (1-alpha)
##
##    (that is, the u-quantile integrated over all u in [alpha, 1] divided by
##    the length of the integration interval 1-alpha). A substitution leads to
##
##       ES_alpha(X) = int_{F^-(alpha)}^Inf x dF(x) / (1-alpha)
##                   = int_{-Inf}^Inf x * I{x > F^-(alpha)} dF(x) / (1-alpha)
##
##    where I{} denotes the indicator function of the given event.
##    We can approximate this by
##
##       hat{ES}_alpha(X) = sum_{i=1}^n x_i * I{x_i > F^-(alpha)} / (1-alpha)
##
##    which we will use as an estimator below. Note that we actually also
##    estimate the 'F^-(alpha)' in this formula (that's typically the case
##    when F^-(alpha) is unknown).


### 0 Setup ####################################################################

n <- 1300 # sample size (~= 5y of daily data)
B <- 1000 # number of bootstrap replications (= number or realizations of VaR, ES)
th <- 2 # parameter of the Pareto distribution


### 1 Auxiliary functions ######################################################

##' @title Quantile Function of F(x) = 1-(1+x)^{-theta}
##' @param p probability (in [0,1])
##' @param theta Pareto distribution parameter
##' @return p-quantile of the Pareto distribution with parameter theta
##' @author Marius Hofert
qPar <- function(p, theta) (1-p)^(-1/theta) - 1

##' @title Valut-at-Risk for a Par(theta) Distribution
##' @param alpha confidence level
##' @param theta Pareto parameter
##' @return Theoretical value of VaR_alpha(L)
##' @author Marius Hofert
##' @note This is just a convenience wrapper for qPar()
VaR_Par <- function(alpha, theta) qPar(alpha, theta)

##' @title Expected Shortfall for a Par(theta) Distribution
##' @param alpha confidence level
##' @param theta Pareto parameter
##' @return Theoretical value of ES_alpha(L)
##' @author Marius Hofert
ES_Par <- function(alpha, theta)
{
    stopifnot(theta > 1)
    (theta/(theta-1)) * (1-alpha)^(-1/theta) - 1
}

##' @title Nonparametric VaR Estimator
##' @param x losses L
##' @param alpha confidence level
##' @param type 'type' used (1 = inverse of empirical df)
##' @return Nonparametric VaR_alpha estimate
##' @author Marius Hofert
##' @note - Vectorized in x and alpha
##'       - Estimate VaR_alpha for different alpha based on the *same* data x
##'         => less variance.
##'       - We use type=1 here as for sufficiently large alpha, type=7
##'         (quantile()'s default) would interpolate between the two largest
##'         losses and thus return a(n even) smaller VaR_alpha estimate.
VaR <- function(x, alpha, type=1)
    quantile(x, probs=alpha, names=FALSE, type=type) # vectorized in x and alpha

##' @title Nonparametric Expected Shortfall Estimator
##' @param x losses L
##' @param alpha confidence level
##' @return Nonparametric ES_alpha estimate
##' @author Marius Hofert
##' @note - Vectorized in x and alpha
##'       - Estimate ES_alpha for different alpha based on the *same* data x
##'         => less variance.
##'       - If F_L is continuous, then ES_alpha(L) = E(L|L > VaR_alpha(L))
##'         from which we get the below estimator based on the SLLN
##'       - ">" : Mathematically correct for discrete dfs, but
##'               produces NaN for alpha > 1-1/n (=> F^-(alpha) = x_{(n)} but
##'               there is no loss beyond x_{(n)}) because mean(numeric(0)) = NaN
##'         ">=": mean() Will always include the largest loss (so no NaN appears),
##'               but might be computed just based on this one loss.
ES <- function(x, alpha, method=c(">", ">="))
{
    stopifnot(0 < alpha, alpha < 1)
    method <- match.arg(method)
    VaR <- VaR(x, alpha=alpha) # length(alpha)-vector
    vapply(VaR, function(v) { # v = VaR value for one alpha
        ind <- if(method == ">") x > v else x >= v
        num <- sum(ind)
        if(num == 0) {
            warning("No loss ",method," VaR; NaN returned instead")
        } else if(num == 1){
            warning("Only ",num," loss ",method," VaR")
        } else if(num <= 5) {
            warning("Only ",num," losses ",method," VaR")
        }
        mean(x[ind])
    }, NA_real_)
}

##' @title Empirically Estimated Confidence Intervals (default: 95%)
##' @param x losses L
##' @param beta significance level
##' @return Estimated (lower, upper) beta confidence interval
##' @author Marius Hofert
CI <- function(x, beta=0.05, na.rm=FALSE)
    quantile(x, probs=c(beta/2, 1-beta/2), na.rm=na.rm, names=FALSE)

##' @title Bootstrap the Nonparametric VaR or ES Estimator for all alpha
##' @param x losses L
##' @param B number of bootstrap replications
##' @param alpha confidence level
##' @param method risk measure used
##' @return (length(alpha), B)-matrix where the bth column contains the risk
##'         measure estimates evaluated at each alpha based on the bth bootstrap
##'         sample of the losses
##' @author Marius Hofert
##' @note - Vectorized in x and alpha
bootstrap <- function(x, B, alpha, method=c("VaR", "ES"))
{
    ## Checks
    stopifnot(is.numeric(x), (n <- length(x)) >= 1, B >= 1) # sanity checks
    method <- match.arg(method) # check and match 'method'

    ## Bootstrap losses (randomly sampling them B times with replacement)
    x.boot <- matrix(sample(x, size=n*B, replace=TRUE), ncol=B) # (n, B)-matrix of bootstrap samples

    ## For each of the B samples of length n, estimate VaR or ES (for all alpha)
    apply(x.boot, 2, if(method=="VaR") VaR else ES, alpha=alpha) # (length(alpha), B)-matrix
}


### 2 Simulate losses and nonparametrically estimate VaR and ES ################

## Simulate losses
## (as we don't have real ones and we want to investigate the estimators)
set.seed(271) # set a seed (for reproducibility)
L <- qPar(runif(n), theta=th) # simulate losses with the 'inversion method'

## Nonparametrically estimates of VaR and ES
alpha <- 0.99
(VaR. <- VaR(L, alpha=alpha))
(ES.  <-  ES(L, alpha=alpha))
## Just the numbers won't tell us much

## More interesting: behavior in alpha
alpha <- 1-1/10^seq(0.5, 5, by=0.05) # alphas we investigate (concentrated near 1)
stopifnot(0 < alpha, alpha < 1)
VaR. <- VaR(L, alpha=alpha) # estimate VaR_alpha for all alpha
ES.  <-  ES(L, alpha=alpha) # estimate ES_alpha for all alpha

## Plot with logarithmic y-axis and true VaR_alpha and ES_alpha values
VaR.Par. <- VaR_Par(alpha, theta=th) # theoretical VaR_alpha values
ES.Par.  <-  ES_Par(alpha, theta=th) # theoretical ES_alpha values
ran <- range(VaR., ES., VaR.Par., ES.Par., na.rm=TRUE)
plot(alpha, VaR., type="l", ylim=ran, log="y", col="royalblue3",
     xlab=expression(alpha),
     ylab=expression("True"~VaR[alpha]~"and"~ES[alpha]~"and estimates"))
lines(alpha, ES., type="l", col="maroon3")
lines(alpha, VaR.Par., type="l", lty=2, col="royalblue3")
lines(alpha, ES.Par., type="l", lty=2, col="maroon3")
legend("topleft", bty="n", y.intersp=1.2, lty=c(1,2,1,2),
       col=rep(c("maroon3", "royalblue3"), each=2),
       legend=c(expression(widehat(ES)[alpha]), expression(ES[alpha]),
                expression(widehat(VaR)[alpha]), expression(VaR[alpha])))
## => We already see that ES_alpha is more difficult to estimate than VaR_alpha

## Plot with logarithmic y- and x-axis and in 1-alpha (this is why we chose
## an exponential sequence (powers of 10) in alpha concentrated near 1,
## so that the x-axis points are equidistant when plotted in log-scale)
plot(1-alpha, VaR., type="l", ylim=ran, log="xy", col="royalblue3",
     xlab=expression(1-alpha),
     ylab=expression("True"~VaR[alpha]~"and"~ES[alpha]~"and estimates"))
lines(1-alpha, ES., type="l", col="maroon3")
lines(1-alpha, VaR.Par., type="l", lty=2, col="royalblue3")
lines(1-alpha, ES.Par., type="l", lty=2, col="maroon3")
legend("topright", bty="n", y.intersp=1.2, lty=c(1,2,1,2),
       col=rep(c("maroon3", "royalblue3"), each=2),
       legend=c(expression(widehat(ES)[alpha]), expression(ES[alpha]),
                expression(widehat(VaR)[alpha]), expression(VaR[alpha])))
## => Already critical for ES_0.99 and certainly for alpha above 1-1e-3.
##    By using nonparametric estimators we underestimate the risk capital
##    and this although we have comparably large 'n' here!

## Q: Why do the true VaR_alpha and ES_alpha seem 'linear' in alpha for
##    large alpha in log-log scale?
## A: - Linearity in log-log scale means that y is a power function in x:
##      If y = x^beta, then log(y) = beta * log(x), so (log(x), log(y))
##      is a line with slope beta.
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

## Bootstrap the VaR estimator
VaR.boot <- bootstrap(L, B=B, alpha=alpha) # (length(alpha), B)-matrix containing the bootstrapped VaR estimators
any(is.na(VaR.boot)) # no NA or NaN

## Compute statistics
VaR.boot. <- rowMeans(VaR.boot) # bootstrapped mean of the VaR estimator; length(alpha)-vector
VaR.boot.var <- apply(VaR.boot, 1, var) # bootstrapped variance of the VaR estimator; length(alpha)-vector
VaR.boot.CI  <- apply(VaR.boot, 1, CI) # bootstrapped 95% CIs; (2, length(alpha))-matrix

## Bootstrap the ES estimator
system.time(ES.boot <- bootstrap(L, B=B, alpha=alpha, method="ES")) # (length(alpha), B)-matrix containing the bootstrapped ES estimators

## Investigate appearing NaNs (due to too few losses exceeding hat(VaR)_alpha)
isNaN <- is.nan(ES.boot) # => contains some NaN
image(x=alpha, y=seq_len(B), z=isNaN,
      col=c("white", "black"), xlab=expression(alpha), ylab="Bootstrap replication")
mtext("Black: NaN; White: OK", side=4, line=1, adj=0)
percNaN <- 100 * apply(isNaN, 1, mean) # % of NaNs for all alpha
plot(1-alpha, percNaN, type="l", log="x",
     xlab=expression(1-alpha), ylab=expression("% of NaN when estimating"~ES[alpha]))
## => Becomes a serious issue for alpha >= 0.999

## Compute statistics
na.rm <- TRUE # only partially helps (if no data available, still NaN results)
ES.boot. <- rowMeans(ES.boot, na.rm=na.rm) # bootstrapped mean of the ES estimator; length(alpha)-vector
ES.boot.var <- apply(ES.boot, 1, var, na.rm=na.rm) # bootstrapped variance of the ES estimator; length(alpha)-vector
ES.boot.CI  <- apply(ES.boot, 1, CI, na.rm=na.rm) # bootstrapped 95% CIs; (2, length(alpha))-matrix
## => Compute it over the remaining non-NaN values (contains NAs then)


### 4 Plots ####################################################################

## Plot as a function in alpha:
## 1) True VaR_alpha
## 2) Nonparametrically estimated VaR_alpha (based on the sample L)
## 3) The bootstrapped nonparametric estimate of VaR_alpha (mean of the
##    nonparametric VaR_alpha estimators computed from the bootstrap samples)
## 4) The bootstrapped 95% confidence intervals
## 5) The bootstrapped variance of the nonparametric estimator of VaR_alpha
##    (variance of the nonparametric VaR_alpha estimators computed from the
##    bootstrap samples)
## ... and everything for ES_alpha as well
library(sfsmisc)
ran <- range(VaR.Par., # true VaR
             VaR., # nonparametric estimate
             VaR.boot., # bootstrapped estimate (variance sig^2/B for sig^2 = Var(VaR.))
             VaR.boot.CI, # bootstrapped confidence intervals
             VaR.boot.var, # bootstrapped variance
             ES.Par., ES., ES.boot., ES.boot.CI, ES.boot.var, na.rm=TRUE) # same for ES_alpha
## VaR_alpha
plot(1-alpha, VaR.Par., type="l", lty="dotdash", log="xy", xaxt="n", yaxt="n",
     col="royalblue3", ylim=ran, xlab=expression(1-alpha), ylab="") # true VaR_alpha
lines(1-alpha, VaR., lty="dashed", lwd=1.4, col="royalblue3") # nonparametrically estimated VaR_alpha
lines(1-alpha, VaR.boot., lty="solid", col="royalblue3") # bootstrapped nonparametric estimate of VaR_alpha
lines(1-alpha, VaR.boot.CI[1,], lty="dotted", col="royalblue3") # bootstrapped 95% CI
lines(1-alpha, VaR.boot.CI[2,], lty="dotted", col="royalblue3")
lines(1-alpha, VaR.boot.var, lty="4C88C488", col="royalblue3") # bootstrapped Var(hat(VaR)_alpha)
## ES_alpha
lines(1-alpha, ES.Par., lty="dotdash", col="maroon3") # true ES_alpha
lines(1-alpha, ES., lty="dashed", lwd=1.4, col="maroon3") # nonparametrically estimated ES_alpha
lines(1-alpha, ES.boot., lty="solid", col="maroon3") # bootstrapped nonparametric estimate of ES_alpha
lines(1-alpha, ES.boot.CI[1,], lty="dotted", col="maroon3") # bootstrapped 95% CI
lines(1-alpha, ES.boot.CI[2,], lty="dotted", col="maroon3")
lines(1-alpha, ES.boot.var, lty="4C88C488", col="maroon3") # bootstrapped Var(hat(ES)_alpha)
## Misc
eaxis(1) # a nicer exponential y-axis
eaxis(2) # a nicer exponential y-axis
mtext(substitute(B==B.~~"replications of size"~~n==n.~~"from Par("*th.*")",
                 list(B.=B, n.=n, th.=th)), side=4, line=1, adj=0) # secondary y-axis label
legend("bottomleft", bty="n", y.intersp=1.2,
       lty=rep(c("dotdash", "dashed", "solid", "dotted", "4C88C488"), times=2),
       col=rep(c("royalblue3", "maroon3"), each=5),
       legend=c(## VaR_alpha
                expression("True"~VaR[alpha]),
                expression(widehat(VaR)[alpha]),
                expression("Bootstrapped"~~widehat(VaR)[alpha]),
                "Bootstrapped 95% CIs",
                expression("Bootstrapped"~~Var(widehat(VaR)[alpha])),
                ## ES_alpha
                expression("True"~ES[alpha]),
                expression(widehat(ES)[alpha]),
                expression("Bootstrapped"~~widehat(ES)[alpha]),
                "Bootstrapped 95% CIs",
                expression("Bootstrapped"~~Var(widehat(ES)[alpha]))))

## Results:
## - hat(VaR)_alpha < hat(ES)_alpha (clear)
## - hat(VaR)_alpha and hat(ES)_alpha are wrong for large alpha, hat(ES)_alpha
##   even a bit more (a non-parametric estimator cannot adequately capture the
##   tail; this especially applies to hat(ES)_alpha which looks even further
##   in the tail)
## - hat(ES)_alpha (and bootstrapped quantities) are not available for large alpha
##   (there is no observation beyond hat(VaR)_alpha anymore)
## - Var(hat(ES)_alpha) > Var(hat(VaR)_alpha)
##   (hat(ES)_alpha requires information from the *whole* tail)
## - Var(hat(VaR)_alpha) and Var(hat(ES)_alpha) increase in alpha
##   (less data to the right of VaR_alpha => higher variance)
