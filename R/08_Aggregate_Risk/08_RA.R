## By Marius Hofert

## How to compute worst and best VaR with the rearrangement algorithm (RA)
## and the adaptive rearrangement algorithm (ARA)


### Setup ######################################################################

library(copula)
library(qrmdata)
library(qrmtools)


### 1 Auxiliary functions ######################################################

##' @title A Basic Version of the RA
##' @param X The (N,d)-matrix of discretized marginals
##' @param tol The tolerance epsilon; if NULL, columns are rearranged
##'        until each column is oppositely ordered to the sum of all others
##' @return The worst VaR
##' @author Marius Hofert
##' @note Already uses relative tolerances
basic_rearrange_worst_VaR <- function(X, tol = NULL)
{
    N <- nrow(X)
    d <- ncol(X)
    m.rs.old <- min(rowSums(X))
    while (TRUE) {
        Y <- X
        for(j in 1:d)
            Y[,j] <- sort(Y[,j], decreasing = TRUE)[rank(rowSums(Y[,-j, drop = FALSE]), ties.method = "first")]
        Y.rs <- rowSums(Y)
        m.rs.new <- min(Y.rs)
        tol. <- abs((m.rs.new - m.rs.old)/m.rs.old) # relative tolerance
        tol.reached <- if(is.null(tol)) {
            identical(Y, X)
        } else { tol. <= tol }
        if(tol.reached) {
            break
        } else {
            m.rs.old <- m.rs.new
            X <- Y
        }
    }
    min(rowSums(Y))
}


### 2 A first example based on basic_rearrange_worst_VaR() #####################

## Define the parameters of three example margins we consider here
th <- 2.5 # Pareto shape parameter
m <- 10 # mean of the log-normal
v <- 20 # variance of the log-normal
s <- 4 # shape of the gamma underlying the log-gamma
r <- 5 # rate of the gamma underlying the log-gamma

## Define list of marginal dfs
qF <- list(qPar = function(p) (1 - p)^(-1/th) - 1,
           qLN  = function(p) qlnorm(p, meanlog = log(m)-log(1+v/m^2)/2, # mean(log(LN)) = mean(N)
                                          sdlog = sqrt(log(1+v/m^2))), # sd(log(LN)) = sd(N))
           qLG  = function(p) exp(qgamma(p, shape = s, rate = r)))

## Compute worst VaR with basic_rearrange_worst_VaR()
alpha <- 0.99 # confidence level
N <- 2^10 # number of discretization points
p <- alpha + (1-alpha)*0:(N-1)/N # probabilities at which to evaluate the marginal quantile functions
X <- sapply(qF, function(qF.) qF.(p)) # input matrix for computing worst VaR
(worst.VaR.basic <- basic_rearrange_worst_VaR(X)) # rearrange until *all* columns are oppositely ordered
stopifnot(all.equal(worst.VaR.basic,
                    basic_rearrange_worst_VaR(X, tol = 0), tol = 2e-4))
## => Almost the same as rearranging until there is no change anymore ('tol = 0')


### 3 The more sophisticated RA() and ARA() ####################################

## Rearrangement algorithm (RA)
## In comparison to the original RA, qrmtools' RA()...
## - ... uses relative tolerances
## - ... deals with ties
## - ... computes more information
## - ... uses various speed-ups (see the vignette "VaR_bounds")
RA.worst <- RA(alpha, qF = qF, N = N)
str(RA.worst)
(worst.VaR.RA <- mean(RA.worst$bounds)) # worst VaR estimate
stopifnot(all.equal(worst.VaR.RA, worst.VaR.basic, tol = 5e-4)) # comparison with basic_...()

## Adaptive rearrangement algorithm (ARA)
## Besides the above, this...
## - ... chooses 'N' adaptively
## - ... uses the notion of a joint relative tolerance to guarantee that the two
##   bounds (on, say, worst VaR) are close
ARA.worst <- ARA(alpha, qF = qF)
str(ARA.worst) # => a lot of information
(worst.VaR.ARA <- mean(ARA.worst$bounds)) # worst VaR estimate
stopifnot(all.equal(worst.VaR.ARA, worst.VaR.RA, tol = 1e-3)) # comparison with RA()

## ... and for best VaR
RA.best <- RA(alpha, qF = qF, N = N, method = "best.VaR")
ARA.best  <- ARA(alpha, qF = qF, method = "best.VaR")
stopifnot(all.equal(mean(RA.best$bounds), mean(ARA.best$bounds), tol = 5e-3))


### 4 How the underlying rearrange() proceeds ##################################

## See the algorithm 'at work'
A <- matrix(c(1:4, 2*(1:4)-1, 2^(0:3)), ncol = 3) # input matrix representing (dummy) quantiles
rearrange(A, tol = NULL, sample = FALSE, is.sorted = TRUE, trace = TRUE)

## An example for which the algorithm does not find the optimal rearrangement
## (but at least, with our stable sort, it terminates!)
B <- matrix(rep(1:3, 3), ncol = 3) # input matrix representing (dummy) quantiles
rearrange(B, tol = NULL, sample = FALSE, is.sorted = TRUE, trace = TRUE)


### 5 Making the worst- and best-case copulas visible ##########################

## Worst-case copula conditional on all margins being in [alpha, 1]
X.worst <- ARA.worst[["X.rearranged"]]$up # extract rearranged matrix (upper bound)
U.worst <- pobs(X.worst) # compute pseudo-observations
pairs2(U.worst) # approx. sample of a copula leading to worst VaR for our marginal dfs

## Best-case copula conditional on all margins being in [alpha, 1]
X.best <- ARA.best[["X.rearranged"]]$up
U.best <- pobs(X.best)
pairs2(U.best, col = adjustcolor("black", alpha.f = 0.2)) # approx. sample of a copula leading to best VaR for our marginal dfs


### 6 Real-data example ########################################################

## Data we work with
data("SP500_const") # load the constituents data of the S&P 500
stocks <- c("GOOGL", "AAPL", "MSFT") # Google, Apple, Microsoft
time <- c("2006-01-03", "2015-12-31") # time period
S <- SP500_const[paste0(time, collapse = "/"), stocks] # data
stopifnot(all(!is.na(S)))
X <- -returns(S) # -log-returns
## Note: One would need to fit time-series models here (deGARCHing), but
##       we omit that step here and treat the -log-returns as iid.

## Determine thresholds for POT method
mean_excess_plot(X[X[,"GOOGL"] > 0, "GOOGL"])
abline(v = 0.012)
mean_excess_plot(X[X[,"AAPL"]  > 0, "AAPL"])
abline(v = 0.021)
mean_excess_plot(X[X[,"MSFT"]  > 0, "MSFT"])
abline(v = 0.017)
u <- c(0.012, 0.021, 0.017)

## Fit GPDs to the excesses (per margin)
d <- ncol(X)
fit <- lapply(1:d, function(j) fit_GPD_MLE(X[X[,j] > u[j],j] - u[j]))
params <- sapply(fit, function(x) x$par)
colnames(params) <- names(X)

## Estimate threshold exceedance probabilities
p.exceed <- sapply(1:d, function(j) mean(X[,j] > u[j]))
stopifnot(alpha >= max(1-p.exceed)) # check validity (alpha >= F(u))

## Define corresponding quantile functions
qF <- lapply(1:d, function(j) function(p)
    qGPDtail(p, threshold = u[j], p.exceed = p.exceed[j],
             shape = params["shape",j], scale = params["scale",j]))

## Compute worst Value-at-Risk
ARA.worst <- ARA(alpha, qF = qF)
mean(ARA.worst$bounds) # worst VaR estimate
U <- pobs(ARA.worst[["X.rearranged"]]$up)
pairs2(U) # worst VaR dependence
plot(U[,1])
plot(U[,2])
plot(U[,3])

## ... and for best VaR (lots of ties here)
ARA.best  <- ARA(alpha, qF = qF, method = "best.VaR")
mean(ARA.best$bounds) # best VaR estimate
if(FALSE) { # degenerate here (ties problem)
    U <- pobs(ARA.best[["X.rearranged"]]$up)
    pairs2(U)
    plot(U[,1])
    plot(U[,2])
    plot(U[,3])
}

