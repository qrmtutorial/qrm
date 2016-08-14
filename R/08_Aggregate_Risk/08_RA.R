## By Marius Hofert

## How to compute best and worst VaR with the rearrangement algorithm (RA)
## and the adaptive rearrangement algorithm (ARA)


### 0 Setup ####################################################################

library(qrmtools)
library(qrmdata)


### 1 Building the marginal quantile functions #################################

## Select the data we work with
data("SP500_const") # load the constituents data of the S&P 500
stocks <- c("GOOGL", "AAPL", "MSFT") # Google, Apple, Microsoft
time <- c("2006-01-03", "2015-12-31") # time period
S <- SP500_const[paste0(time, collapse = "/"), stocks] # data
stopifnot(all(!is.na(S)))

## Build negative log-returns and corresponding marginal empirical quantile
## functions
X <- -log_returns(S) # -log-returns
qF <- lapply(1:ncol(X), function(j)
    function(p) quantile(X[,j], probs = p, names = FALSE, type = 1))

## Append two parametric marginal quantile functions (as obtained by estimation,
## for example, if the number of losses is too small to use empirical dfs)
qF <- c(qF, function(p) qPar(p, theta = 2),
            function(p) qlnorm(p, meanlog = 5, sdlog = 0.5)) # mean(log(LN)) = mean(N), sd(log(LN)) = sd(N)


### 2 A basic implementation of the RA #########################################

##' @title A basic version of the RA
##' @param X The (N,d)-matrix of discretized marginals
##' @param tol The tolerance epsilon; if NULL, columns are rearranged
##'        until each column is oppositely ordered to the sum of all others
##' @return The worst VaR
##' @author Marius Hofert
##' @note Already uses relative errors
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
        tol. <- abs((m.rs.new - m.rs.old)/m.rs.old)
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

## Example
alpha <- 0.99 # confidence level
N <- 2^10 # number of discretization points
set.seed(271) # for reproducibility
p <- alpha + (1-alpha)*0:(N-1)/N # probabilities at which to evaluate the marginal quantile functions
X <- sapply(qF, function(qF.) qF.(p)) # input matrix for computing worst VaR
(worst.VaR <- basic_rearrange_worst_VaR(X))
stopifnot(all.equal(worst.VaR, basic_rearrange_worst_VaR(X, tol = 0)))


### 3 The more sophisticated RA() and ARA() ####################################

## Rearrangement algorithm (RA)
## In comparison to the original RA, RA()...
## - ... uses relative tolerances
## - ... deals with ties
## - ... computes more information
## - ... uses various speed-ups (see the vignette "VaR_bounds")
RA.worst <- RA(alpha, qF = qF, N = N)
str(RA.worst)
(worst.VaR.RA <- mean(RA.worst$bounds)) # the worst VaR estimate
stopifnot(all.equal(worst.VaR.RA, worst.VaR, tol = 5e-4))
(best.VaR.RA  <- mean(RA(alpha, qF = qF, N = N, method = "best.VaR")$bounds))

## Adaptive rearrangement algorithm (ARA)
## Besides the above, this...
## - ... chooses 'N' adaptively
## - ... uses the notion of a joint relative tolerance to guarantee that the two
##   bounds (on, say, worst VaR) are close
ARA.worst <- ARA(alpha, qF = qF)
ARA.best  <- ARA(alpha, qF = qF, method = "best.VaR")
str(ARA.worst) # => a lot of information
(worst.VaR.ARA <- mean(ARA.worst$bounds)) # the worst VaR estimate
(best.VaR.ARA <- mean(ARA.best$bounds)) # bounds on the best VaR
stopifnot(all.equal(worst.VaR.ARA, worst.VaR.RA, tol = 5e-4))
stopifnot(all.equal(best.VaR.ARA,  best.VaR.RA,  tol = 5e-2)) # => similar results for (A)RA()

## For more examples, see the vignette "VaR_bounds".


### 4 A (reproducing) example with Pareto margins ##############################

## Note: This can be used to reproduce the numbers of McNeil et al. (2015,
##       Table 8.1 and 8.2)

## Setup
alpha <- 0.95
d <- 8
theta <- 3
qF <- rep(list(function(p) qPar(p, theta = theta)), d)

## Worst VaR
N <- 1e5
set.seed(271)
system.time(RA.worst.VaR <- RA(alpha, qF = qF, N = N, method = "worst.VaR"))
RA.worst.VaR$bounds # The resulting bounds on worst VaR
stopifnot(RA.worst.VaR$converged, # => RA() converged
          all.equal(RA.worst.VaR$bounds[["low"]],
                    RA.worst.VaR$bounds[["up"]], tol = 1e-4)) # => bounds sufficiently close

## Best VaR
N <- 1e5
set.seed(271)
system.time(RA.best.VaR <- RA(alpha, qF = qF, N = N, method = "best.VaR"))
RA.best.VaR$bounds # The resulting bounds on best VaR
stopifnot(RA.best.VaR$converged, # => RA() converged
          all.equal(RA.best.VaR$bounds[["low"]],
                    RA.best.VaR$bounds[["up"]], tol = 1e-4)) # => bounds sufficiently close

## Best ES
N <- 2e5 # actually, we need a (much larger) N here (but that's time consuming)
set.seed(271)
system.time(RA.best.ES <- RA(alpha, qF = qF, N = N, method = "best.ES"))
RA.best.ES$bounds # The resulting bounds on best ES
stopifnot(RA.best.ES$converged, # => RA() converged
          all.equal(RA.best.ES$bounds[["low"]],
                    RA.best.ES$bounds[["up"]], tol = 5e-1)) # => bounds 'close'


### 5 How the underlying rearrange() proceeds ##################################

## See the algorithm 'at work'
A <- matrix(c(1:4, 2*(1:4)-1, 2^(0:3)), ncol = 3) # input matrix representing (dummy) quantiles
rearrange(A, tol = NULL, sample = FALSE, is.sorted = TRUE, trace = TRUE)

## An example for which the algorithm does not find the optimal rearrangement
## (but at least, with our stable sort, it terminates!)
B <- matrix(rep(1:3, 3), ncol = 3) # input matrix representing (dummy) quantiles
rearrange(B, tol = NULL, sample = FALSE, is.sorted = TRUE, trace = TRUE)
