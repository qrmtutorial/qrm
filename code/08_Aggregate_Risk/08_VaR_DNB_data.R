## By Marius Hofert and Alexander J. McNeil

## Note: The code computes the worst/best Value-at-Risk for the sum of six
##       losses for the DNB portfolio as described in


### Setup ######################################################################

library(copula)
library(qrmtools)
library(qrmdata)


### Data preparation ###########################################################

## Load the data
data(DNB) # 1% largest DNB losses

## Tail losses (above alpha quantile) for margins 1 to 3 (market, credit, asset risk)
alpha <- 0.9997 # confidence level
n <- nrow(DNB)
n.orig <- 100 * n
ii.orig <- ceiling(n.orig * alpha) : n.orig # indices of the largest 1-alpha losses in the original data
ii <- tail(1:n, n = length(ii.orig)) # number of indices => grab out those of DNB
X123 <- apply(DNB, 2, function(x) sort(x)[ii]) # losses of the first three margins
str(X123) # (751, 3)-matrix

## Specification of margins 4, 5, 6
N <- nrow(X123) # number of discretization points for each of the margins 4, 5, 6
mlog <- c(6.4741049, 6.445997, 6.0534537) # mean of underlying normal distribution (mean(log(.)))
slog <- c(0.7213475, 0.574740, 0.2489763) # sd of underlying normal distribution (sd(log(.)))


### 1 Compute worst VaR via input matrix X #####################################

### Lower bound to worst VaR

## Discretization probabilities (for the lower bound on worst VaR)
p.low <- alpha + (1-alpha)*(0:(N-1))/N

## Tail discretization of margin 4 (operational risk)
X4.low <- qlnorm(p.low, meanlog = mlog[1], sdlog = slog[1])
## Tail discretization of margin 5 (business risk)
X5.low <- qlnorm(p.low, meanlog = mlog[2], sdlog = slog[2])
## Tail discretization of margin 6 (insurance risk)
X6.low <- qlnorm(p.low, meanlog = mlog[3], sdlog = slog[3])
## The matrix X to be rearranged
X.low <- cbind(X123, X4.low, X5.low, X6.low)

## Call rearrange
set.seed(271)
rearr.low <- rearrange(X.low)
stopifnot(rearr.low$converged)
(worst.VaR.low <- rearr.low$bound)


### Upper bound to worst VaR

## Discretization probabilities (for the upper bound on worst VaR)
p.up <- alpha + (1-alpha)*(1:N)/N
p.up[length(p.up)] <- alpha + (1-alpha)*(1-1/(2*N)) # adjust upper quantiles (otherwise worst VaR is NaN)

## Tail discretizaton of margin 4 (operational risk)
X4.up <- qlnorm(p.up, meanlog = mlog[1], sdlog = slog[1])
## Tail discretizaton of margin 5 (business risk)
X5.up <- qlnorm(p.up, meanlog = mlog[2], sdlog = slog[2])
## Tail discretizaton of margin 6 (insurance risk)
X6.up <- qlnorm(p.up, meanlog = mlog[3], sdlog = slog[3])
## The matrix X to be rearranged
X.up <- cbind(X123, X4.up, X5.up, X6.up)

## Call rearrange
set.seed(271)
rearr.up <- rearrange(X.up)
stopifnot(rearr.up$converged)
(worst.VaR.up <- rearr.up$bound)


### 2 Compute worst VaR via quantile functions #################################

## Define the list of quantile functions
## Note: We artificially fill the DNB data to the original (unknown) size so
##       that quantiles are computed correctly based on the 'original' distribution
##       rather than a conditional one. How we fill the data doesn't matter as
##       these data are discarded anyways.
qF <- list(function(p) quantile(c(rep(0, 99 * n), DNB[,"Market"]), probs = p, names = FALSE, type = 1),
           function(p) quantile(c(rep(0, 99 * n), DNB[,"Credit"]), probs = p, names = FALSE, type = 1),
           function(p) quantile(c(rep(0, 99 * n), DNB[,"Asset"]),  probs = p, names = FALSE, type = 1),
           function(p) qlnorm(p, meanlog = mlog[1], sdlog = slog[1]),
           function(p) qlnorm(p, meanlog = mlog[2], sdlog = slog[2]),
           function(p) qlnorm(p, meanlog = mlog[3], sdlog = slog[3]))

## Call rearrange
set.seed(271)
RA.worst <- RA(alpha, qF = qF, N = N)
stopifnot(RA.worst$converged)
(worst.VaR <- RA.worst$bound)
stopifnot(all.equal(worst.VaR[["low"]], worst.VaR.low, tol = 1e-3))
stopifnot(all.equal(worst.VaR[["up"]],  worst.VaR.up,  tol = 1e-3))

## Plot of the 'worst VaR' copula among all losses above their alpha-quantile
pairs(pobs(RA.worst$X.rearranged$up), gap = 0, pch = ".",
      labels = as.expression(sapply(1:length(qF), function(j) bquote(L[.(j)]))))

## Note: For 'best VaR' we would need all DNB losses.