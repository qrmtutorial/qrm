## By Marius Hofert and Alexander J. McNeil

## Note: The code computes the worst/best Value-at-Risk for the sum of six
##       losses for the DNB portfolio as described in Section 3.2 of Aas, K.
##       and G. Puccetti (2014). Bounds for total economic capital: the DNB
##       case study. Extremes, 17(4), 693--715.


### Setup ######################################################################

library(copula)
library(qrmtools)

## Load the data
load("DNBdata.RData")
x <- DNBdata
str(x)
n <- nrow(x)

## Setup
alpha <- 0.9997 # confidence level


### 1 Compute worst VaR ########################################################

### 1.1 Via the input matrix X #################################################

### 1.1.1 Lower bound to worst VaR #############################################

## Tail losses (above alpha quantile) for margins 1 to 3 (market, credit,
## ownership risk)
X123 <- apply(x, 2, function(x.) sort(x.)[ceiling(n*alpha):n])
str(X123)
## Tail discretization of margin 4 (operational risk; obtained from a LN distribution)
N <- nrow(X123) # the number of discretization points in the tail
p.low <- alpha + (1-alpha)*(0:(N-1))/N # discretization probabilities for the lower bound on worst VaR
X4.low <- qlnorm(p.low, meanlog = 6.4741049, sdlog = 0.7213475)
## Tail discretization of margin 5 (business risk; obtained from a LN distribution)
X5.low <- qlnorm(p.low, meanlog = 6.445997, sdlog = 0.574740)
## Tail discretization of margin 6 (insurance risk; obtained from a LN distribution)
X6.low <- qlnorm(p.low, meanlog = 6.0534537, sdlog = 0.2489763)
## The final matrix X
X.low <- cbind(X123, X4.low, X5.low, X6.low)

## Call rearrange
set.seed(271)
rearr.low <- rearrange(X.low)
stopifnot(rearr.low$converged)
(worst.VaR.low <- rearr.low$bound)


### 1.1.2 Upper bound to worst VaR #############################################

## Tail discretizaton of margin 4 (operational risk; obtained from a LN distribution)
p.up <- alpha + (1-alpha)*(1:N)/N # discretization probabilities for the upper bound on worst VaR
p.up[length(p.up)] <- alpha + (1-alpha)*(1-1/(2*N)) # adjust upper quantile (otherwise worst VaR is NaN)
X4.up <- qlnorm(p.up, meanlog = 6.4741049, sdlog = 0.7213475)
## Tail discretizaton of margin 5 (business risk; obtained from a LN distribution)
X5.up <- qlnorm(p.up, meanlog = 6.445997, sdlog = 0.574740)
## Tail discretizaton of margin 5 (insurance risk; obtained from a LN distribution)
X6.up <- qlnorm(p.up, meanlog = 6.0534537, sdlog = 0.2489763)
## The final matrix X
X.up <- cbind(X123, X4.up, X5.up, X6.up)

## Call rearrange
set.seed(271)
rearr.up <- rearrange(X.up)
stopifnot(rearr.up$converged)
(worst.VaR.up <- rearr.up$bound)


### 1.2 Via provided quantile functions ########################################

## Define the list of quantile functions
qF <- list(function(p) quantile(x[,"market"],    probs = p, names = FALSE, type = 1),
           function(p) quantile(x[,"credit"],    probs = p, names = FALSE, type = 1),
           function(p) quantile(x[,"ownership"], probs = p, names = FALSE, type = 1),
           function(p) qlnorm(p, meanlog = 6.4741049, sdlog = 0.7213475),
           function(p) qlnorm(p, meanlog = 6.445997,  sdlog = 0.574740),
           function(p) qlnorm(p, meanlog = 6.0534537, sdlog = 0.2489763))

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


### 2 Compute best VaR #########################################################

### 2.1 Via the input matrix X #################################################

### 2.1.1 Lower bound to best VaR ##############################################

## Tail losses (below alpha quantile) for margins 1 to 3 (market, credit,
## ownership risk)
X123 <- apply(x, 2, function(x.) sort(x.)[1:ceiling(n*alpha)])
## Tail discretization of margin 4 (operational risk; obtained from a LN distribution)
N <- nrow(X123) # the number of discretization points in the tail
p.low <- alpha*(0:(N-1))/N # discretization probabilities for the lower bound on best VaR
X4.low <- qlnorm(p.low, meanlog = 6.4741049, sdlog = 0.7213475)
## Tail discretization of margin 5 (business risk; obtained from a LN distribution)
X5.low <- qlnorm(p.low, meanlog = 6.445997, sdlog = 0.574740)
## Tail discretization of margin 5 (insurance risk; obtained from a LN distribution)
X6.low <- qlnorm(p.low, meanlog = 6.0534537, sdlog = 0.2489763)
## The final matrix X
X.low <- cbind(X123, X4.low, X5.low, X6.low)

## Call rearrange
set.seed(271)
system.time(rearr.low <- rearrange(X.low, method = "best.VaR")) # takes longer
stopifnot(rearr.low$converged)
(best.VaR.low <- rearr.low$bound)


### 2.1.2 Upper bound to best VaR ##############################################

## Tail discretizaton of margin 4 (operational risk; obtained from a LN distribution)
p.up <- alpha*(1:N)/N # discretization probabilities for the upper bound on best VaR
X4.up <- qlnorm(p.up, meanlog = 6.4741049, sdlog = 0.7213475)
## Tail discretizaton of margin 5 (business risk; obtained from a LN distribution)
X5.up <- qlnorm(p.up, meanlog = 6.445997, sdlog = 0.574740)
## Tail discretizaton of margin 5 (insurance risk; obtained from a LN distribution)
X6.up <- qlnorm(p.up, meanlog = 6.0534537, sdlog = 0.2489763)
## The final matrix X
X.up <- cbind(X123, X4.up, X5.up, X6.up)

## Call rearrange
set.seed(271)
system.time(rearr.up <- rearrange(X.up, method = "best.VaR")) # takes longer
stopifnot(rearr.up$converged)
(best.VaR.up <- rearr.up$bound)


### 2.2 Via provided quantile functions ########################################

## Call rearrange
set.seed(271)
system.time(RA.best <- RA(alpha, qF = qF, N = N, method = "best.VaR")) # takes longer
stopifnot(RA.best$converged)
(best.VaR <- RA.best$bound)
stopifnot(all.equal(best.VaR[["low"]], best.VaR.low, tol = 1e-3))
stopifnot(all.equal(best.VaR[["up"]],  best.VaR.up,  tol = 1e-3))
