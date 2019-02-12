## By Alexander McNeil and Marius Hofert

## Univariate stylized facts:
## (U1) Return series are not iid although they show little serial correlation;
## (U2) Series of absolute (or squared) returns show profound serial correlation;
## (U3) Conditional expected returns are close to zero;
## (U4) Volatility (conditional standard deviation) appears to vary over time;
## (U5) Extreme returns appear in clusters;
## (U6) Return series are leptokurtic or heavy-tailed (power-like tail).


### Setup ######################################################################

library(xts)
library(qrmdata)
library(qrmtools)


### 1 Index data ###############################################################

## S&P 500, FTSE, SMI
data("SP500")
data("FTSE")
data("SMI")

## Build log-returns
SP500.X <- returns(SP500)
FTSE.X <- returns(FTSE)
SMI.X <- returns(SMI)
## Remove zero returns (caused almost certainly by market closure)
SP500.X <- SP500.X[SP500.X != 0]
SMI.X <- SMI.X[SMI.X != 0]
FTSE.X <- FTSE.X[FTSE.X != 0]

## Merge and aggregate
X <- merge(SP500 = SP500.X, FTSE = FTSE.X, SMI = SMI.X, all = FALSE)
X.w <- apply.weekly(X, FUN = colSums) # weekly log-returns (summing within weeks)


### 2 (U1) and (U2) ############################################################

## Auto- and crosscorrelations of returns and absolute returns
acf(X)
acf(abs(X))
acf(X.w)
acf(abs(X.w))


### 3 (U3) and (U4) ############################################################

## Plot of returns
plot.zoo(X, xlab = "Time", main = "Log-returns") # this (partly) confirms (U3)

## Plot of volatilities
X.vols <- apply.monthly(X, FUN = function(x) apply(x, 2, sd)) # componentwise volatilities
plot.zoo(X.vols, xlab = "Time", main = "Volatility estimates")


### 3 (U5) #####################################################################

## Plot the largest losses (already partly confirms (U5))
L <- -SP500.X # consider -log-returns of S&P 500 here
r <- 0.01 # probability to determine large losses
u <- quantile(L, probs = 1 - r) # determine empirical (1-r)-quantile
xtr.L <- L[L > u] # exceedances over the chosen empirical quantile = largest so-many losses
plot(as.numeric(xtr.L), type = "h", xlab = "Time",
     ylab = substitute("Largest"~r.*"% of losses ("*n.~"losses)",
                       list(r. = 100 * r, n. = length(xtr.L))))

## Now consider spacings (certainly not exponentially distributed)
spcs <- as.numeric(diff(time(xtr.L)))
qq_plot(spcs, FUN = function(p) qexp(p, rate = r)) # r = exceedance probability

## Compare with a Poisson process (simulated iid data)
set.seed(271)
L. <- rt(length(L), df = 3) # simulate iid data from a t_3 distribution
u. <- quantile(L., probs = 1 - r) # determine empirical (1-r)-quantile
xtr.L. <- L.[L. > u.] # exceedances
plot(xtr.L., type = "h", xlab = "Time",
     ylab = substitute("Largest"~r.*"% of losses ("*n.~"losses)",
                       list(r. = 100 * r, n. = length(xtr.L.))))

## Sapcings
spcs. <- diff(which(L. > u.))
qq_plot(spcs., FUN = function(p) qexp(p, rate = r))


### 4 (U6) #####################################################################

## Show leptokurtosis (heavier tails than a normal) by Q-Q plots against a normal
layout(t(1:3)) # (1,3)-layout
for (i in 1:3)
    qq_plot(X.w[,i], FUN = qnorm, method = "empirical", # ... to account for the unknown location (mu) and scale (sig) of qnorm()
            main = names(X.w)[i]) # equivalent to qqnorm(); qqline()
layout(1) # restore layout
## One can also look at various formal tests
