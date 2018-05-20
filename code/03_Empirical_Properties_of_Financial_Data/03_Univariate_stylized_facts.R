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

## Plot the 100 largest losses (already partly confirms (U5))
L <- -SP500.X # consider -log-returns of S&P 500 here
xtr.L <- L[rank(as.numeric(-L)) <= 100] # 100 smallest gains = largest losses
plot(as.numeric(xtr.L), type = "h", xlab = "Time", ylab = "100 largest losses")

## Now consider spacings (certainly not exponentially distributed)
spcs <- as.numeric(diff(time(xtr.L)))
qq_plot(spcs, FUN = qexp, method = "empirical")
if(FALSE) {
    ## This is equivalent to:
    qqplot(qexp(ppoints(length(spcs))), y = spcs, xlab = "Theoretical quantiles",
           ylab = "Sample quantiles")
    qqline(spcs, distribution = qexp)
}

## Compare with a Poisson process (simulated iid data)
L. <- xts(rt(length(L), df = 5), time(L)) # simulate from a t_5
xtr.L. <- L.[rank(as.numeric(-L.)) <= 100]
plot(as.numeric(xtr.L.), type = "h", xlab = "Time", ylab = "100 largest losses")
spcs. <- as.numeric(diff(time(xtr.L.)))
qq_plot(spcs., FUN = qexp, method = "empirical")


### 4 (U6) #####################################################################

## Show leptokurtosis (heavier tails than a normal) by Q-Q plots against a normal
layout(t(1:3)) # (1,3)-layout
for (i in 1:3)
    qq_plot(X.w[,i], FUN = qnorm, method = "empirical", main = names(X.w)[i])
layout(1) # restore layout
