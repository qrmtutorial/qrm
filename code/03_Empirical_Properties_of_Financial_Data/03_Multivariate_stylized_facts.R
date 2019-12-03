## By Alexander McNeil and Marius Hofert

## Multivariate stylized facts:
## (M1) Return series show little cross-correlation, except for contemporaneous
##      returns;
## (M2) Absolute (or squared) returns show profound cross-correlation;
## (M3) Correlations between contemporaneous returns vary over time;
## (M4) Extreme returns often coincide with extreme returns in (several) other
##      time series.


### Setup ######################################################################

library(xts)
library(nvmix) # for rNorm(), rStudent()
library(qrmdata)
library(qrmtools)


### 1 Index data ###############################################################

## FTSE and SMI index data
data("FTSE")
data("SMI")

## Build log-returns and merge
FTSE.X <- returns(FTSE)
SMI.X <- returns(SMI)
X <- merge(FTSE = FTSE.X, SMI = SMI.X, all = FALSE)


### 2 (M1) and (M2) ############################################################

## Autocorrelations of returns and absolute returns
acf(X) # raw values
acf(abs(X)) # absolute values


### 3 (M3) #####################################################################

## Plot cross-correlation (and volatilities)
X.cor  <- apply.monthly(X, FUN = cor)[,2] # cross-correlations estimated based on monthly non-overlapping windows
X.vols <- apply.monthly(X, FUN = function(x) apply(x, 2, sd)) # componentwise volatilities
X.cor.vols <- merge(X.cor, X.vols)
names(X.cor.vols) <- c("Cross-correlation", "Volatility FTSE", "Volatility SMI")
plot.zoo(X.cor.vols, xlab = "Time", main = "Cross-correlation and volatility estimates")
## => Cross-correlation varies over time; below we further investigate the
##    relationship between volatility and cross-correlation

## Extract volatilities and compute Fisher-transformed cross-correlations
FTSE.sig <- as.numeric(X.vols[,"FTSE"])
SMI.sig  <- as.numeric(X.vols[,"SMI"])
fisher <- function(r) log((1 + r)/(1 - r))/2 # Fisher (z-)transformation
rho <- fisher(X.cor) # Fisher transform cross correlations

## Plot Fisher-transformed cross-correlation against volatility with regression lines
## Note: If (X,Y) are jointly normal with correlation rho, fisher(<sample correlation>)
##       is approximately N(log((1+rho)/(1-rho))/2, 1/(n-3)) distributed (n = sample size).

## FTSE
plot(FTSE.sig, rho, xlab = "Estimated volatility", ylab = "Estimated cross-correlation")
reg <- lm(rho ~ FTSE.sig)
summary(reg)
abline(reg)

## SMI
plot(SMI.sig, rho, xlab = "Estimated volatility", ylab = "Estimated cross-correlation")
reg <- lm(rho ~ SMI.sig)
summary(reg)
abline(reg)

## EXERCISE: Now try and compare with these SWN (strict white noise) data
set.seed(271)
X.t <- xts(rStudent(n = nrow(X), df = 3, scale = cor(X)), time(X))
X.N <- xts(rNorm(n = nrow(X), scale = cor(X)), time(X))


### 4 (M4) #####################################################################

plot.zoo(X, xlab = "Time", main = "Log-returns") # seems like extremes occur together

X. <- apply(X, 2, rank) / (nrow(X) + 1)
plot(X., main = "Componentwise scaled ranks") # now better visible (bottom-left, top-right corners)
## Note: More on that ("pseudo-observations") in Chapter 7
