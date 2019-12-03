## By Alexander McNeil and Marius Hofert

## Loading, exploring, merging, aggregating and presenting financial time series


### Setup ######################################################################

library(xts)
library(qrmdata)
library(qrmtools)


### Stock prices ###############################################################

## Dow Jones stock price data
data("DJ_const")
str(DJ_const)

## We extract a time period and take the first 10 stocks
DJdata <- DJ_const['2006-12-29/2015-12-31', 1:10]

## Use plot for zoo objects to get multiple plots
plot.zoo(DJdata, xlab = "Time", main = "DJ (10 component series)")
X <- returns(DJdata) # or diff(log(DJdata))[-1,]
head(X)
plot.zoo(X, xlab = "Time", main = "Log-returns of 10 DJ component series")

## Aggregating log returns by summation for each column
## Weekly
X.w <- apply.weekly(X, FUN = colSums)
dim(X.w)
plot.zoo(X.w, type = "h", xlab = "Time", main = "Weekly log-returns") # 'h' = histogram = vertical lines only
## Monthly
X.m <- apply.monthly(X, FUN = colSums)
dim(X.m)
plot.zoo(X.m, type = "h", xlab = "Time", main = "Monthly log-returns")
## Quarterly
X.q <- apply.quarterly(X, FUN = colSums)
dim(X.q)
plot.zoo(X.q, type = "h", xlab = "Time", main = "Quarterly log-returns")


### Stock indexes ##############################################################

## Load stock index data
data("SP500") # S&P 500
data("FTSE") # FTSE
data("SMI") # SMI
plot.zoo(SP500, xlab = "Time", ylab = "S&P 500")

## Merge the three time series
all <- merge(SP500, FTSE, SMI)
nms <- c("S&P 500", "FTSE", "SMI")
colnames(all) <- nms
plot.zoo(all, xlab = "Time", main = "All")

## Merge and retain only days where all three time series are available
all.avail <- merge(SP500, FTSE, SMI, all = FALSE)
colnames(all.avail) <- nms
plot.zoo(all.avail, xlab = "Time", main = "All available")

## Compute returns
SP500.X <- returns(SP500)
FTSE.X  <- returns(FTSE)
SMI.X   <- returns(SMI)
X <- merge(SP500.X, FTSE.X, SMI.X, all = FALSE)
colnames(X) <- nms
plot.zoo(X, xlab = "Time", main = "Log-returns")
pairs(as.zoo(X), gap = 0, cex = 0.4)

## Aggregating
## By week
X.w <- apply.weekly(X, FUN = colSums)
plot.zoo(X.w, xlab = "Time", main = "Weekly log-returns", type = "h")
## By month
X.m <- apply.monthly(X, FUN = colSums)
plot.zoo(X.m, xlab = "Time", main = "Monthly log-returns", type = "h")


### Exchange rates #############################################################

## Load exchange rate data
data("GBP_USD") # 1 GBP in USD
data("EUR_USD") # 1 EUR in USD
data("JPY_USD") # 1 JPY in USD
data("CHF_USD") # 1 CHF in USD
FX <- merge(GBP_USD, EUR_USD, JPY_USD, CHF_USD)
head(FX)
plot.zoo(FX, xlab = "Time", main = "Exchange rates to USD")
X <- returns(FX)
plot.zoo(X, xlab = "Time", main = "Log-returns of exchange rates to USD",
         col = c("black", "royalblue3", "maroon3", "darkorange2"))


### Zero-coupon bond yields ####################################################

## Load zero-coupon bond yield data (in USD)
## Note: Typically, yield = ((face value / current bond price)^(1/maturity) - 1) * 100%
##       (as face value = current price * (1 + yield)^maturity
data("ZCB_USD")
dim(ZCB_USD) # => 30-dimensional; each dimension is a maturity
head(ZCB_USD)
ZCB <- ZCB_USD['2002-01-02/2011-12-30']
plot.zoo(ZCB, xlab = "Time", main = "Percentage yields")

## Compute differences (first row is removed)
X <- returns(ZCB, method = "diff")

## Pick 3 maturities
X3 <- X[, c(1, 5, 30)]
plot.zoo(X3, xlab = "Time", main = "Differences (3 maturities)")
pairs(as.zoo(X3), gap = 0, cex = 0.4)

## Plot the corresponding "pseudo-observations" (componentwise scaled ranks)
U3 <- apply(X3, 2, rank) / (ncol(X3) + 1)
pairs(as.zoo(U3), gap = 0, cex = 0.4)