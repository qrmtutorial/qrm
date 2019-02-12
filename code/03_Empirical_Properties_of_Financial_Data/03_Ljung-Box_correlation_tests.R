## By Alexander McNeil and Marius Hofert


### Setup ######################################################################

library(xts)
library(qrmdata)
library(qrmtools)


### 1 Constituent data #########################################################

## Dow Jones constituent data
data("DJ_const")

## We extract a time period and remove the time series 'Visa' (as it only has
## a very short history in the index)
DJdata <- DJ_const['2006-12-29/2015-12-31',-which(names(DJ_const) == "V")]

## Use plot for zoo objects to get multiple plots
plot.zoo(DJdata, xlab = "Time", main = "DJ component series (without Visa)")

## Build log-returns and aggregate to obtain monthly log-returns
X <- returns(DJdata) # could also work with negative log-returns
X.m <- apply.monthly(X, FUN = colSums)


### 2 Ljung--Box tests of serial independence of stationary data ###############

## Compute (lists of) Ljung--Box tests
LB.raw   <- apply(X,        2, Box.test, lag = 10, type = "Ljung-Box")
LB.abs   <- apply(abs(X),   2, Box.test, lag = 10, type = "Ljung-Box") # could also work with squared log-returns
LB.raw.m <- apply(X.m,      2, Box.test, lag = 10, type = "Ljung-Box")
LB.abs.m <- apply(abs(X.m), 2, Box.test, lag = 10, type = "Ljung-Box")

## Extract p-values
p.LB.raw   <- sapply(LB.raw,   `[[`, "p.value")
p.LB.abs   <- sapply(LB.abs,   `[[`, "p.value")
p.LB.raw.m <- sapply(LB.raw.m, `[[`, "p.value")
p.LB.abs.m <- sapply(LB.abs.m, `[[`, "p.value")
round(cbind(p.LB.raw, p.LB.abs, p.LB.raw.m, p.LB.abs.m), 2)


### 3 Reproduce Ljung--Box tests from the QRM book #############################

## Up to minor differences (see below), this reproduces Table 3.1 in the
## QRM book (2015) (see also Table 4.1 in the QRM book (2005)):

## Note: This uses older DJ data from 'QRM'
DJ.QRM <- QRM::DJ
DJ.old <- as.xts(DJ.QRM)
DJ.old <- DJ.old['1993-01-01/2000-12-31']
X.old <- returns(DJ.old)
X.old.m <- apply.monthly(X.old, FUN = colSums)

## Compute (lists of) Ljung--Box tests
LB.raw   <- apply(X.old,        2, Box.test, lag = 10, type = "Ljung-Box")
LB.abs   <- apply(abs(X.old),   2, Box.test, lag = 10, type = "Ljung-Box")
LB.raw.m <- apply(X.old.m,      2, Box.test, lag = 10, type = "Ljung-Box")
LB.abs.m <- apply(abs(X.old.m), 2, Box.test, lag = 10, type = "Ljung-Box")

## Extract p-values
p.LB.raw   <- sapply(LB.raw,   `[[`, "p.value")
p.LB.abs   <- sapply(LB.abs,   `[[`, "p.value")
p.LB.raw.m <- sapply(LB.raw.m, `[[`, "p.value")
p.LB.abs.m <- sapply(LB.abs.m, `[[`, "p.value")
(res <- round(cbind(p.LB.raw, p.LB.abs, p.LB.raw.m, p.LB.abs.m), 2))

## Note: The minor differences to the tables in the book come from a different
##       approach. The tables in the book were produces without 'xts' objects:
library(timeSeries) # Caution: returns() now from timeSeries
X.old. <- returns(DJ.QRM)
X.old. <- window(X.old., start = timeDate("1993-01-01"), end = timeDate("2000-12-31"))
X.old.m. <- aggregate(X.old., by = unique(timeLastDayInMonth(time(X.old.))), sum)
## Compute (lists of) Ljung--Box tests
LB.raw.   <- apply(X.old.,        2, Box.test, lag = 10, type = "Ljung-Box")
LB.abs.   <- apply(abs(X.old.),   2, Box.test, lag = 10, type = "Ljung-Box")
LB.raw.m. <- apply(X.old.m.,      2, Box.test, lag = 10, type = "Ljung-Box")
LB.abs.m. <- apply(abs(X.old.m.), 2, Box.test, lag = 10, type = "Ljung-Box")
## Extract p-values
p.LB.raw.   <- sapply(LB.raw.,   `[[`, "p.value")
p.LB.abs.   <- sapply(LB.abs.,   `[[`, "p.value")
p.LB.raw.m. <- sapply(LB.raw.m., `[[`, "p.value")
p.LB.abs.m. <- sapply(LB.abs.m., `[[`, "p.value")
## Result
(res. <- round(cbind(p.LB.raw., p.LB.abs., p.LB.raw.m., p.LB.abs.m.), 2))
## Differences
summary(res-res.)

