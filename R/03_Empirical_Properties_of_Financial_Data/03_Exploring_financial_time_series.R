
# REQUIRED LIBRARIES
require(QRM)

# DATA
load("INDEXES-2000-2012.RData")
load("Canadian-ZCB-Yields.RData")

# Risk Factor Data (Daily)

# Stock Prices (in QRM)
data(DJ)
class(DJ)
dim(DJ)
stocks <- c("KO","MSFT","INTC","DIS")
DJ.select <- DJ[,stocks]
plot(DJ.select)
?returns
X.DJ <- window(returns(DJ.select), "1993-01-01", "2000-12-31")
plot(X.DJ)
summary(X.DJ)
pairs(X.DJ)

# Aggregating by week, month, quarter
by <- timeSequence(from = start(X.DJ),  to = end(X.DJ), by = "week")
X.DJ.w <- aggregate(X.DJ, by, sum)
head(X.DJ.w)
dim(X.DJ.w)
by <- unique(timeLastDayInMonth(time(X.DJ)))
X.DJ.m <- aggregate(X.DJ, by, sum)
head(X.DJ.m)
dim(X.DJ.m)
by <- unique(timeLastDayInQuarter(time(X.DJ)))
X.DJ.q <- aggregate(X.DJ, by, sum)
head(X.DJ.q)
dim(X.DJ.q)


# Index Data (new data)

class(INDEXES0012)
plot(INDEXES0012)
X.INDEXES <- returns(INDEXES0012)
plot(X.INDEXES)
pairs(X.INDEXES)

# Aggregating by week, month
by <- timeSequence(from = start(X.INDEXES),  to = end(X.INDEXES), by = "week")
X.INDEXES.w <- aggregate(X.INDEXES, by, sum)
head(X.INDEXES.w)
dim(X.INDEXES.w)
by <- unique(timeLastDayInMonth(time(X.INDEXES)))
X.INDEXES.m <- aggregate(X.INDEXES, by, sum)
head(X.INDEXES.m)
dim(X.INDEXES.m)


# Exchange Rates (in QRM)
data(FXGBP)
class(FXGBP)
head(FXGBP)
tail(FXGBP)
plot(FXGBP)
FX.select <- window(FXGBP, "1996-01-01", "2003-12-31")
plot(FX.select)
X.FX <- returns(FX.select)
plot(X.FX)

# Aggregating by month
by <- unique(timeLastDayInMonth(time(X.FX)))
X.FX.m <- aggregate(X.FX, by, sum)
head(X.FX.m)
dim(X.FX.m)


# Zero-Coupon-Bond Yields (new data)
class(ZCB)
dim(ZCB)
head(ZCB)
ZCB <- window(100*ZCB, "2002-01-02", "2011-12-30")
plot(ZCB[,1:10])
X.ZCB <- diff(ZCB)
X.ZCB3 <- X.ZCB[,c(1,8,40)]
plot(X.ZCB3)
pairs(X.ZCB3)

