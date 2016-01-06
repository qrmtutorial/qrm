# by Alexander McNeil & Marius Hofert
# REQUIRED LIBRARIES
require(xts)
require(qrmdata)

# STOCK PRICES
# Dow Jones stock price data
data("DJ_const")
class(DJ_const)
# We extract a time period and take the first 10 stocks
DJdata <- DJ_const['2006-12-29/2015-12-31',1:10]
# Use plot for zoo objects to get multiple plots
plot.zoo(DJdata)
DJ.X <- diff(log(DJdata))[-1,]
head(DJ.X)
plot.zoo(DJ.X)

# aggregating log returns by summation
DJ.X.w <- apply.weekly(DJ.X, FUN=colSums)
dim(DJ.X.w)
plot.zoo(DJ.X.w,type="h")
DJ.X.m <- apply.monthly(DJ.X, FUN=colSums)
dim(DJ.X.m)
plot.zoo(DJ.X.m,type="h")
DJ.X.q <- apply.quarterly(DJ.X, FUN=colSums)
dim(DJ.X.q)
plot.zoo(DJ.X.q,type="h")

# STOCK INDEXES
data("SP500")
data("FTSE")
data("SMI")
class(SP500)
plot(SP500)
# merge all the data
INDEXESall <- merge(SP500,FTSE,SMI,all=TRUE)
plot.zoo(INDEXESall)
# merge and retain only days where all indexes have values
INDEXES <- merge(SP500,FTSE,SMI,all=FALSE)
plot.zoo(INDEXES)

# compute returns
SP500.X <- diff(log(SP500))[-1]
FTSE.X <- diff(log(FTSE))[-1]
SMI.X <- diff(log(SMI))[-1]
INDEXES.X <- merge(SP500.X,FTSE.X,SMI.X,all=FALSE)
plot.zoo(INDEXES.X)
pairs(as.zoo(INDEXES.X))

# Aggregating by week, month
INDEXES.X.w <- apply.weekly(INDEXES.X,FUN=colSums)
plot.zoo(INDEXES.X.w,type="h")
INDEXES.X.m <- apply.monthly(INDEXES.X, FUN=colSums)
plot.zoo(INDEXES.X.m,type="h")


# EXCHANGE RATES 
data("GBP_USD")
data("EUR_USD")
data("JPY_USD")
data("CHF_USD")
FX <- merge(GBP_USD,EUR_USD,JPY_USD,CHF_USD,all=TRUE)
head(FX)
plot.zoo(FX)
FX.X <- diff(log(FX))[-1,]
plot.zoo(FX.X,col=1:4)



# ZERO-COUPON BOND YIELDS
data("ZCB_CA")
dim(ZCB_CA)
head(ZCB_CA)
# change to percentages and select period
ZCB <- 100*ZCB_CA['2002-01-02/2011-12-30']
plot.zoo(ZCB[,1:10],col=1:10)
# compute simple returns and remove first row
ZCB.X <- diff(ZCB)[-1,]
# pick 3 maturities
ZCB.X.3 <- ZCB.X[,c(1,8,40)]
plot.zoo(ZCB.X.3,col=1:3)
pairs(as.zoo(ZCB.X.3))

