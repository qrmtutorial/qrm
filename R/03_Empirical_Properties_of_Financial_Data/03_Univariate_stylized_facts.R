## by Alexander McNeil
require(xts)
require(qrmdata)

data("SP500")
data("FTSE")
data("SMI")

SP500.X <- diff(log(SP500))[-1]
FTSE.X <- diff(log(FTSE))[-1]
SMI.X <- diff(log(SMI))[-1]
# remove zero returns (caused almost certainly by market closure)
SP500.X <- SP500.X[SP500.X!=0]
SMI.X <- SMI.X[SMI.X!=0]
FTSE.X <- FTSE.X[FTSE.X!=0]

INDEXES.X <- merge(SP500.X,FTSE.X,SMI.X,all=FALSE)
# compute weekly log-returns by summing within weeks
INDEXES.X.w <- apply.weekly(INDEXES.X,FUN=colSums)

plot.zoo(INDEXES.X)
pairs(as.zoo(INDEXES.X))


# Stylized Facts

# Correlation
acf(INDEXES.X)
acf(abs(INDEXES.X))
acf(INDEXES.X.w)
acf(abs(INDEXES.X.w))

# Leptokurtosis - heavy tails
par(mfrow=c(1,3))
for (i in 1:3){
  data <- INDEXES.X.w[,i]
  nm <- names(INDEXES.X.w)[i]
  qqnorm(data,main=nm)
  qqline(data,col=2)
}
par(mfrow=c(1,1))

# Clustered extreme values
tsdata <- INDEXES.X[,1]
# pick the 100 largest negative log returns
extremedata <- tsdata[rank(-as.numeric(tsdata))<=100]
plot(extremedata,type="h")
spaces <- as.numeric(diff(time(extremedata)))
qqplot(qexp(ppoints(length(spaces))), spaces)
qqline(spaces, distribution = function(p){qexp(p)},col=2)

# Poisson process for simulated iid data 
tsdata <- xts(rt(length(tsdata),df=5),time(tsdata))
extremedata <- tsdata[rank(-as.numeric(tsdata))<=100]
plot(extremedata,type="h")
spaces <- as.numeric(diff(time(extremedata)))
qqplot(qexp(ppoints(length(spaces))), spaces)
qqline(spaces, distribution = function(p){qexp(p)},col=2)

