require(QRM)
load("INDEXES-2000-2012.RData")

# Some stock index Data 
class(INDEXES0012)
plot(INDEXES0012)
X.INDEXES <- returns(INDEXES0012)
plot(X.INDEXES)
pairs(X.INDEXES)

# Illustration of ggregation by week and month
by <- timeSequence(from = start(X.INDEXES),  to = end(X.INDEXES), by = "week")
X.INDEXES.w <- aggregate(X.INDEXES, by, sum)
head(X.INDEXES.w)
dim(X.INDEXES.w)
by <- unique(timeLastDayInMonth(time(X.INDEXES)))
X.INDEXES.m <- aggregate(X.INDEXES, by, sum)
head(X.INDEXES.m)
dim(X.INDEXES.m)

# Plotting
plot(X.INDEXES)
dim(X.INDEXES)
plot(X.INDEXES.w,type="h")
dim(X.INDEXES.w)

# Stylized Facts

# Correlation
acf(X.INDEXES)
acf(abs(X.INDEXES))
acf(X.INDEXES.w)
acf(abs(X.INDEXES.w))

# Leptokurtosis - heavy tails
par(mfrow=c(1,3))
for (i in 1:3){
  data <- X.INDEXES.w[,i]
  nm <- names(X.INDEXES.w)[i]
  qqnorm(data,main=nm)
  qqline(data,col=2)
}
par(mfrow=c(1,1))

# Clustered extreme values
tsdata <- X.INDEXES[,1]
u <- findthreshold(tsdata,ne=100)
extremedata <- tsdata[tsdata>u,]
plot(extremedata,type="h")
spaces <- as.numeric(diff(time(extremedata)))
qqplot(qexp(ppoints(length(spaces))), spaces)
qqline(spaces, distribution = function(p){qexp(p)},col=2)

# Poisson process for simulated strict white noise data 
tsdata <- timeSeries(rt(length(tsdata),df=5),time(tsdata))
u <- findthreshold(tsdata,ne=100)
extremedata <- tsdata[tsdata>u,]
plot(extremedata,type="h")
spaces <- as.numeric(diff(time(extremedata)))
qqplot(qexp(ppoints(length(spaces))), spaces)
qqline(spaces, distribution = function(p){qexp(p)},col=2)

