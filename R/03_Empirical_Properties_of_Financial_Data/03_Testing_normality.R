# REQUIRED LIBRARIES
require(xts)
require(qrmdata)

# Prepare risk factor data
data("DJ_const")
DJ.X.all <- diff(log(DJ_const))[-1,]
DJ.X <- DJ.X.all['1993-01-01/2000-12-31',c("KO","MSFT","INTC","DIS")]
DJ.X.w <- apply.weekly(DJ.X, FUN=colSums)
DJ.X.m <- apply.monthly(DJ.X, FUN=colSums)
DJ.X.q <- apply.quarterly(DJ.X, FUN=colSums)

# Are Asset Returns Normally Distributed ?
apply(DJ.X,2,shapiro.test)
qqnorm(DJ.X[,1])
qqline(DJ.X[,1],col=2)
apply(DJ.X.w,2,shapiro.test)
qqnorm(DJ.X.w[,1])
qqline(DJ.X.w[,1],col=2)
apply(DJ.X.m,2,shapiro.test)
qqnorm(DJ.X.m[,4])
qqline(DJ.X.m[,4],col=2)
apply(DJ.X.q,2,shapiro.test)
qqnorm(DJ.X.q[,1])
qqline(DJ.X.q[,1],col=2)

# EXERCISE. Try the Jarque Bera test as an alternative. It can be found in "tseries" package.

# How real normal data should behave (repeat many times)
data.normal <- rnorm(1000)
hist(data.normal)
qqnorm(data.normal)
qqline(data.normal,col=2)
shapiro.test(data.normal)


