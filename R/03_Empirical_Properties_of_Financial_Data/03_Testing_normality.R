# REQUIRED LIBRARIES
require(QRM)

# DATA
load("INDEXES-2000-2012.RData")


# Prepare risk factor data
stocks <- c("KO","MSFT","INTC","DIS")
DJ.select <- DJ[,stocks]
X.DJ <- window(returns(DJ.select), "1993-01-01", "2000-12-31")
by <- timeSequence(from = start(X.DJ),  to = end(X.DJ), by = "week")
X.DJ.w <- aggregate(X.DJ, by, sum)
by <- unique(timeLastDayInMonth(time(X.DJ)))
X.DJ.m <- aggregate(X.DJ, by, sum)
by <- unique(timeLastDayInQuarter(time(X.DJ)))
X.DJ.q <- aggregate(X.DJ, by, sum)

# Are Asset Returns Normally Distributed ?
apply(X.DJ,2,shapiro.test)
qqnorm(X.DJ[,1])
qqline(X.DJ[,1],col=2)
apply(X.DJ.w,2,shapiro.test)
qqnorm(X.DJ.w[,1])
qqline(X.DJ.w[,1],col=2)
apply(X.DJ.m,2,shapiro.test)
qqnorm(X.DJ.m[,4])
qqline(X.DJ.m[,4],col=2)
apply(X.DJ.q,2,shapiro.test)
qqnorm(X.DJ.q[,1])
qqline(X.DJ.q[,1],col=2)

# EXERCISE. Try the Jarque Bera test as analternative. It can be found in "tseries" package.

# EXERCISE. Test normality of components of X.INDEXES

X.INDEXES <- returns(INDEXES0012)
by <- timeSequence(from = start(X.INDEXES),  to = end(X.INDEXES), by = "week")
X.INDEXES.w <- aggregate(X.INDEXES, by, sum)
by <- unique(timeLastDayInMonth(time(X.INDEXES)))
X.INDEXES.m <- aggregate(X.INDEXES, by, sum)


# How real normal data should behave (repeat many times)
data.normal <- rnorm(1000)
hist(data.normal)
qqnorm(data.normal)
qqline(data.normal,col=2)
shapiro.test(data.normal)


