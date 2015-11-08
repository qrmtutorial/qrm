library(lattice)

# Take a subset of Dow Jones Data
stocks <- c("AXP","EK","BA","C","KO","MSFT","HWP","INTC","JPM","DIS")
DJ.select <- DJ[,stocks]
X.daily <- window(returns(DJ.select), "1993-01-01", "2000-12-31")

by <- timeSequence(from = start(X.daily),  to = end(X.daily), by = "week")
Xdata <- aggregate(X.daily, by, sum)
dim(Xdata)
head(Xdata)

plot(Xdata)
times <- time(Xdata)
Xdata <- series(Xdata)


PCA.analysis <- princomp(Xdata)
summary(PCA.analysis)
plot(PCA.analysis)
loadings(PCA.analysis)
mu <- PCA.analysis$center
Y <- PCA.analysis$scores
PCA.factors <- Y[,1:3]
ignored.factors <- Y[,4:10]
PCA.factors <- timeSeries(PCA.factors,times)
plot(PCA.factors)
cor(PCA.factors)

Gamma <- unclass(loadings(PCA.analysis))
Gamma %*% t(Gamma)
Gamma1 <- Gamma[,1:3]
Gamma2 <- Gamma[,4:10]

tmp <- mu + Gamma1 %*% t(PCA.factors) + Gamma2 %*% t(ignored.factors)
head(t(tmp))
head(Xdata)

# Could fit a univariate GARCH model to each PCA.factor. This is called PC-GARCH model

