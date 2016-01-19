library(lattice)
library(xts)
library(qrmdata)

data("DJ_const")
# Take a subset of Dow Jones Data
# We will take stocks with complete record of returns since 1990
tickers <- names(DJ_const)[-c(5,10,25,26,27)]
tickers
DJ_const.select <- DJ_const['1990-01-01/',tickers]

# Compute log returns
DJ_const.r <- diff(log(DJ_const.select))[-1,]

# Compute monthly log returns
Xdata <- apply.monthly(DJ_const.r,FUN=colSums)
dim(Xdata)
d <- dim(Xdata)[2]
head(Xdata)
plot.zoo(Xdata,type="h")


PCA.analysis <- princomp(Xdata)
summary(PCA.analysis)
plot(PCA.analysis)
loadings(PCA.analysis)
mu <- PCA.analysis$center
Y <- PCA.analysis$scores
PCA.factors <- Y[,1:3]
ignored.factors <- Y[,4:d]
PCA.factors <- xts(PCA.factors,time(Xdata))
plot.zoo(PCA.factors,type="h")
cor(PCA.factors)

Gamma <- unclass(loadings(PCA.analysis))
levelplot(Gamma %*% t(Gamma))
Gamma1 <- Gamma[,1:3]
Gamma2 <- Gamma[,4:d]

tmp <- mu + Gamma1 %*% t(PCA.factors) + Gamma2 %*% t(ignored.factors)
head(t(tmp))
head(Xdata)

# Could fit a univariate GARCH model to each PCA.factor. This is called PC-GARCH model

