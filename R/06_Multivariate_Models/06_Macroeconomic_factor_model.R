library(lattice)
library(xts)
library(qrmdata)

data("DJ_const")
data("DJ")
# Take a subset of Dow Jones Data
# We will take stocks with complete record of returns since 1990
tickers <- names(DJ_const)[-c(5,10,25,26,27)]
tickers
DJ_const.select <- DJ_const['1990-01-01/',tickers]
DJ.select <- DJ['1990-01-01/']

# Compute log returns
DJ_const.r <- diff(log(DJ_const.select))[-1,]
DJ.r <- diff(log(DJ.select))[-1]

# Compute monthly log returns
Xdata <- apply.monthly(DJ_const.r,FUN=colSums)
Fdata <- apply.monthly(DJ.r,FUN=colSums)
dim(Xdata)
length(Fdata)
head(Xdata)
head(Fdata)

plot.zoo(Xdata,type="h")
plot(Fdata,type="h")


# Fit Multivariate Regression
out <- lm(Xdata ~ Fdata)
names(out)
summary(out)
par.ests <- coef(out)
par.ests

# Get parameter estimates
a <- par.ests[1,]
B <- par.ests[-1,]

# Get residuals and estimate covariance matrix
epsilon <- resid(out)
# Are off-diagonal correlations small ?
round(cor(epsilon),2)
# Are the errors uncorrelated with factors
round(cor(epsilon,Fdata),2)

# Colour plot of correlation matrix
levelplot(cor(epsilon))
col.l <- colorRampPalette(c('green', 'orange', 'red','purple'))(256)
levelplot(cor(epsilon),col.regions=col.l)


# Construct implied covariance and correlation matrix
Psi <- var(epsilon)
Omega <- as.matrix(var(Fdata))
Sigma <- B %*% Omega %*% t(B) + diag(diag(Psi))
dimnames(Sigma) <- list(dimnames(par.ests)[[2]],dimnames(par.ests)[[2]])
R.fact <- cov2cor(Sigma)

# look at discrepancies between factor model correlation matrix
# and sample correlation matrix
levelplot(R.fact-cor(Xdata))


