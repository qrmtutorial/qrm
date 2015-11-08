library(lattice)

# Take a subset of Dow Jones Data
stocks <- c("AXP","EK","BA","C","KO","MSFT","HWP","INTC","JPM","DIS")
DJ.select <- DJ[,stocks]
X.daily <- window(returns(DJ.select), "1993-01-01", "2000-12-31")
F.daily <- window(returns(dji), "1993-01-01", "2000-12-31")

by <- timeSequence(from = start(X.daily),  to = end(X.daily), by = "week")
Xdata <- aggregate(X.daily, by, sum)
Fdata <- aggregate(F.daily, by, sum)
dim(Xdata)
length(Fdata)
head(Xdata)
head(Fdata)

plot(Xdata)
plot(Fdata)
times <- time(Xdata)
Xdata <- series(Xdata)
Fdata <- series(Fdata)

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
cor.data <- cor(Xdata)
round(R.fact,2)
round(cor.data,2)



# To get the book example rerun analyses with these data
# stocks <- c("MO","KO","EK","HWP","INTC","MSFT","IBM","MCD","WMT","DIS")
# DJ.select <- DJ[,stocks]
# Xdata <- window(returns(DJ.select), "1992-01-01", "1998-12-31")
# Fdata <- window(returns(dji), "1992-01-01", "1998-12-31")
