# Library required
require(QRM)

# The bivariate normal distribution
BiDensPlot(func = dmnorm, xpt=c(-4,4), ypts=c(-4,4),mu = c(0, 0), Sigma = equicorr(2, -0.7))

# Some trivariate simulated data
S <- equicorr(d = 3, rho = 0.7)
data <- rmnorm(1000, Sigma = S)

# Re-estimate the parameters
fit.norm(data)

# Some 10-d simulated data
S <- equicorr(d = 10, rho = 0.6)
data <- rmnorm(1000, Sigma = S)

# Test of marginal normality
apply(data,2,shapiro.test)

# Tests of joint normality
MardiaTest(data)
jointnormalTest(data)

# Dow Jones Data
data(DJ)
r <- returns(DJ)
stocks <- c("AXP","EK","BA","C","KO","MSFT","HWP","INTC","JPM","DIS")
ss <- window(r[, stocks], "1993-01-01", "2000-12-31")
plot(ss)

# Test of marginal normality
apply(ss,2,shapiro.test)

# Test for joint normality
MardiaTest(ss)
jointnormalTest(ss)

# Talking point: data which are only marginally normal
data = rAC("gumbel", n=1000, d=3, theta=3)
data = apply(data,2,qnorm)
pairs(data)
apply(data,2,shapiro.test)
MardiaTest(ss)
jointnormalTest(ss)
