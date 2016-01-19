## by Alexander McNeil

require(QRM)

plot(danish)
MEplot(danish)
abline(v=10)
abline(v=20)

# Digression: a little comparison of mean excess plots
data1 <- abs(rt(5000,df=4))
MEplot(data1)
data2 <- rexp(5000)
data3 <- abs(rnorm(5000))
MEplot(data2)
MEplot(data3)
data4 <- runif(5000)
MEplot(data4)

# Fit GPD model above 10
?fit.GPD
mod1 <- fit.GPD(danish,threshold=10)
mod1
plotFittedGPDvsEmpiricalExcesses(danish, threshold = 10)

RMs1 <- RiskMeasures(mod1, c(0.99,0.995))
RMs1
plotTail(mod1)
showRM(mod1, alpha = 0.99, RM = "VaR", method = "BFGS")
abline(h=0.01)
showRM(mod1, alpha = 0.99, RM = "ES", method = "BFGS")

# Effect of changing threshold on xi
xiplot(danish)

# Fit GPD model above 20
mod2 <- fit.GPD(danish,threshold=20)
mod2
plotFittedGPDvsEmpiricalExcesses(danish, threshold = 20)

RMs2 <- RiskMeasures(mod2, c(0.99,0.995))
RMs2
plotTail(mod2)
showRM(mod2, alpha = 0.99, RM = "VaR", method = "BFGS")


