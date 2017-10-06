## By Alexander McNeil

library(xts)
library(qrmdata)
library(QRM)

data("DJ_const")
# create log returns
DJ.r <- diff(log(DJ_const))[-1,]
# aggregate by week
DJ.w <- apply.weekly(DJ.r, FUN = colSums)
# convert to weekly percentage losses
DJ.wp <- -(exp(DJ.w)-1)*100
# take data from start of millenium
DJ.wp <- DJ.wp['2000-01-01/']
plot.zoo(DJ.wp[,1:12], type = "h")

# Show mean excess plots
op <- par(mfrow = c(3,4), mar = c(3,3,2,1))
for (i in 1:12){
  data <- DJ.wp[,i]
  MEplot(data[data>0], main = names(DJ_const)[i],
         xlab = "", ylab = "")
}
par(op)

# Select Apple share price
data <- DJ.wp[,"AAPL"]
u <- 6

# Mean excess plot of losses omitting gains
MEplot(data[data>0])
abline(v = u)

# Fit GPD model above u
mod1 <- fit.GPD(data,u)
mod1
plotFittedGPDvsEmpiricalExcesses(data, threshold = u)

# Compute VaR and expected shortfall
RMs1 <- RiskMeasures(mod1, c(0.99,0.995))
RMs1
plotTail(mod1)
# Illustrate 99% quantile calculation
abline(h = 0.01)
abline(v = RMs1[1,"quantile"])

# Illustration of confidence intervals
showRM(mod1, alpha = 0.99, RM = "VaR", method = "BFGS")
showRM(mod1, alpha = 0.99, RM = "ES", method = "BFGS")

# Effect of changing threshold on xi
xiplot(data,start = 15, end = 200)
# models above line have infinite kurtosis
abline(h = 0.25)


