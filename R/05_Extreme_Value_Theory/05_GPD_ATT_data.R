## by Alexander McNeil

require(QRM)

# Get ATT share price and plot it
ATT.price <- DJ[,"T"]
plot(ATT.price)
# Compute daily log returns
ATT.d <- returns(ATT.price,method="continuous")
plot(ATT.d)

# Aggregate to get weekly log returns
by <- timeSequence(from = start(ATT.d),  to = end(ATT.d), by = "week")
ATT.w <- aggregate(ATT.d, by, sum)

# Transform to percentage losses
ATT.w <- -(exp(ATT.w)-1)*100
plot(ATT.w,type="h")

# Mean excess plot of losses omitting gains
MEplot(ATT.w[ATT.w>0])
abline(v=2.75)

# Fit GPD model above 2.75%
# Note that numbers are a little different to those reported 
# in the book. 

mod1 <- fit.GPD(ATT.w,2.75)
mod1
plotFittedGPDvsEmpiricalExcesses(ATT.w, threshold = 2.75)

# Compute VaR and expected shortfall
RMs1 <- RiskMeasures(mod1, c(0.99,0.995))
RMs1
plotTail(mod1)
showRM(mod1, alpha = 0.99, RM = "VaR", method = "BFGS")
abline(h=0.01)
showRM(mod1, alpha = 0.99, RM = "ES", method = "BFGS")

# Effect of changing threshold on xi
xiplot(ATT.w,start=15,end=200)


