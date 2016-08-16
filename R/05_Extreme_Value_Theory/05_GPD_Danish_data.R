## By Alexander McNeil and Marius Hofert


### Setup ######################################################################

library(QRM)


### EVT related plots ##########################################################

## Sample mean excess plot
plot(danish)
MEplot(danish) # sample mean excess plot
u10 <- 10 # threshold 1
u20 <- 20 # threshold 2
abline(v = u10)
abline(v = u20)

## Fit GPD model
str(fit.GPD)
(mod.u10 <- fit.GPD(danish, threshold = u10))

## Risk measure estimates
(RM.u10 <- RiskMeasures(mod.u10, c(0.99, 0.995)))

## Plot the empirical excess df and the theoretical excess df (GPD)
plotFittedGPDvsEmpiricalExcesses(danish, threshold = u10) # hat{F}_{u,n}(x-u) vs G_{hat{xi},hat{beta}}(x-u)
legend("bottomright", bty = "n", lty = 0:1, pch = c(19, NA), col = c("blue", "black"),
       legend = c("Empirical excess df", "Theoretical excess df (GPD)"))

## Semi-parametric Smith/tail estimator and implied risk measures
plotTail(mod.u10, main = "Semi-parametric Smith/tail estimator") # hat(bar(F))(x), x >= u
showRM(mod.u10, alpha = 0.99, RM = "VaR", method = "BFGS") # with VaR estimate and CIs
showRM(mod.u10, alpha = 0.99, RM = "ES",  method = "BFGS") # with ES estimate and CIs

## Effect of changing the threshold on xi
xiplot(danish) # => variance increases for larger thresholds


### 2nd threshold analysis #####################################################

## Fit GPD model
(mod.u20 <- fit.GPD(danish, threshold = u20))

## Risk measure estimates
(RM.u20 <- RiskMeasures(mod.u20, c(0.99, 0.995)))

## Plot the sample excess df and the theoretical excess df (GPD)
plotFittedGPDvsEmpiricalExcesses(danish, threshold = u20)
legend("bottomright", bty = "n", lty = 0:1, pch = c(19, NA), col = c("blue", "black"),
       legend = c("Empirical excess df", "Theoretical excess df (GPD)"))

## Semi-parametric Smith/tail estimator and implied risk measures
plotTail(mod.u20, main = "Semi-parametric Smith/tail estimator") # hat(bar(F))(x), x >= u
showRM(mod.u20, alpha = 0.99, RM = "VaR", method = "BFGS") # with VaR estimate and CIs
showRM(mod.u20, alpha = 0.99, RM = "ES",  method = "BFGS") # with ES estimate and CIs
