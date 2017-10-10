## By Alexander McNeil and Marius Hofert


### Setup ######################################################################

library(QRM)


### 1 Find the threshold(s) and fitting GPD(s) #################################

## Plot
plot(danish)

## Sample mean excess plot (with two plausible thresholds)
MEplot(danish) # sample mean excess plot; omits largest three points ('omit = 3')
u10 <- 10 # threshold 1
u20 <- 20 # threshold 2
abline(v = u10)
abline(v = u20)

## Effect of changing the threshold on xi
xiplot(danish) # => variance increases for larger thresholds

## Fit GPD models
str(fit.GPD) # => gets the original data as 'data' and 'threshold' or 'nextremes'
(mod.u10 <- fit.GPD(danish, threshold = u10))
(mod.u20 <- fit.GPD(danish, threshold = u20))


### 2 Plot the sample excess df hat{F}_{u,n} and the fitted excess df (GPD) at x-u

## Understanding plotFittedGPDvsEmpiricalExcesses()
plotFittedGPDvsEmpiricalExcesses # => compare with GPD-based tail bar{F}; note: mod$data = data[data > u]
## Note: F_u(x-u) = P(X-u <= x-u | X > u) = P(X <= x | X > u)
##       => hat{F}_{u,n}(x-u) = hat{F}_n(x) where hat{F}_n is based on all X > u
##       => that's what's internally used, see 'edf(mod$data)'
stopifnot(mod.u10$data == danish[danish > u10])

## Plot the sample excess df hat{F}_{u,n} at x-u and the fitted excess df (GPD)
## at x-u for x >= u
plotFittedGPDvsEmpiricalExcesses(danish, threshold = u10) # hat{F}_{u,n}(x-u) vs G_{hat{xi},hat{beta}}(x-u)
legend("bottomright", bty = "n", lty = 0:1, pch = c(19, NA), col = c("blue", "black"),
       legend = c("Empirical excess df", "Theoretical excess df (GPD)"))

## Threshold 2
plotFittedGPDvsEmpiricalExcesses(danish, threshold = u20)
legend("bottomright", bty = "n", lty = 0:1, pch = c(19, NA), col = c("blue", "black"),
       legend = c("Empirical excess df", "Theoretical excess df (GPD)"))


### 3 Compute semi-parametric risk measure estimates ###########################

## Risk measure estimates
RiskMeasures # => compare with GPD-based VaR and ES formulas
(RM.u10 <- RiskMeasures(mod.u10, p = c(0.99, 0.995)))

## Threshold 2
(RM.u20 <- RiskMeasures(mod.u20, p = c(0.99, 0.995)))


### 4 Compute semi-parametric tail estimators and implied risk measures ########

## Semi-parametric Smith/tail estimator and implied risk measures
opar <- par(mar = c(5, 4, 4, 2) + c(0, 1, 0, 0))
plotTail(mod.u10, main = "Semi-parametric Smith/tail estimator") # hat(bar(F))(x), x >= u
showRM(mod.u10, alpha = 0.99, RM = "VaR", method = "BFGS") # with VaR estimate and CIs
showRM(mod.u10, alpha = 0.99, RM = "ES",  method = "BFGS") # with ES  estimate and CIs
par(opar)

## Threshold 2
opar <- par(mar = c(5, 4, 4, 2) + c(0, 1, 0, 0))
plotTail(mod.u20, main = "Semi-parametric Smith/tail estimator") # hat(bar(F))(x), x >= u
showRM(mod.u20, alpha = 0.99, RM = "VaR", method = "BFGS") # with VaR estimate and CIs
showRM(mod.u20, alpha = 0.99, RM = "ES",  method = "BFGS") # with ES  estimate and CIs
par(opar)
## => CIs quite a bit wider (as we have less data)
