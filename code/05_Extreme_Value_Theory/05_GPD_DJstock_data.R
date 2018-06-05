## By Alexander McNeil and Marius Hofert

library(xts)
library(qrmdata)
library(qrmtools)


### 1 Exploring the data #######################################################

data("DJ_const")
## Create log returns
DJ.r <- diff(log(DJ_const))[-1,]
## Aggregate by week
DJ.w <- apply.weekly(DJ.r, FUN = colSums)
## Convert to weekly percentage losses
DJ.wp <- -(exp(DJ.w)-1)*100
## Take data from start of millenium
DJ.wp <- DJ.wp['2000-01-01/']
plot.zoo(DJ.wp[,1:12], type = "h")

## Show mean excess plots
op <- par(mfrow = c(3, 4), mar = c(2,2,2,1))
for (i in 1:12){
    data <- DJ.wp[,i]
    mean_excess_plot(data[data>0], xlab = "", ylab = "",
                     main = names(DJ_const)[i])
}
par(op)

## Select Apple share price
data <- DJ.wp[,"AAPL"]


### 2 Find the threshold #######################################################

## Mean excess plot of losses omitting gains
mean_excess_plot(data[data>0])
u <- 6
abline(v = u)

## Effect of changing threshold on xi
GPD_shape_plot(data)
abline(v = u)
abline(h = 0.25) # models above line have infinite kurtosis (E(X^k) = Inf <=> xi >= 1/k; k = 4 here)


### 3 Fit a GPD to the excesses ################################################

## Fit GPD model to excesses via MLE
exceed <- data[data > u] # exceedances
excess <- exceed - u # excesses
(fit <- fit_GPD_MLE(excess)) # MLE
shape.u <- fit$par[["shape"]]
scale.u <- fit$par[["scale"]]


### 4 Visually compare the sample excess df and the fitted (GPD) excess df #####

## Plot empirical excess df vs fitted GPD
res <- edf_plot(excess, log = "x")
z <- tail(res$t, n = -1)
lines(z, pGPD(z, shape = shape.u, scale = scale.u)) # fitted GPD

## Plot empirical exceedance df vs shifted fitted GPD
res <- edf_plot(exceed, log = "x")
z <- tail(res$t, n = -1)
lines(z, pGPD(z-u, shape = shape.u, scale = scale.u)) # shifted fitted GPD

## Corresponding Q-Q plots (more meaningful)
qf <- function(p) # quantile function of df
    qGPD(p, shape = shape.u, scale = scale.u)
qq_plot(excess, FUN = qf)
qq_plot(exceed, FUN = function(p) u + qf(p))


### 5 Compute semi-parametric risk measure estimators ##########################

## VaR_alpha, ES_alpha for two alphas
alpha <- c(0.99, 0.995)
(VaR <- VaR_GPDtail(alpha, threshold = u, p.exceed = mean(data > u),
                    shape = shape.u, scale = scale.u))
(ES  <-  ES_GPDtail(alpha, threshold = u, p.exceed = mean(data > u),
                    shape = shape.u, scale = scale.u))


### 6 Semi-parametric Smith estimator including VaR_0.99 #######################

## Empirical tail probabilities with Smith estimator overlaid
tail_plot(data, threshold = u, shape = shape.u, scale = scale.u)
abline(h = 1-alpha[1], v = VaR[1], lty = 2) # 0.99, VaR_0.99
abline(h = 1-alpha[1], v = ES[1],  lty = 2) # 0.99, ES_0.99

## A version including confidence intervals for VaR_0.99 and ES_0.99
fit. <- QRM::fit.GPD(data, u) # need 'QRM' for that
QRM::showRM(fit., alpha = 0.99, RM = "VaR", method = "BFGS")
QRM::showRM(fit., alpha = 0.99, RM = "ES",  method = "BFGS")
