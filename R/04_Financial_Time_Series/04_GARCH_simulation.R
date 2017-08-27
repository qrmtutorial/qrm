## By Alexander McNeil and Marius Hofert

## Simulating paths, extracting components and plotting of a GARCH(1,1) model


### Setup ######################################################################

library(rugarch)
library(qrmtools)


### 1 Specify the GARCH(1,1) model #############################################

## GARCH(1,1) model specification
(uspec <- ugarchspec(variance.model = list(model = "sGARCH", # standard GARCH
                                           garchOrder = c(1, 1)), # GARCH(1,1)
                     mean.model = list(armaOrder = c(0, 0), # no ARMA part
                                       include.mean = FALSE), # no mean included
                     distribution.model = "norm", # normal innovations
                     fixed.pars = list(omega = 0.02, # cond. var. constant (= alpha_0)
                                       alpha1 = 0.15, # alpha_1
                                       beta1 = 0.8))) # beta_1


### 2 Generate paths from a specified model ####################################

## Generate two realizations of length 2000
m <- 2 # number of paths
n <- 2000 # sample size
set.seed(271) # set seed (for reproducibility)
(paths <- ugarchpath(uspec, n.sim = n, n.start = 50, # burn-in sample size
                     m.sim = m))
str(paths) # S4 object of class 'uGARCHpath'
str(paths@path) # simulated processes X_t, volatilities sigma_t and unstandardized residuals epsilon_t

## We can use getClass() to see more information about such objects
getClass("uGARCHpath")

## We can also use getSlots() to see the composition
getSlots("uGARCHpath")

## EXERCISE: Try simulating paths with different innovation distributions (e.g. Student t)


### 3 Plotting uGARCHpath objects ##############################################

### Extract the simulated series
X <- fitted(paths) # simulated process X_t = mu_t + epsilon_t for epsilon_t = sigma_t * Z_t
sig <- sigma(paths) # volatilities sigma_t (conditional standard deviations)
eps <- paths@path$residSim # unstandardized residuals epsilon_t = sigma_t * Z_t
## Note: There are no extraction methods for the unstandardized residuals epsilon_t
## Sanity checks (=> fitted() and sigma() grab out the right quantities)
stopifnot(all.equal(X,   paths@path$seriesSim, check.attributes = FALSE),
          all.equal(sig, paths@path$sigmaSim,  check.attributes = FALSE))

## Plot the simulated paths (X_t)
plot(X[,1], type = "l", ylab = expression(X[1]))
plot(X[,2], type = "l", ylab = expression(X[2]))

## Plot the corresponding volatilities (sigma_t)
plot(sig[,1], type = "h", ylab = expression(sigma[1]))
plot(sig[,2], type = "h", ylab = expression(sigma[2]))

## Plot the corresponding unstandardized residuals (epsilon_t)
plot(eps[,1], type = "l", ylab = expression(epsilon[1]))
plot(eps[,2], type = "l", ylab = expression(epsilon[2]))

## Compute the standardized residuals (Z_t)
Z <- eps/sig
plot(Z[,1], type = "p", ylab = expression(Z[1]))
plot(Z[,2], type = "p", ylab = expression(Z[2]))
## Check their N(0,1) distribution
qq_plot(Z[,1])
qq_plot(Z[,2])
## Check their joint distribution
plot(Z, xlab = expression(Z[1]), ylab = expression(Z[2]))
## Note: By doing all of this backwards, one can put different joint distributions
##       on the standard residuals. This will become clear after Chapters 6 and 7.

## There are also special plotting functions for 'uGARCHpath' objects:
plot(paths, which = 1) # plots sig[,1] (no matter how many 'paths')
plot(paths, which = 2) # plots X[,1]
plot(paths, which = 3) # plots kernel density estimate of sig[,1]
plot(paths, which = 4) # plots kernel density estimate of X[,1]
## How to see the documentation of the plot function
showMethods("plot")
getMethod("plot", signature = c(x = "uGARCHpath", y = "missing"))


### 4 Does the simulated series conform to our univariate stylized facts? ######

## We (only) consider the first path here
acf(X[,1], main = expression(X[1])) # => no serial correlation
acf(abs(X[,1]), main = expression(group("|",X[1],"|"))) # => profound serial correlation
qq_plot(X[,1], method = "empirical", main = "Normal Q-Q plot") # => heavier tailed than normal
shapiro.test(X[,1]) # Normal rejected
