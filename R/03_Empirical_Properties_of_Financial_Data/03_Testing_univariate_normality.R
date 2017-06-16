## By Marius Hofert

## Testing univariate normality


### Setup ######################################################################

library(nortest) # for cvm.test()
library(ADGofTest) # for ad.test()
library(moments) # for agostino.test(), jarque.test()
library(qrmtools)
library(qrmdata)
library(xts) # for time-series related functions


### 1 Generate data from N(mu, sig^2) and t_nu(mu, sig^2) ######################

n <- 1000 # sample size
mu <- 1 # location
sig <- 2 # scale
nu <- 3 # degrees of freedom
set.seed(271) # set seed (for reproducibility)
X.norm <- rnorm(n, mean = mu, sd = sig) # sample from N(mu, sig^2)
X.t <- mu + sig * rt(n, df = nu) * sqrt((nu-2)/nu) # sample from t_nu(mu, sqrt((nu-2)/nu)*sig^2) (same variance as N(mu, sig^2))
if(FALSE) {
    var(X.norm)
    var(X.t)
    ## => quite apart (but closer for larger n)
}


### 2 Testing for N(mu, sig^2) #################################################

## We treat mu and sig^2 as unknown and estimate them
mu.norm <- mean(X.norm)
mu.t <- mean(X.t)
sig.norm <- sd(X.norm)
sig.t <- sd(X.t)


### 2.1 Tests for general distributions (here applied to the normal) ###########

## Applied to normal data (with estimated mu and sig^2)
(ks <- ks.test(X.norm, y = "pnorm", mean = mu.norm, sd = sig.norm)) # Kolmogorov--Smirnov
(cvm <- cvm.test(X.norm)) # Cramer--von Mises for normality
(ad <- ad.test(X.norm, distr.fun = pnorm, mean = mu.norm, sd = sig.norm)) # Anderson--Darling
stopifnot(min(ks$p.value, cvm$p.value, ad$p.value) >= 0.05) # => no rejections

## Applied to t data (with estimated mu and sig^2)
(ks <- ks.test(X.t, y = "pnorm", mean = mu.t, sd = sig.t)) # Kolmogorov--Smirnov
(cvm <- cvm.test(X.t)) # Cramer--von Mises for normality
(ad <- ad.test(X.t, distr.fun = pnorm, mean = mu.t, sd = sig.t)) # Anderson--Darling
stopifnot(max(ks$p.value, cvm$p.value, ad$p.value) < 0.05) # => all rejections based on 5%


### 2.2 Tests specifically for the normal distribution #########################

## Applied to normal data
(sh <- shapiro.test(X.norm)) # Shapiro--Wilk
(ag <- agostino.test(X.norm)) # D'Agostino's test
(jb <- jarque.test(X.norm)) # Jarque--Bera test
stopifnot(min(sh$p.value, ag$p.value, jb$p.value) >= 0.05) # => no rejections

## Applied to t data
(sh <- shapiro.test(X.t)) # Shapiro--Wilk
(ag <- agostino.test(X.t)) # D'Agostino's test
(jb <- jarque.test(X.t)) # Jarque--Bera test
stopifnot(max(sh$p.value, ag$p.value, jb$p.value) < 0.05) # => all rejections based on 5%


### 2.3 Graphical tests of normality ###########################################

## Applied to normal data
pp_plot(X.norm, FUN = function(q) pnorm(q, mean = mu.norm, sd = sig.norm))
qq_plot(X.norm, FUN = function(p) qnorm(p, mean = mu.norm, sd = sig.norm))

## Applied to t data
pp_plot(X.t, FUN = function(q) pnorm(q, mean = mu.t, sd = sig.t))
qq_plot(X.t, FUN = function(p) qnorm(p, mean = mu.t, sd = sig.t))
## => S-shape => Data seems to come from heavier-tailed distribution than
##    the normal distribution.


### 3 Data application #########################################################

## Prepare risk factor data
data(DJ_const) # constituents data
margin <- c("KO", "MSFT", "INTC", "DIS") # margins considered
DJ.const <- DJ_const['1993-01-01/2000-12-31', margin]
X <- -returns(DJ.const) # compute -log-returns; daily risk-factor changes
X. <- list(daily = X, # daily risk-factor changes
           weekly = apply.weekly(X, FUN = colSums), # weekly risk-factor changes
           monthly = apply.monthly(X, FUN = colSums), # monthly risk-factor changes
           quarterly = apply.quarterly(X, FUN = colSums)) # quarterly risk-factor changes

## Are the risk-factor changes normally distributed?
## Formal tests
pvals <- sapply(X., function(x) apply(x, 2, function(data) jarque.test(data)$p.value))
## => Daily and weekly risk-factor changes are not (univariate) normal

## Q-Q plots
time <- c("daily", "weekly", "monthly", "quarterly")
opar <- par(ask = TRUE) # ask after each plot
for(j in seq_len(ncol(X))) { # for all margins, do...
    for(k in seq_along(time)) { # for daily, weekly, monthly and quarterly data, do...
        qq_plot(X.[[k]][,j], FUN = function(p)
            qnorm(p, mean = mean(X.[[k]][,j]), sd = sd(X.[[k]][,j])),
            main = paste("Q-Q plot for margin",margin[j],"based on",time[k],"data"))
        mtext(paste("p-value:", round(pvals[j,k], 4)), side = 4, line = 1, adj = 0)
    }
}
par(opar)


### 4 Simulation ###############################################################

## Under H0, p-values are uniformly distributed. Let's check that for N(0,1) data
set.seed(271)
pvals.H0 <- sapply(1:200, function(b) jarque.test(rnorm(n))$p.value)
qq_plot(pvals.H0, FUN = qunif)

## This does not hold under H1
set.seed(271)
pvals.H1 <- sapply(1:200, function(b) jarque.test(rt(n, df = nu))$p.value)
qq_plot(pvals.H1, FUN = qunif)
all(pvals.H1 == 0) # all p-values 0
