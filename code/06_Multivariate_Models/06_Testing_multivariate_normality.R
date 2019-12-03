## By Marius Hofert

## Tests of multivariate normality


### Setup ######################################################################

library(qrmtools)
library(nvmix) # for rNorm() and rStudent()
library(moments) # for jarque.test()
library(xts) # for time-series related functions


### 1 Generate data ############################################################

n <- 1000 # sample size
d <- 3 # dimension
mu <- 1:d # location vector
sd <- sqrt(1:d) # standard deviations for scale matrix
tau <- 0.5 # (homogeneous) Kendall's tau
rho <- sin(tau * pi/2) # (homogeneous) correlation coefficient
Sigma <- matrix(rho, nrow = d, ncol = d) # scale matrix
diag(Sigma) <- rep(1, d) # => now Sigma is a correlation matrix
Sigma <- Sigma * (sd %*% t(sd)) # covariance matrix
nu <- 3 # degrees of freedom
set.seed(271) # set seed (for reproducibility)
X.N <- rNorm(n, loc = mu, scale = Sigma) # sample from N(mu, Sigma)
X.t <- rStudent(n, df = nu, loc = mu, scale = ((nu-2)/nu) * Sigma) # sample from t_d(nu, mu, ((nu-2)/nu) * Sigma) (same covariance as N(mu, Sigma))
X.t.N <- sapply(1:d, function(j) qnorm(rank(X.t[,j]) / (nrow(X.t) + 1),
                                       mean = mu[j], sd = sqrt(Sigma[j, j]))) # t dependence with normal margins; see Chapter 7
if(FALSE) {
    cov(X.N)
    cov(X.t)
    ## => quite apart, but closer for larger n
}


### 2 Testing for N(mu, Sigma) #################################################

## We treat mu and Sigma as unknown and estimate them
mu.N <- mean(X.N)
mu.t <- mean(X.t)
Sig.N <- cov(X.N)
Sig.t <- cov(X.t)


### 2.1 For N(mu, Sigma) data ##################################################

### 2.1.1 Formal tests

## Univariate normality; see also 03_Testing_univariate_normality.R
stopifnot(apply(X.N, 2, function(x) shapiro.test(x)$p.value) > 0.05) # Shapiro--Wilk
stopifnot(apply(X.N, 2, function(x) jarque.test(x) $p.value) > 0.05) # Jarque--Bera
## Note: Careful with 'multiple testing'

## Joint normality
ADmaha2 <- maha2_test(X.N) # Anderson--Darling test of squared Mahalanobis distances being chi_d^2
mardiaK <- mardia_test(X.N) # Mardia's kurtosis test of joint normality based on squared Mahalanobis distances
mardiaS <- mardia_test(X.N, type = "skewness") # Mardia's skewness test of joint normality based on Mahalanobis angles
stopifnot(mardiaS$p.value > 0.05, mardiaK$p.value > 0.05, ADmaha2$p.value > 0.05)


### 2.1.2 Graphical tests

## Univariate normality
qq_plot(X.N[,1], FUN = qnorm, method = "empirical") # margin 1
qq_plot(X.N[,2], FUN = qnorm, method = "empirical") # margin 2
qq_plot(X.N[,3], FUN = qnorm, method = "empirical") # margin 3

## Joint normality
pairs(X.N, gap = 0, pch = ".") # visual assessment
D2.d <- mahalanobis(X.N, center = colMeans(X.N), cov = cov(X.N)) # squared Mahalanobis distances
qq_plot(D2.d, FUN = function(p) qchisq(p, df = d)) # => no significant departure visible


### 2.2 For t_d(nu, mu, Sigma) data ############################################

### 2.2.1 Formal tests

## Univariate normality
stopifnot(apply(X.t, 2, function(x) shapiro.test(x)$p.value) <= 0.05) # Shapiro--Wilk
stopifnot(apply(X.t, 2, function(x) jarque.test(x) $p.value) <= 0.05) # Jarque--Bera
## Note: Careful with 'multiple testing'

## Joint normality
ADmaha2 <- maha2_test(X.t) # Anderson--Darling test of squared Mahalanobis distances being chi_d^2
mardiaK <- mardia_test(X.t) # Mardia's kurtosis test of joint normality based on squared Mahalanobis distances
mardiaS <- mardia_test(X.t, type = "skewness") # Mardia's skewness test of joint normality based on Mahalanobis angles
stopifnot(mardiaS$p.value <= 0.05, mardiaK$p.value <= 0.05, ADmaha2$p.value <= 0.05)


### 2.2.2 Graphical tests

## Univariate normality
qq_plot(X.t[,1], FUN = qnorm, method = "empirical") # margin 1
qq_plot(X.t[,2], FUN = qnorm, method = "empirical") # margin 2
qq_plot(X.t[,3], FUN = qnorm, method = "empirical") # margin 3

## Joint normality
pairs(X.t, gap = 0, pch = ".") # visual assessment
D2.d <- mahalanobis(X.t, center = colMeans(X.t), cov = cov(X.t)) # squared Mahalanobis distances
qq_plot(D2.d, FUN = function(p) qchisq(p, df = d)) # => departure clearly visible


### 2.3 For t_d(nu, 0, P) dependence but N(0, 1) margins #######################

### 2.3.1 Formal tests

## Univariate normality
stopifnot(apply(X.t.N, 2, function(x) shapiro.test(x)$p.value) > 0.05) # Shapiro--Wilk
stopifnot(apply(X.t.N, 2, function(x) jarque.test(x) $p.value) > 0.05) # Jarque--Bera

## Joint normality
ADmaha2 <- maha2_test(X.t.N) # Anderson--Darling test of squared Mahalanobis distances being chi_d^2
mardiaK <- mardia_test(X.t.N) # Mardia's kurtosis test of joint normality based on squared Mahalanobis distances
mardiaS <- mardia_test(X.t.N, type = "skewness") # Mardia's skewness test of joint normality based on Mahalanobis angles
stopifnot(mardiaS$p.value <= 0.05, mardiaK$p.value <= 0.05, ADmaha2$p.value <= 0.05)


### 2.3.2 Graphical tests

## Univariate normality
qq_plot(X.t.N[,1], FUN = qnorm, method = "empirical") # margin 1
qq_plot(X.t.N[,2], FUN = qnorm, method = "empirical") # margin 2
qq_plot(X.t.N[,3], FUN = qnorm, method = "empirical") # margin 3

## Joint normality
pairs(X.t.N, gap = 0, pch = ".") # visual assessment
D2.d <- mahalanobis(X.t.N, center = colMeans(X.t.N), cov = cov(X.t.N)) # squared Mahalanobis distances
qq_plot(D2.d, FUN = function(p) qchisq(p, df = d)) # => departure clearly visible
