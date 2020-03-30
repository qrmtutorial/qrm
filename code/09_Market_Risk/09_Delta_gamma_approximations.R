## By Alexander J. McNeil and Marius Hofert

## Illustration of the quality of delta (linear) and delta-gamma (quadratic) approximations
## for typical risk-factor changes


### Setup ######################################################################

library(xts)
library(qrmtools)
library(qrmdata)


### 1 Data preparation #########################################################

## Load some index return data and implied volatility data
data(SP500)
data(VIX)

## Visually check
plot(SP500)
plot(VIX)

## Now compute log-returns and make a bivariate dataset
X1 <- returns(SP500)['1990-01-01/']
X2 <- returns(VIX/100, method = "diff")['1990-01-01/'] # note: VIX is in % (matters for 'diff')
X. <- merge(X1, X2)
any(is.na(X.))
X <- as.matrix(na.fill(X., fill = "extend"))
colnames(X) <- c("SP500", "VIX")
plot(X)


### 2 Losses and their approximations ##########################################

### 2.1 Option at the money ####################################################

## Option characteristics
T <- 1 # maturity one year
K <- 100 # strike
S <- 100 # stock price now => option currently "at the money"
r <- 0.03 # interest rate 3%
sigma <- 0.25 # annualized volatility 25%

## Full revaluation with risk factor changes
Delta.t <- 1/250
V.t0 <- Black_Scholes(0, S = S, r = r, sigma = sigma, K = K, T = T, type = "call") # single value
V.t1 <- Black_Scholes(Delta.t, S = exp(log(S) + X[,"SP500"]), # current value and vector of risk factor changes
                      r = r,
                      sigma = exp(log(sigma) + X[,"VIX"]), # current value and vector of risk factor changes
                      K = K, T = T, type = "call") # time series
L <- -(V.t1 - V.t0) # losses (implicit loss operator)

## Visual check
hist(L, probability = TRUE); box()

## Delta approximation (linear loss operator)
(Greeks <- Black_Scholes_Greeks(0, S = S, r = r, sigma = sigma, K = K, T = T))
L.Delta <- -(Greeks[["delta"]] * X[,"SP500"] * S +
             Greeks[["theta"]] * Delta.t +
             Greeks[["vega"]]  * X[,"VIX"]   * sigma)

## Visual comparison with full revaluation
plot(L, L.Delta)
abline(0, 1, col = "royalblue3")

## Delta-Gamma approximation (quadratic loss operator)
L.Delta.Gamma <- L.Delta - 0.5 * (Greeks[["gamma"]] * (X[,"SP500"] * S)^2 +
                                  2 * Greeks[["vanna"]] * X[,"SP500"] * X[,"VIX"] * S * sigma +
                                  Greeks[["vomma"]] * (X[,"VIX"] * sigma)^2)

## Visual comparison with full revaluation
plot(L, L.Delta.Gamma)
abline(0, 1, col = "royalblue3")


### 2.2 Option out of the money ################################################

## Option characteristics
S <- 70 # stock price now => option currently "out of the money"

## Full revaluation with risk factor changes
V.t0 <- Black_Scholes(0, S = S, r = r, sigma = sigma, K = K, T = T, type = "call") # single value
V.t1 <- Black_Scholes(Delta.t, S = exp(log(S) + X[,"SP500"]), # current value and vector of risk factor changes
                      r = r,
                      sigma = exp(log(sigma) + X[,"VIX"]), # current value and vector of risk factor changes
                      K = K, T = T, type = "call") # time series
L <- -(V.t1 - V.t0) # losses (implicit loss operator)

## Delta approximation (linear loss operator)
(Greeks <- Black_Scholes_Greeks(0, S = S, r = r, sigma = sigma, K = K, T = T))
L.Delta <- -(Greeks[["delta"]] * X[,"SP500"] * S +
             Greeks[["theta"]] * Delta.t +
             Greeks[["vega"]]  * X[,"VIX"]   * sigma)

## Visual comparison with full revaluation
plot(L, L.Delta)
abline(0, 1, col = "royalblue3")

## Delta-Gamma approximation (quadratic loss operator)
L.Delta.Gamma <- L.Delta - 0.5 * (Greeks[["gamma"]] * (X[,"SP500"] * S)^2 +
                                  2 * Greeks[["vanna"]] * X[,"SP500"] * X[,"VIX"] * S * sigma +
                                  Greeks[["vomma"]] * (X[,"VIX"] * sigma)^2)

## Visual comparison with full revaluation
plot(L, L.Delta.Gamma)
abline(0, 1, col = "royalblue3")


### 3 MFE (2015, Example 9.1) ##################################################

## Option characteristics
t <- 0 # time now
T <- 1 # maturity one year
K <- 100 # strike
S <- 110 # stock price now => option currently in the money
r <- 0.02 # interest rate 2%
sigma <- 0.2 # annualized volatility 20%

## Ingredients for loss operators
delta.S   <- 0.05 # change in stock price = X_{t+1,1}
delta.sig <- 0.02 # change in implied volatility = X_{t+1,2}
Delta.t  <- 1/250
Delta.rho <- 0.001
(Greeks <- Black_Scholes_Greeks(0, S = S, r = r, sigma = sigma, K = K, T = T))

## Full revaluation with one new risk factor change (provided by specifying 'S.new')
V.t0 <- S * Greeks[["delta"]] - Black_Scholes(t, S = S, r = r, sigma = sigma, K = K, T = T)
S.new <- exp(log(S) + delta.S)
V.t1 <- S.new * Greeks[["delta"]] - Black_Scholes(t + Delta.t, S = S.new, r = r,
                                                  sigma = sigma + delta.sig, K = K, T = T)
L <- -(V.t1 - V.t0) # losses (implicit loss operator)

## Delta approximation
delta.term <- Greeks[["delta"]] * S * delta.S
vega.term  <- Greeks[["vega"]]  * delta.sig
theta.term <- Greeks[["theta"]] * Delta.t
rho.term   <- Greeks[["rho"]]   * Delta.rho
## => L.Delta = theta.term + rho.term + vega.term (but we ignore the rho term)
(L.Delta <- theta.term + vega.term)

## Delta-Gamma approximation
gamma.term <- 0.5 * Greeks[["gamma"]] * S^2 * delta.S^2
vanna.term <- Greeks[["vanna"]] * S * delta.S * delta.sig
vomma.term <- 0.5 * Greeks[["vomma"]] * delta.sig^2
L.Delta + gamma.term
L.Delta + vanna.term
(L.Delta.Gamma <- L.Delta + gamma.term + vanna.term + vomma.term)

## Comparison
c(L, L.Delta, L.Delta.Gamma)
abs(L - L.Delta)/abs(L) # relative error of L.Delta
abs(L - L.Delta.Gamma)/abs(L) # relative error of L.Delta.Gamma
