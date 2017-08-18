## By Alexander J. McNeil

library(xts)
library(qrmtools)
library(qrmdata)

## Now we want to examine the quality of delta and delta-gamma approximations for typical risk-factor changes

## First load some index return and implied volatility data
data(SP500)
data(VIX)
plot(SP500)
plot(VIX)
## Now compute log returns since 2000 and make bivariate dataset
SP500.X <- diff(log(SP500))['1990-01-01/']
VIX.X <- diff(log(VIX))['1990-01-01/']
length(SP500.X)
length(VIX.X)
X <- as.matrix(merge(SP500.X, VIX.X))
plot(X)

## Option characteristics
K <- 100
## Option currently "at the money"
S <- 100
## Interest rate fixed at 3%
r <- 0.03
## Annualized volatility 25%
sigma <- 0.25
## Maturity one year
T <- 1


## Full revaluation (loss operator implicit)
deltat <- 1/250
V_t0 <- Black_Scholes(0, S, r, sigma, K, T, "call")
V_t1 = Black_Scholes(deltat, exp(log(S)+X[,1]), r, exp(log(sigma)+X[,2]), K, T, "call")
L = - (V_t1 - V_t0)
hist(L)

## Delta Approximation; linear loss operator
(tmp <- Black_Scholes_Greeks(0, S, r, sigma, K, T))
delta <- tmp[,"delta"]
theta <- tmp[,"theta"]
vega <- tmp[,"vega"]

LDelta <- -(delta*X[,1]*S + theta*(deltat) + vega*X[,2]*sigma)
plot(L,LDelta)
abline(0, 1, col = 2)

## Delta-Gamma Approximation; quadratic loss operator
gamma <- tmp[,"gamma"]
vanna <- tmp[,"vanna"]
vomma <- tmp[,"vomma"]


LDeltaGamma <- LDelta -0.5*(gamma*(X[,1]*S)^2 + 2*vanna*X[,1]*X[,2]*S*sigma + vomma*(X[,2]*sigma)^2)
plot(L, LDeltaGamma)
abline(0, 1, col = 2)

## Now for an option that is out of the money
## Option characteristics
K <- 100
S <- 70
r <- 0.03
sigma <- 0.25
T <- 1
(tmp <- Black_Scholes_Greeks(0, S, r, sigma, K, T))
delta <- tmp[,"delta"]
theta <- tmp[,"theta"]
vega <- tmp[,"vega"]
gamma <- tmp[,"gamma"]
vanna <- tmp[,"vanna"]
vomma <- tmp[,"vomma"]

V_t0 <- Black_Scholes(0, S, r, sigma, K, T, "call")
V_t1 <- Black_Scholes(deltat, exp(log(S)+X[,1]), r, exp(log(sigma)+X[,2]), K, T, "call")
L <- -(V_t1 - V_t0)
LDelta <- -(delta*X[,1]*S + theta*(deltat) + vega*X[,2]*sigma)
plot(L, LDelta)
abline(0, 1, col = 2)
LDeltaGamma <- LDelta -0.5*(gamma*(X[,1]*S)^2 + 2*vanna*X[,1]*X[,2]*S*sigma + vomma*(X[,2]*sigma)^2)
plot(L, LDeltaGamma)
abline(0, 1, col = 2)


## THE EXAMPLE IN TEXTBOOK (Example 9.1)

t <- 0
T <- 1
S <- 110
K <- 100
r <- 0.02
sigma <- 0.2

delta.S <- 0.05
delta.sigma <- 0.02
delta.t <- 1/250
delta.r <- 0.001

(tmp <- Black_Scholes_Greeks(0, S, r, sigma, K, T))
delta <- tmp[,"delta"]
theta <- tmp[,"theta"]
vega <- tmp[,"vega"]
rho <- tmp[,"rho"]
gamma <- tmp[,"gamma"]
vanna <- tmp[,"vanna"]
vomma <- tmp[,"vomma"]

(delta.term <- delta*S*delta.S)
(vega.term <- vega*delta.sigma)
(theta.term <- theta*delta.t)
(rho.term <- rho*delta.r)

## ignore rho from now on
## linear.loss = theta.term + rho.term + vega.term
linear.loss <- theta.term + vega.term
linear.loss

Vt <- S*delta - Black_Scholes(t, S, r, sigma, K, T, "call")
Snew <- exp(log(S)+delta.S)
Vt1 <- Snew*delta - Black_Scholes((t+1/250), Snew, r, sigma+delta.sigma, K, T, "call")
Vt
Vt1
trueloss <- (Vt-Vt1)
trueloss


(gamma.term <- 0.5*gamma*S^2*delta.S^2)
linear.loss + gamma.term

(vanna.term <- vanna*S*delta.S*delta.sigma)
linear.loss + gamma.term + vanna.term

(vomma.term <- 0.5*vomma*delta.sigma^2)
quadratic.loss <- linear.loss + gamma.term + vanna.term + vomma.term
quadratic.loss

c(trueloss, linear.loss, quadratic.loss)

(trueloss-linear.loss) / trueloss
(trueloss-quadratic.loss) / trueloss


