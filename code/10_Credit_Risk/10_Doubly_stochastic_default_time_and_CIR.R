## By Ruediger Frey

## Illustration for doubly stochastic random times using the CIR model

## We use MC simulation to generate default times via threshold simulation.
## Moreover, we sample (from the loss distribution of a corporate bond) the loss
## distribution of a corporate zero coupon bond and we do bond-option
## pricing. For the latter two applications we have to distinguish carefully
## between real-world quantities (under P) and risk-neutral quantities (under Q).
## Theory and notation are as in McNeil, Frey, Embrechts (2015, Chapter 10).


### Setup ######################################################################

library(sde)
library(qrmtools)


### 1 Simulate default times ###################################################

### 1.1 Auxiliary functions ####################################################

##' @title Computing E[exp(-int_0^tau rho.0 + rho.1* Psi_s ds) | Psi_0 = psi.0]
##         for Psi CIR with kappa, theta, sigma
##' @param psi.0 numeric providing the initial value
##' @param tau numeric providing the time to maturity
##' @param kappa CIR spead of mean reversion
##' @param theta CIR level of mean reversion
##' @param sigma CIR volatility
##' @param rho.0 see above
##' @param rho.1 see above
##' @return E[exp(-int_0^tau rho.0 + rho.1* Psi_s ds) | Psi_0 = psi.0]
##' @author Ruediger Frey
affine_transformation <- function(psi.0, tau, kappa, theta, sigma, rho.0 = 0, rho.1 = 1)
{
    gamma <- sqrt(kappa^2 + 2 * sigma^2 * rho.1)
    denom <- gamma - kappa + exp(gamma * tau) * (gamma + kappa)
    num.for.A <- 2 * gamma * exp(tau * (gamma + kappa) / 2)
    alpha <- -rho.0 * tau + 2 * (kappa * theta / sigma^2) * log(num.for.A / denom)
    beta <- (-2 * rho.1 * (exp(gamma * tau) - 1)) / denom
    exp(alpha + beta * psi.0)
}

##' @title Threshold simulation for doubly stochastic default times (under P)
##' @param gamma.path vector containing the path of hazard rate
##' @param Delta.t the time step used when generating gamma.path
##' @return default indicator (TRUE if tau < T, FALSE otherwise) and
##'         default time (if tau < T)
##' @author Ruediger Frey
threshold_sim <- function(gamma.path, Delta.t)
{
    n.steps <- length(gamma.path)-1
    E <- rexp(1)
    t <- 1
    Gamma <- 0
    while((t <= n.steps) && (Gamma < E)) {
        Gamma <- Gamma + Delta.t * mean(gamma.path[t], gamma.path[t+1])
        t <- t + 1
    }
    if (Gamma >= E) {
        tau <- (t-0.5) * Delta.t
        def.ind <- TRUE
    } else {
        tau <- Inf
        def.ind <- FALSE
    }
    list(indicator = def.ind, time = tau, threshold = E)
}


### 1.2 Recursive default time simulation ######################################

## Set 1 of model parameters
## We consider a constant risk-free rate and gamma.P as rho.P * Psi (to allow
## for different bonds with varying credit quality). Under P, Psi is CIR with
## parameters:
psi.0 <- 0.05 # initial value
theta.P <- 0.05 # long-term level
kappa.P <- 0.5 # speed of mean reversion
sigma <- 0.2 # volatility (identical under P and Q)
rho.P <- 1 # risk loading
stopifnot(2 * kappa.P * theta.P - sigma^2 > 0) # condition for strict positivity of CIR

## Simulation of tau
horizon <- 1
Delta.t <- 1/250 # daily time step
n.steps <- round(horizon/Delta.t)

## First we draw a trajectory/threshold which leads to default
set.seed(19) # for reproducibility
psi.path <- sde.sim(X0 = psi.0, theta = c(kappa.P * theta.P, kappa.P, sigma),
                    model = "CIR", T = horizon, N = n.steps)
(sim.result <- threshold_sim(psi.path, Delta.t))

Gamma <- rep(0, n.steps + 1)
for (t in 1:n.steps)
    Gamma[t+1] <- Gamma[t] + Delta.t * mean(psi.path[t], psi.path[t+1])
times <- (0:(n.steps)) * Delta.t
plot(times, Gamma, type = "l", xlab = "t", ylab = expression(Gamma[t]),
     ylim = range(Gamma), main = "Hazard process for the simulated CIR")
abline(h = sim.result$threshold)
abline(v = sim.result$time)

## Now a trajectory which does not lead to default (essentially due to large threshold)
set.seed(45)
psi.path <- sde.sim(X0 = psi.0, theta = c(kappa.P*theta.P, kappa.P, sigma),
                    model = "CIR", T = horizon, N = n.steps)
(sim.result <- threshold_sim(psi.path, Delta.t))

Gamma <- rep(0,n.steps+1)
for (t in 1:n.steps)
    Gamma[t+1] <- Gamma[t] + Delta.t * mean(psi.path[t], psi.path[t+1])
times <- (0:(n.steps)) * Delta.t
plot(times, Gamma, type = "l", ylim = c(0, sim.result$threshold),
     xlab = "t", ylab = expression(Gamma[t]),
     main = "Hazard process for the simulated CIR")
abline(h = sim.result$threshold)
abline(v = sim.result$time)

## Compare simulated and computed default probabilities
Delta.t <- 1/25 # biweekly time step
n.steps <- round(horizon/Delta.t)
n.paths <- 5000

tau <- vector("numeric", n.paths)
default <- vector("logical", n.paths)
for (i in 1:n.paths) {
    psi.path <- sde.sim(X0 = psi.0, theta = c(kappa.P*theta.P, kappa.P, sigma),
                        model = "CIR", T = horizon, N = n.steps)
    sim.result <- threshold_sim(psi.path, Delta.t)
    if(i %% 500 == 0) cat("path ", i, "\n")
    default[i] <- sim.result$indicator
    tau[i] <- sim.result$time
}
mean(default) # simulated default probability
affine_transformation(psi.0, horizon, kappa.P, theta.P, sigma) # correct value


### 2 Issues related to bond prices ############################################

## Set 2 of model parameters
## For simplicity: LGD = 1 (pure survival claim), dynamics under Q (risk-neutral measure).
## We assume gamma.Q = rho.Q * Psi (typically rho.Q > rho.Psi); the Brownian
## motions are related by dWQ.t = dWP.t + lambda sqrt(Psi_t) dt for risk premium
## lambda. This implies kappa.Q = kappa.P+ sigma *lambda, theta.Q = theta.P *
## (kappa.P/kappa.Q). The price in t of a survival claim with maturity date T.mat
## is thus 1_{tau >t} E.Q(exp(- int_0^T.mat r + rho.Q Psi_s ds | \Psi_t))
psi.0 <- 0.02
theta.P <- 0.02
kappa.P <- 0.75
sigma <- 0.15 # identical under P and Q
lambda <- -0.1 # lambda < 0 means higher mean reversion level and lower mean reversion speed under Q
kappa.Q <- kappa.P + lambda*sigma
theta.Q <- theta.P * (kappa.P / kappa.Q)
rho.P <- 1
rho.Q <- rho.P * 1.5 # default risk premium
stopifnot(2*kappa.Q*theta.Q - sigma^2 > 0,
          2*kappa.P*theta.P - sigma^2 > 0) # condition for strict positivity of CIR
T.mat <- 2 # maturity of bond
r <- 0.01 # risk-free interest rate


### 2.1 VaR of the bond under various setups ###################################

n.paths <- 5000 # number of paths for VaR simulation
p.0 <- affine_transformation(psi.0, tau = T.mat, kappa = kappa.Q, theta = theta.Q,
                             sigma, rho.0 = r, rho.1 = rho.Q) # current bond price

## Parameter set 1: 1 year (long time horizon)
horizon <- 1 # risk management horizon
Delta.t <- 1/25 # biweekly time step for simulation of paths
n.steps <- round(horizon/Delta.t)
with.default <- TRUE # take default risk into account when computing loss distribution

## Simulate loss trajectories of gamma.Q = rho.Q * psi under P (using P-parameters)
Loss <- vector("numeric", n.paths)
set.seed(271)
for(i in 1:n.paths) {
    psi.path <- sde.sim(X0 = psi.0, theta = c(kappa.P * theta.P, kappa.P, sigma),
                        model = "CIR", T = horizon, N = n.steps)
    gammaP.path <- rho.P * psi.path
    p.horizon <- affine_transformation(psi.path[n.steps + 1], tau = (T.mat - horizon),
                                       kappa = kappa.Q, theta = theta.Q, sigma,
                                       rho.0 = r, rho.1 = rho.Q)
    if(with.default)
        p.horizon <- (1 - threshold_sim(gammaP.path, Delta.t)$indicator) * p.horizon
    Loss[i] <- -(p.horizon - p.0)
}
hist(Loss, breaks = "FD", xlim = c(-0.1, 0.95)); box()

## Parameter set 2: 10 days (short time horizon)
horizon <- 10/360 # risk management horizon
Delta.t <- 1/180 # time step of 2 days for simulation of paths
n.steps <- round(horizon/Delta.t)
with.default <- FALSE

## Simulate loss trajectories of gamma.Q = rho.Q * psi under P (using P-parameters)
Loss. <- vector("numeric", n.paths)
set.seed(271)
for(i in 1:n.paths) {
    psi.path <- sde.sim(X0 = psi.0, theta = c(kappa.P * theta.P, kappa.P, sigma),
                        model = "CIR", T = horizon, N = n.steps)
    gammaP.path <- rho.P * psi.path
    p.horizon <- affine_transformation(psi.path[n.steps + 1], tau = (T.mat - horizon),
                                       kappa = kappa.Q, theta = theta.Q, sigma,
                                       rho.0 = r, rho.1 = rho.Q)
    if(with.default)
        p.horizon <- (1 - threshold_sim(gammaP.path, Delta.t)$indicator) * p.horizon
    Loss.[i] <- -(p.horizon - p.0)
}
hist(Loss., breaks = "FD", xlim = c(-0.1, 0.95)); box()

## Risk measures
alpha <- c(0.95, 0.975, 0.99) # confidence level alpha
VaR_np(Loss, alpha = alpha, names = TRUE) # VaR
ES_np(Loss, alpha = alpha, names = TRUE) # ES
## Remark: For a long horizon (1 year), default risk is highly relevant;
##         for a short time horizon it matters only for ES_alpha for large alpha


### 2.2 Call Option pricing via Monte Carlo ####################################

## Here we can work entirely under Q

K <- p.0 # option is at the money
option.maturity <- 1 # one-year option
n.paths <- 1000
Delta.t <- 1/25 # biweekly time step for simulation of paths
n.steps <- round(option.maturity/Delta.t)
survival <- TRUE
## If survival == TRUE, the call option is treated as a survival claim. This gives
## lower variance (possible here since payoff is zero on {tau <= option.maturity})
option.payoff <- vector("numeric", n.paths)
for(i in 1:n.paths) {
    psi.path <- sde.sim(X0 = psi.0, theta = c(kappa.Q * theta.Q, kappa.Q, sigma),
                        model = "CIR", T = option.maturity, N = n.steps)
    gammaQ.path <- rho.Q * psi.path # note that we use gamma Q now
    p.horizon <- affine_transformation(psi.path[n.steps+1], tau = (T.mat - option.maturity),
                                       kappa = kappa.Q, theta = theta.Q, sigma,
                                       rho.0 = r, rho.1 = rho.Q)
    if(survival){
        Gamma.Q <- Delta.t*(sum(gammaQ.path) - 0.5*(gammaQ.path[1] + gammaQ.path[n.steps + 1])) # integral of gammaQ.path
        option.payoff[i] <- exp(-Gamma.Q) * max(p.horizon-K, 0)
    } else { # direct computation
        p.horizon <- (1- threshold_sim(gammaQ.path, Delta.t)$indicator)*p.horizon
        option.payoff[i] <- max(p.horizon-K, 0)
    }
}
exp(-r * option.maturity) * mean(option.payoff) # call price
