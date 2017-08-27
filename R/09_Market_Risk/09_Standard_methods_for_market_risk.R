## By Marius Hofert

## Estimating the risk measures VaR and ES with standard methods based on
## (BMW, Siemens) data from 2000-01-03 to 2009-12-31.


### Setup ######################################################################

library(mvtnorm) # for sampling from a multivariate t distribution
library(QRM) # for fit.mst(), fit.GPD()
library(xts) # for na.fill()
library(qrmdata) # for the data
library(qrmtools) # for the data analysis


### 1 Working with the (BMW, Siemens) data #####################################

## Load the data
data(EURSTX_const)
db <- EURSTX_const
str(db) # check the *str*ucture of 'db'

## Pick out the sub-database of stock prices we work with and fill gaps
time <- c("2000-01-03", "2009-12-31")
S <- db[paste0(time, collapse = "/"), c("BMW.DE", "SIE.DE")] # pick out data (data.frame)
str(S)
head(S) # show the beginning
tail(S) # show the end
S <- na.fill(S, fill = "extend") # fill NAs (see also NA_plot(S))
colnames(S) <- c("BMW", "SIE")

## Use scatter plots of each time series to check if anything is 'suspicious'
plot.zoo(S[,"BMW"], main = "BMW stock data",
         xlab = "Date t", ylab = expression(Stock~price~S[t]))
plot.zoo(S[,"SIE"], main = "SIEMENS stock data",
         xlab = "Date t", ylab = expression(Stock~price~S[t]))

## Compute the risk-factor changes and plot them (here: against each other)
X <- as.matrix(returns(S)) # compute log-returns and range (sign-adjustment in loss operator below)
plot(X, cex = 0.4)


### 2 Implement (and document) some auxiliary functions ########################

##' @title Compute the Loss Operator
##' @param x Matrix of risk-factor changes
##' @param weights Weights w = lambda_j * S_{t,j}
##' @return Losses
##' @author Marius Hofert
loss_operator <- function(x, weights)
    -rowSums(expm1(x) * matrix(weights, nrow = nrow(x), ncol = length(weights), byrow = TRUE))

##' @title Estimate VaR and ES
##' @param S Stock data, an (n, d)-matrix
##' @param lambda Number of shares of each stock
##' @param alpha Confidence level for VaR and ES
##' @param method A character string specifying the estimator
##' @param ... Additional arguments passed to the various methods
##' @return A list containing the estimated risk measures VaR and ES, and
##'         possibly other results (depending on the estimator)
##' @author Marius Hofert
risk_measure <- function(S, lambda, alpha,
                         method = c("Var.Cov", "hist.sim", "MC.N", "MC.t", "GPD"),
                         ...)
{
    ## Input checks and conversions
    if(!is.matrix(S)) S <- rbind(S, deparse.level = 0L) # guarantees that ncol() works
    stopifnot(0 < alpha, alpha < 1, # check whether alpha is in (0,1)
              length(lambda) == ncol(S), lambda > 0) # check length and sign of lambda
    method <- match.arg(method) # match correct method if not fully provided

    ## Ingredients required for *all* methods
    X <- returns(S) # compute risk-factor changes
    if(!length(X)) stop("'S' should have more than just one line") # check
    S. <- as.numeric(tail(S, n = 1)) # pick out last available stock prices ("today")
    w <- lambda * S. # weights w today

    ## Method switch (now consider the various methods)
    switch(method,
           "Var.Cov" = { # variance-covariance method
               ## Estimate a multivariate normal distribution
               mu.hat <- colMeans(X) # estimate the mean vector mu
               Sigma.hat  <- var(X) # estimate the covariance matrix Sigma
               L.delta.mean <- -sum(w * mu.hat) # mean of the approx. normal df of the linearized loss L^{\Delta}
               L.delta.sd <- sqrt(t(w) %*% Sigma.hat %*% w) # corresponding standard deviation
               ## Compute VaR and ES and return
               qa <- qnorm(alpha)
               list(VaR = L.delta.mean + L.delta.sd * qa,
                    ES  = L.delta.mean + L.delta.sd * dnorm(qa) / (1-alpha))
               ## => We could just return a bivariate vector here, but
               ##    for other methods, we might want to return additional
               ##    auxiliary results, and we should *always* return similar
               ##    objects (here: lists)
           },
           "hist.sim" = { # historical simulation method
               ## Using nonparametrically estimated risk measures
               L <- loss_operator(X, weights = w) # compute historical losses
               ## Nonparametrically estimate VaR and ES and return
               list(VaR = VaR_np(L, alpha),
                    ES  =  ES_np(L, alpha))
           },
           "MC.N" = { # Monte Carlo based on a fitted multivariate normal
               stopifnot(hasArg(N)) # check if the number 'N' of MC replications has been provided (via '...')
               N <- list(...)$N # pick out N from '...'
               mu.hat <- colMeans(X) # estimate the mean vector mu
               Sigma.hat  <- var(X) # estimate the covariance matrix Sigma
               X. <- rmvnorm(N, mean = mu.hat, sigma = Sigma.hat) # simulate risk-factor changes
               L <- loss_operator(X., weights = w) # compute corresponding (simulated) losses
               ## Compute VaR and ES and return
               list(VaR = VaR_np(L, alpha), # nonparametrically estimate VaR
                    ES  =  ES_np(L, alpha), # nonparametrically estimate ES
                    ## Additional quantities returned here
                    mu    = mu.hat, # fitted mean vector
                    Sigma = Sigma.hat) # fitted covariance matrix
           },
           "MC.t" = { # Monte Carlo based on a fitted multivariate t
               stopifnot(hasArg(N)) # check if the number 'N' of MC replications has been provided (via '...')
               N <- list(...)$N # pick out N from '...'
               fit <- fit.mst(X, method = "BFGS") # fit a multivariate t distribution
               X. <- rmvt(N, sigma = as.matrix(fit$Sigma), df = fit$df, delta = fit$mu) # simulate risk-factor changes
               L <- loss_operator(X., weights = w) # compute corresponding (simulated) losses
               ## Compute VaR and ES and return
               list(VaR = VaR_np(L, alpha), # nonparametrically estimate VaR
                    ES  =  ES_np(L, alpha), # nonparametrically estimate ES
                    ## Additional quantities returned here
                    mu    = fit$mu, # fitted location vector
                    sigma = fit$Sigma, # fitted dispersion matrix
                    Sigma = fit$covariance, # fitted covariance matrix
                    df    = fit$df) # fitted degrees of freedom
           },
           "GPD" = { # simulate losses from a fitted Generalized Pareto distribution (GPD); this is underlying the peaks-over-threshold method
               stopifnot(hasArg(q)) # check if the quantile-threshold 'q' has been provided
               L. <- loss_operator(X, weights = w) # historical losses
               u <- quantile(L., probs = list(...)$q, names = FALSE) # determine the threshold as the q-quantile of the historical losses
               fit <- fit.GPD(L., threshold = u) # fit a GPD to the excesses (note: fit.GPD() requires the losses)
               xi <- fit$par.ests[["xi"]] # fitted xi
               beta <- fit$par.ests[["beta"]] # fitted beta
               if(xi <= 0) stop("Risk measures only implemented for xi > 0.")
               ## Now compute semi-parametric VaR and ES estimates
               ## G_{xi,beta}(x) = 1-(1+xi*x/beta)^{-1/xi} if xi != 0
               L.. <- L.[L. > u] - u # compute the excesses over u
               Fbu <- length(L..) / length(L.) # number of excesses / number of losses = N_u / n
               VaR <- u + (beta/xi)*(((1-alpha)/Fbu)^(-xi)-1) # see McNeil, Frey, Embrechts (2015, Section 5.2.3)
               ES <- (VaR + beta-xi*u) / (1-xi) # see McNeil, Frey, Embrechts (2015, Section 5.2.3)
               if(xi >= 1) ES <- Inf # adjust to be Inf if xi >= 1 (i.e., ES < 0); see Coles (2001, p. 79)
               ## Return
               list(VaR = VaR, # parametrically estimate VaR
                    ES  = ES, # parametrically estimate ES
                    ## Additional quantities returned here
                    xi     = xi, # fitted xi
                    beta   = beta, # fitted beta
                    converged = fit$converged, # did the fitting algorithm converge?
                    u      = u, # threshold
                    excess = L..) # excesses over u
           },
           stop("Wrong 'method'"))
}


### 3 Compute VaR and ES for the standard methods ##############################

lambda <- c(1, 10) # (example) number of shares of the two stocks
alpha <- 0.99 # confidence levels for computing the risk measures
N <- 1e4 # Monte Carlo sample size

## Estimate VaR and ES with the various methods
set.seed(271) # set a seed so that all simulation results are reproducible; see ?set.seed
var.cov  <- risk_measure(S, lambda = lambda, alpha = alpha, method = "Var.Cov")
hist.sim <- risk_measure(S, lambda = lambda, alpha = alpha, method = "hist.sim")
MC.N     <- risk_measure(S, lambda = lambda, alpha = alpha, method = "MC.N", N = N)
GPD      <- risk_measure(S, lambda = lambda, alpha = alpha, method = "GPD", N = N, q = 0.9)
MC.t     <- risk_measure(S, lambda = lambda, alpha = alpha, method = "MC.t", N = N)

## Pick out VaR and ES for all methods
(rm <- rbind("MC (normal)"    = unlist(MC.N[c("VaR", "ES")]),
             "Var.-cov."      = unlist(var.cov),
             "Hist. sim."     = unlist(hist.sim),
             "GPD"            = unlist(GPD [c("VaR", "ES")]),
             "MC (Student t)" = unlist(MC.t[c("VaR", "ES")])))

### Graphical goodness-of-fit check for the GPD

## Transform the excesses with the fitted GPD distribution function.
## The resulting sample should roughly follow a standard uniform distribution.
excess <- GPD$excess # excesses over the threshold
xi.hat <- GPD$xi # estimated xi
beta.hat <- GPD$beta # estimated beta
z <- pGPD(excess, xi = xi.hat, beta = beta.hat) # should be U[0,1]
plot(z, ylab = "Fitted GPD applied to the excesses") # looks fine

## We can also consider a (more sophisticated) Q-Q plot for this task.
qq_plot(excess, FUN = function(p) qGPD(p, xi = xi.hat, beta = beta.hat),
        main = paste0("Q-Q plot for the fitted GPD(", round(xi.hat, 2),", ",
                      round(beta.hat, 2),") distribution")) # looks fine


### 4 Graphical analysis #######################################################

## Compute historical losses (for creating a histogram); see risk.measure()
S. <- as.numeric(tail(S, n = 1)) # pick out last available stock prices ("today")
w <- lambda * S. # weights w
L <- loss_operator(X, weights = w) # historical losses
summary(L) # get important statistics about the losses

## Plot a histogram of the losses including the VaR and ES estimates
hist(L, breaks = "Scott", probability = TRUE, xlim = c(0, max(L, rm)), main = "",
     xlab = substitute("Losses L > 0 from"~sd~"to"~ed~"with"~
                                         widehat(VaR)[a] <= widehat(ES)[a],
                       list(sd = time[1], ed = time[2], a = alpha)), col = "gray90") # histogram
box() # box around histogram
lty <- c(3, 2, 1, 4, 5)
lwd <- c(1.2, 1.6, 1, 1.2, 1.2)
cols <- c("maroon3", "black", "black", "royalblue3", "darkorange2")
for(k in seq_len(nrow(rm))) abline(v = rm[k,], lty = lty[k], lwd = lwd[k],
                                   col = cols[k]) # colored vertical lines indicating VaR and ES
legend("topright", bty = "n", inset = 0.02, lty = lty, lwd = lwd,
       col = cols, legend = rownames(rm),
       title = as.expression(substitute(widehat(VaR)[a] <= widehat(ES)[a], list(a = alpha)))) # legend

## Results:
## - The Monte Carlo method based on a fitted multivariate normal distribution
##   and the variance-covariance method lead to similar results
##   (they both assume multivariate normal distributed risk-factor changes),
##   both underestimating VaR and ES in comparison to the historical simulation
##   method.
## - The historical simulation method, however, implies that the loss
##   distribution is more heavy-tailed. This is captured quite well by
##   GPD-based method and (possibly too well by) the Monte Carlo method for
##   multivariate t distributed risk-factor changes
##   (degrees of freedom are MC.t$df ~= 2.4).
## - It's overall reassuring that several methods (historical simulation,
##   GPD-based EVT approach and -- with slight departure for ES --
##   MC for a Student t) lead to similar results. Somewhere in this range,
##   one can then determine an adequate amount of risk capital.
