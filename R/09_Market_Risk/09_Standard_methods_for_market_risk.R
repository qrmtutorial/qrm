## By Marius Hofert

## Estimating the risk measures VaR and ES with standard methods based on
## (BMW, Siemens) data from 1985-01-02 to 1994-12-30.


### Setup ######################################################################

library(mvtnorm)
library(QRM) # for fit.mst(), fit.GPD()


### 1 Working with the (BMW, Siemens) data #####################################

## Download the data; if this fails, download it manually
file <- "DAX.RAW.txt"
if(!file.exists(file))
    download.file(paste0("http://www.ma.hw.ac.uk/~mcneil/ftp/", file),
                  destfile=file, method="auto") # download the data
db <- read.table(file, header=TRUE, row.names="Positions", # read the *d*ata *b*ase
                 encoding="UTF-8") # for 'didactical' reasons
str(db) # check the *str*ucture of 'db'
str(rn <- rownames(db)) # row names contain dates

## Convert dates to *proper* date objects
date <- as.Date(rn, format="%m/%d/%Y") # convert the dates to *proper* date objects
str(date)
rownames(db) <- date # use the newly formatted dates as row names

## Pick out the sub-database of stock prices we work with here
start.date <- "1985-01-02"
end.date   <- "1994-12-30"
S <- db[start.date <= date & date <= end.date, c("BMW", "SIEMENS")] # pick out data (data.frame)
str(S)
head(S) # show the beginning
S[1:6,] # the same
tail(S) # show the end

## Use scatter plots of each time series to check if anything is 'suspicious'
rns <- rownames(S) # row names of S (character)
## BMW
plot(as.Date(rns), S[,"BMW"], type="l",
     main="BMW stock data", # title
     xlab="Date t", # x-axis label
     ylab=expression(Stock~price~S[t])) # y-axis label
## Siemens
plot(as.Date(rns), S[,"SIEMENS"], type="l",
     main="Siemens stock data", # title
     xlab="Date t", # x-axis label
     ylab=expression(Stock~price~S[t])) # y-axis label

## Compute the risk-factor changes and plot them (here: against each other)
ran <- range(X <- apply(log(S), 2, diff)) # compute log-returns and range
plot(X, xlim=ran, ylim=ran, main="Risk-factor changes", cex=0.2)


### 2 Implement (and document) some auxiliary functions ########################

##' @title Compute the Loss Operator
##' @param x Matrix of risk-factor changes
##' @param weights Weights w = lambda_j * S_{t,j}
##' @return Losses
##' @author Marius Hofert
loss_operator <- function(x, weights)
    -rowSums(expm1(x) * matrix(weights, nrow=nrow(x), ncol=length(weights), byrow=TRUE))

##' @title Non-parametric VaR estimator
##' @param L Losses
##' @param alpha Confidence level
##' @return Non-parametric estimate of VaR at level alpha
##' @author Marius Hofert
VaR_hat <- function(L, alpha) quantile(L, probs=alpha, names=FALSE)

##' @title Non-parametric ES estimator
##' @param L Losses
##' @param alpha Confidence level
##' @return Non-parametric estimate of ES at level alpha
##' @author Marius Hofert
ES_hat  <- function(L, alpha) mean(L[L > VaR_hat(L, alpha=alpha)])

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
    if(!is.matrix(S)) S <- rbind(S, deparse.level=0L) # guarantees that ncol() works
    stopifnot(0 < alpha, alpha < 1, # check whether alpha is in (0,1)
              length(lambda) == ncol(S), lambda > 0) # check length and sign of lambda
    method <- match.arg(method) # match correct method if not fully provided

    ## Ingredients required for *all* methods
    X <- apply(log(S), 2, diff) # compute risk-factor changes
    if(!length(X)) stop("'S' should have more than just one line") # check
    S. <- as.numeric(tail(S, n=1)) # pick out last available stock prices ("today")
    w <- lambda * S. # weights w

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
               L <- loss_operator(X, weights=w) # compute historical losses
               ## Nonparametrically estimate VaR and ES and return
               list(VaR = VaR_hat(L, alpha),
                    ES =   ES_hat(L, alpha))
           },
           "MC.N" = { # Monte Carlo based on a fitted multivariate normal
               stopifnot(hasArg(N)) # check if the number 'N' of MC replications has been provided (via '...')
               N <- list(...)$N # pick out N from '...'
               mu.hat <- colMeans(X) # estimate the mean vector mu
               Sigma.hat  <- var(X) # estimate the covariance matrix Sigma
               X. <- rmvnorm(N, mean=mu.hat, sigma=Sigma.hat) # simulate risk-factor changes
               L <- loss_operator(X., weights=w) # compute corresponding (simulated) losses
               ## Compute VaR and ES and return
               list(VaR = VaR_hat(L, alpha), # nonparametrically estimate VaR
                    ES  =  ES_hat(L, alpha), # nonparametrically estimate ES
                    ## Additional quantities returned here
                    mu    = mu.hat, # fitted mean vector
                    Sigma = Sigma.hat) # fitted covariance matrix
           },
           "MC.t" = { # Monte Carlo based on a fitted multivariate t
               stopifnot(hasArg(N)) # check if the number 'N' of MC replications has been provided (via '...')
               N <- list(...)$N # pick out N from '...'
               fit <- fit.mst(X, method = "BFGS") # fit a multivariate t distribution
               X. <- rmvt(N, sigma=as.matrix(fit$Sigma), df=fit$df, delta=fit$mu) # simulate risk-factor changes
               L <- loss_operator(X., weights=w) # compute corresponding (simulated) losses
               ## Compute VaR and ES and return
               list(VaR = VaR_hat(L, alpha), # nonparametrically estimate VaR
                    ES =   ES_hat(L, alpha), # nonparametrically estimate ES
                    ## Additional quantities returned here
                    mu    = fit$mu, # fitted location vector
                    sigma = fit$Sigma, # fitted dispersion matrix
                    Sigma = fit$covariance, # fitted covariance matrix
                    df    = fit$df) # fitted degrees of freedom
           },
           "GPD" = { # simulate losses from a fitted Generalized Pareto distribution (GPD); this is underlying the peaks-over-threshold method
               stopifnot(hasArg(q)) # check if the quantile-threshold 'q' has been provided
               L. <- loss_operator(X, weights=w) # historical losses
               u <- quantile(L., probs=list(...)$q, names=FALSE) # determine the threshold as the q-quantile of the historical losses
               fit <- fit.GPD(L., threshold=u) # fit a GPD to the excesses
               xi <- fit$par.ests[["xi"]] # fitted xi
               beta <- fit$par.ests[["beta"]] # fitted beta
               if(xi <= 0) stop("Risk measures only implemented for xi > 0.")
               ## Now compute semi-parametric VaR and ES estimates
               ## G_{xi,beta}(x) = 1-(1+xi*x/beta)^{-1/xi} if xi != 0
               L.. <- L.[L. > u] - u # compute the excesses over u
               Fbu <- length(L..) / length(L.) # = N_u/n
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
var.cov  <- risk_measure(S, lambda=lambda, alpha=alpha, method="Var.Cov")
hist.sim <- risk_measure(S, lambda=lambda, alpha=alpha, method="hist.sim")
MC.N     <- risk_measure(S, lambda=lambda, alpha=alpha, method="MC.N", N=N)
GPD      <- risk_measure(S, lambda=lambda, alpha=alpha, method="GPD", N=N, q=0.9)
MC.t     <- risk_measure(S, lambda=lambda, alpha=alpha, method="MC.t", N=N)

## Pick out VaR and ES for all methods
(rm <- rbind("Var.-cov."   = unlist(var.cov),
             "MC (normal)" = unlist(MC.N[c("VaR", "ES")]),
             "Hist. sim."  = unlist(hist.sim),
             "GPD"         = unlist(GPD [c("VaR", "ES")]),
             "MC (Student t)"      = unlist(MC.t[c("VaR", "ES")])))

### Graphical goodness-of-fit check for the GPD

## Transform the excesses with the fitted GPD distribution function.
## The resulting sample should roughly follow a standard uniform distribution.
excess <- GPD$excess # excesses over the threshold
xi.hat <- GPD$xi # estimated xi
beta.hat <- GPD$beta # estimated beta
z <- pGPD(excess, xi=xi.hat, beta=beta.hat) # should be U[0,1]
plot(z, ylab="Fitted GPD applied to the excesses") # looks fine

## We can also consider a (more sophisticated) Q-Q plot for this task.
excess. <- sort(excess) # sorted data
qF <- function(p) qGPD(p, xi=xi.hat, beta=beta.hat)
qF. <- qF(ppoints(length(excess.))) # theoretical quantiles
plot(qF., excess., xlab="Theoretical quantiles",
     ylab="Sample quantiles", main=paste0("Q-Q plot for the fitted GPD(",
     round(xi.hat, 2),", ",round(beta.hat, 2),") distribution"))
qqline(y=excess., distribution=qF) # looks fine


### 4 Graphical analysis #######################################################

## Compute historical losses (for creating a histogram); see risk.measure()
S. <- as.numeric(tail(S, n=1)) # pick out last available stock prices ("today")
w <- lambda * S. # weights w
L <- loss_operator(X, weights=w) # historical losses
summary(L) # get important statistics about the losses

## Plot a histogram of the losses including the VaR and ES estimates
hist(L, breaks="Scott", probability=TRUE, xlim=c(0, max(L, rm)), main="",
     xlab=substitute("Losses L > 0 from"~sd~"to"~ed~"with"~
                     widehat(VaR)[a]<=widehat(ES)[a],
                     list(sd=start.date, ed=end.date, a=alpha)), col="gray90") # histogram
box() # box around histogram
lty <- c(3,1,2,4,5)
lwd <- c(2.2,1,1.4,1,1.2)
for(k in seq_len(nrow(rm))) abline(v=rm[k,], lty=lty[k], lwd=lwd[k]) # colored vertical lines indicating VaR and ES
legend("topright", bty="n", inset=0.02, lty=lty, lwd=lwd, legend=rownames(rm),
       title=as.expression(substitute(widehat(VaR)[a]<=widehat(ES)[a], list(a=alpha)))) # legend

## Results:
## 1) The variance-covariance method and the Monte Carlo method based
##    on a fitted multivariate normal distribution lead to similar results
##    (they both assume multivariate normal distributed risk-factor changes).
## 2) The historical simulation method, however, implies that the loss
##    distribution is more heavy-tailed. This is captured quite well by the
##    risk measure estimates based on the Monte Carlo method for multivariate t
##    distributed risk-factor changes (degrees of freedom are MC.t$df ~= 3.02);
## 3) Also the GPD-based EVT approach leads results comparable to the Monte Carlo
##    method based on a fitted t distribution and the historical simulation method.
## It's overall reassuring that several methods (historical simulation,
## MC (Student t) and GPD-based EVT approach) lead to similar results. Somewhere
## in this range, upper management can then determine an adequate amount of
## risk capital.
