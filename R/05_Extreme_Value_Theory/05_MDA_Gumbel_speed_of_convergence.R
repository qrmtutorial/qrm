## By Marius Hofert

## Investigating the speed of convergence for various distributions in the
## Gumbel MDA. Parametrically (based on the density) and nonparametrically
## (based on kernel density estimates).


### Setup ######################################################################

library(qrmtools)

## Block sizes and number of blocks
n <- 10^c(1,3,5) # block sizes

## Distributional parameters
pmExp <- 2 # Exp(lambda) parameter lambda
pmG  <- c(2, 2) # Gamma(alpha, beta) parameters alpha, beta (rate)
pmN  <- c(1, 2) # N(mu, sigma^2) parameter mu, sigma
pmLN <- c(1, 2) # LN(mu, sigma^2) parameter mu, sigma

## Evaluation points
x <- seq(-3, 10, length.out = 257)


### 1 Auxiliary functions ######################################################

## Distribution functions
pFun <- list(function(q, ...) pexp(q, rate = pmExp, ...),
             function(q, ...) pgamma(q, shape = pmG[1], rate = pmG[2], ...),
             function(q, ...) pnorm(q, mean = pmN[1], sd = pmN[2], ...),
             function(q, ...) plnorm(q, meanlog = pmLN[1], sdlog = pmLN[2], ...))
ndf <- length(pFun)
## Strings for plotting
dfs <- as.expression(c(paste0("Exp(",pmExp,")"),
                       substitute(Gamma(a,b), list(a = pmG[1], b = pmG[2])),
                       paste0("N(",pmN[1],", ",pmN[2],")"),
                       paste0("LN(",pmLN[1],", ",pmLN[2],")")))

## Corresponding densities
dFun <- list(function(x, ...) dexp(x, rate = pmExp, ...),
             function(x, ...) dgamma(x, shape = pmG[1], rate = pmG[2], ...),
             function(x, ...) dnorm(x, mean = pmN[1], sd = pmN[2], ...),
             function(x, ...) dlnorm(x, meanlog = pmLN[1], sdlog = pmLN[2], ...))

## Corresponding quantile functions
qFun <- list(function(p, ...) qexp(p, rate = pmExp, ...),
             function(p, ...) qgamma(p, shape = pmG[1], rate = pmG[2], ...),
             function(p, ...) qnorm(p, mean = pmN[1], sd = pmN[2], ...),
             function(p, ...) qlnorm(p, meanlog = pmLN[1], sdlog = pmLN[2], ...))

## Normalizing sequences (c_n) and (d_n)
## see Embrechts, Klueppelberg, Mikosch (1997, pp. 155)
d_n <- function(n, k)
    switch(k,
    { log(n)/pmExp }, # for exponential
    { (log(n) + (pmG[1] - 1) * log(log(n)) - lgamma(pmG[1])) / pmG[2] }, # for gamma
    { pmN[1] + pmN[2] * (sqrt(2*log(n)) - (log(4*pi) + log(log(n))) / (2*sqrt(2*log(n)))) }, # for normal (can be derived based on p. 156)
    { exp(pmLN[1] + pmLN[2] * (sqrt(2*log(n)) - (log(4*pi) + log(log(n))) / (2*sqrt(2*log(n))))) }, # for lognormal
    stop("Wrong case"))
c_n <- function(n, k)
    switch(k,
    { 1/pmExp }, # for exponential
    { 1/pmG[2] }, # for gamma
    { pmN[2] / sqrt(2*log(n)) }, # for normal (can be derived based on p. 156)
    { pmLN[2] * d_n(n, k = k) / sqrt(2*log(n)) }, # for lognormal
    stop("Wrong case"))


### 2 Theoretical density of the location-scale transformed block maxima #######

## Block maxima method (for different distributions and sample sizes)
res.p <- vector("list", length = ndf)
for(k in seq_len(ndf)) { # distribution
    res.p[[k]] <- lapply(seq_along(n), function(i) { # sample size
        c.n <- c_n(n[i], k = k)
        d.n <- d_n(n[i], k = k)
        ## Evaluate the density of the location-scaled transformed block maxima
        ## Note: 1) F_{(M_n - d_n)/c_n} (x) = F^n(d_n + c_n x) and thus
        ##          f_{(M_n - d_n)/c_n} (x) = n c_n f(d_n + c_n x) F^{n-1}(d_n + c_n x)
        ##       2) We use a numerically more robust version of
        ##          n[i] * c.n * dFun[[k]](d.n + c.n * x) * (pFun[[k]](d.n + c.n * x))^(n[i]-1)
        ly <- log(n[i]) + log(c.n) + dFun[[k]](d.n + c.n * x, log = TRUE) +
                                             (n[i]-1) * pFun[[k]](d.n + c.n * x, log.p = TRUE)
        cbind(x = x, y = exp(ly))
    })
}

## Plot
dLambda <- exp(-exp(-x)-x) # density of Lambda(x) = exp(-exp(-x)) at x
xlim <- range(x)
ylim <- range(dLambda, unlist(lapply(res.p, function(x) lapply(x, function(x.) x.[,"y"]))))
opar <- par(mar = c(5, 4, 4, 2) + 0.1 - c(1, 1, 3, 0)) # bottom, left, top, right
layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))
k <- 1
for(i in 1:2) {
    for(j in 1:2) {
        r <- res.p[[k]]
        plot(x, dLambda, type = "l", xlim = xlim, ylim = ylim,
             xlab = "x", ylab = "") # limiting density
        lines(r[[1]][,"x"], r[[1]][,"y"], col = "royalblue3") # n = n[1]
        lines(r[[2]][,"x"], r[[2]][,"y"], col = "maroon3") # n = n[2]
        lines(r[[3]][,"x"], r[[3]][,"y"], col = "darkorange2") # n = n[3]
        legend("topright", bty = "n", lty = rep(1, 4),
               col = c("black", "royalblue3", "maroon3", "darkorange2"),
               legend = as.expression(c("Gumbel",
                                        lapply(1:3, function(i) substitute(n == n., list(n. = n[i]))))))
        mtext(dfs[k], side = 4, adj = 0, line = 0.75) # label
        k <- k+1
    }
}
layout(1)
par(opar)


### 3 Nonparametrically (based on kernel density estimates) ####################

## Generate data (here: with inversion method in order to recycle)
m <- 100 # number of blocks (each of size n/m) for each n; results for m = 200 don't seem much better
ns <- n * m
ns.max <- max(ns) # maximal sample size (for recycling)
set.seed(271)
U <- runif(ns.max)
system.time(dat <- lapply(qFun, function(qF) qF(U))) # use same U's; about 1min

## Block maxima method (for different distributions and sample sizes)
res.np <- vector("list", length = ndf)
for(k in seq_len(ndf)) { # distribution
    res.np[[k]] <- lapply(seq_along(n), function(i) { # sample size
        c.n <- c_n(n[i], k = k)
        d.n <- d_n(n[i], k = k)
        ## Build location-scaled transformed block maxima (transformed with
        ## known normalizing sequences c_n, d_n)
        dat. <- head(dat[[k]], n = ns[i]) # grab out the data we work with
        dat.blocks <- split(dat., f = rep(1:m, each = n[i])) # split it in blocks
        M <- sapply(dat.blocks, function(x) (max(x) - d.n) / c.n)
        ## Fit density to the location-scaled transformed block maxima
        fit.dens <- density(M) # kernel density estimate
        cbind(x = fit.dens$x, y = fit.dens$y)
    })
}

## Plot
ylim <- range(dLambda, unlist(lapply(res.np, function(x) lapply(x, function(x.) x.[,"y"]))))
opar <- par(mar = c(5, 4, 4, 2) + 0.1 - c(1, 1, 3, 0)) # bottom, left, top, right
layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))
k <- 1
for(i in 1:2) {
    for(j in 1:2) {
        r <- res.np[[k]]
        plot(x, dLambda, type = "l", xlim = xlim, ylim = ylim,
             xlab = "x", ylab = "") # limiting density
        lines(r[[1]][,"x"], r[[1]][,"y"], col = "royalblue3") # n = n[1]
        lines(r[[2]][,"x"], r[[2]][,"y"], col = "maroon3") # n = n[2]
        lines(r[[3]][,"x"], r[[3]][,"y"], col = "darkorange2") # n = n[3]
        legend("topright", bty = "n", lty = rep(1, 4),
               col = c("black", "royalblue3", "maroon3", "darkorange2"),
               legend = as.expression(c("Gumbel",
                                        lapply(1:3, function(i) substitute(n == n., list(n. = n[i]))))))
        mtext(dfs[k], side = 4, adj = 0, line = 0.75) # label
        k <- k+1
    }
}
layout(1)
par(opar)
