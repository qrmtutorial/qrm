## By Marius Hofert

## Pseudo-random numbers vs quasi-random numbers


library(qrng)
library(copula)

n <- 1000 # sample size


### 1 Examples #################################################################

### 1.1 Independence copula ####################################################

## Sample
set.seed(271) # for reproducibility
U <- matrix(runif(n*2), ncol = 2) # pseudo-random numbers
U. <- ghalton(n, d = 2) # quasi-random numbers (faster alternative: sobol())

## Plot
layout(rbind(1:2))
plot(U,  xlab = expression(U[1]), ylab = expression(U[2]))
plot(U., xlab = expression(U[1]), ylab = expression(U[2]))


### 1.2 t copula ###############################################################

## Define the copula
nu <- 3 # degrees of freedom
th <- iTau(tCopula(df = nu), tau = 0.5) # correlation parameter
cop.t3 <- tCopula(param = th, df = nu) # t copula

## Sample
U.t <-  cCopula(U,  copula = cop.t3, inverse = TRUE)
U.t. <- cCopula(U., copula = cop.t3, inverse = TRUE)

## Plot
layout(rbind(1:2))
plot(U.t,  xlab = expression(U[1]), ylab = expression(U[2]))
plot(U.t., xlab = expression(U[1]), ylab = expression(U[2]))
layout(1) # reset layout

## 3d plot
U.3d. <- ghalton(n, d = 3)
cop <- tCopula(param = th, dim = 3, df = nu)
U.t.3d. <- cCopula(U.3d., copula = cop, inverse = TRUE)
cloud2(U.t.3d., xlab = expression(U[1]), ylab = expression(U[2]),
       zlab = expression(U[3])) # not much visible
pairs2(U.t.3d., labels.null.lab = "U")
## Structure doesn't *look* very convincing, but it is a low-discrepancy sequence


### 1.3 Clayton copula #########################################################

## Define the copula
th <- iTau(claytonCopula(), tau = 0.5)
cop.C <- claytonCopula(th)

## Sample
U.C <-  cCopula(U,  copula = cop.C, inverse = TRUE)
U.C. <- cCopula(U., copula = cop.C, inverse = TRUE)

## Plot
layout(rbind(1:2))
plot(U.C,  xlab = expression(U[1]), ylab = expression(U[2]))
plot(U.C., xlab = expression(U[1]), ylab = expression(U[2]))
layout(1) # reset layout


### 2 Why considering quasi-random numbers? ####################################

## Computing P(U_1 > u_1,..., U_d > u_d) by MC based on PRNG, QRNG
tail_prob <- function(n, copula, u) # sample size, copula, lower-left endpoint
{
    d <- length(u)
    stopifnot(n >= 1, inherits(copula, "Copula"), 0 < u, u < 1,
              d == dim(copula))
    umat <- rep(u, each = n)

    ## Pseudo-random numbers
    clock <- proc.time() # start watch
    U <- rCopula(n, copula = copula) # unfair comparison but *still* slower than Sobol' + cCopula()
    prob.PRNG <- mean(rowSums(U > umat) == d) # = sum(<rows with all entries TRUE>) / n
    usr.time.PRNG <- 1000 * (proc.time() - clock)[[1]] # time in ms

    ## (Randomized) quasi-random numbers
    clock <- proc.time() # start watch
    U. <- cCopula(sobol(n, d = d, randomize = TRUE), # needed for unbiasedness
                  copula = copula, inverse = TRUE)
    prob.QRNG <- mean(rowSums(U. > umat) == d) # = sum(<rows with all entries TRUE>) / n
    usr.time.QRNG <- 1000 * (proc.time() - clock)[[1]] # time in ms

    ## Return
    list(PRNG = c(prob = prob.PRNG, rt = usr.time.PRNG),
         QRNG = c(prob = prob.QRNG, rt = usr.time.QRNG))
}

## Simulate the probabilities of falling in (u_1, 1] x ... x (u_d, 1]
N <- 200 # number of replications
n <- 2000 # sample size
d <- 2 # dimension
u <- rep(0.99, d) # lower-left endpoint of the considered cube
nu <- 3 # degrees of freedom
tau <- 0.5 # Kendall's tau
rho <- iTau(tCopula(df = nu), tau = 0.5) # correlation parameter
cop <- tCopula(param = rho, dim = d, df = nu) # t copula

## Run
set.seed(271) # for reproducibility
system.time(res <- lapply(1:N, function(N.) tail_prob(n, copula = cop, u = u)))

## Grab out the results
prob.PRNG <- sapply(res, function(x) x[["PRNG"]][["prob"]])
prob.QRNG <- sapply(res, function(x) x[["QRNG"]][["prob"]])
rt.PRNG   <- sapply(res, function(x) x[["PRNG"]][["rt"]])
rt.QRNG   <- sapply(res, function(x) x[["QRNG"]][["rt"]])

## Boxplot of computed exceedance probabilities
boxplot(list(PRNG = prob.PRNG, QRNG = prob.QRNG),
        main = substitute("Simulated"~
                          P(bold(U) > bold(u))~~"for a"~t[nu.]~"copula",
                          list(nu. = nu)))
mtext(sprintf("N = %d replications with n = %d and d = %d", N, n, d),
      side = 4, line = 1, adj = 0, las = 0)
## QRNGs provide a smaller variance

## Variance reduction factor and % improvement
(vrf <- var(prob.PRNG)/var(prob.QRNG)) # estimated variance reduction factor

## Boxplot of measured run times (in milliseconds)
boxplot(list(PRNG = rt.PRNG, QRNG = rt.QRNG), outline = FALSE,
        main = substitute("Run times for a" ~ t[nu.]~"copula",
                          list(nu. = nu)))
mtext(sprintf("N = %d replications with n = %d and d = %d", N, n, d),
      side = 4, line = 1, adj = 0, las = 0)
## Clearly (the way we sampled here!), QRNG is slower here in this comparison

## Run-time reduction factor
(rrf <- mean(rt.PRNG)/mean(rt.QRNG)) # estimated run-time factor PRNG w.r.t. QRNG
## < 1 (PRNG faster) here because rCopula() is *way* faster than cCopula()

## Effort comparison (variance per unit of time)
vrf * rrf # = (Var('PRNG') * time('PRNG')) / (Var('QRNG') * time('QRNG'))
## Still better to use QRNG! (if > 1, QRNG is preferable)

## Remark:
## - QRNGs can estimate tail probabilities with smaller variance than PRNGs
##   (=> need less random variates to obtain the same precision => good for
##    memory/storage-intensive methods).
## - QRNGs can be faster than PRNGs (but it depends: Sobol': yes; generalized
##   Halton: no); this may depend on the dimension, too.
## - If more efficient sampling methods for copula-QRNGs are used (so similar
##   to rCopula() instead of cCopula()), run time for QRNGs is significantly
##   reduced; see Cambou, Hofert, Lemieux ("Quasi-random numbers for copula models").
## - Even if slower, it can still be advantageous to use QRNGs instead of PRNGs
##   (because of memory/storage limitations).

## Just comparing a PRNG and QRNGs
n. <- 2e7 # 20 Mio
set.seed(271)
system.time(runif(n.)) # PRNG
system.time(ghalton(n., d = 1)) # QRNG
system.time(sobol(n.,   d = 1, randomize = TRUE)) # Faster QRNG
## Run time also depends on the *type* of QRNG

if(FALSE) {
    ## Profiling: See where run time is spent
    Rprof(profiling <- tempfile(), line.profiling = TRUE) # enable profiling
    res <- lapply(1:N, function(N.) tail_prob(n, copula = cop, u = u))
    Rprof(NULL) # disable profiling
    (profile <- summaryRprof(profiling, lines = "both")) # get a summary
    ## => by.total => most of the run time is spent in cCopula()
}
