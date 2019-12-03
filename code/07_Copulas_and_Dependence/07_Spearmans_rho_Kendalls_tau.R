## By Marius Hofert

## Various experiments with Spearman's rho and Kendall's tau


### Setup ######################################################################

library(copula)


### 1 How close do we get to the maximum log-likelihood for estimating P of a
###   Gaussian copula if Spearman's rho is used to estimate P? #################

## Generate pseudo-observations
set.seed(271)
nc <- normalCopula(0.6, dim = 3) # define a Gauss copula
U <- pobs(rCopula(1000, copula = nc)) # sample pseudo-observations from it

## Maximum pseudo-likelihood estimation of P
fit.N <- fitCopula(normalCopula(, dim = 3, dispstr = "un"), data = U)
fit.N@loglik # log-likelihood

## Estimation of P via Spearman's rho
P. <- cor(U, method = "spearman")
sum(dCopula(U, copula = normalCopula(P2p(P.), dim = 3, dispstr = "un"), log = TRUE)) # log-likelihood


### 2 How close do we get to the maximum log-likelihood for estimating P of a
###   t copula if Kendall's tau is used to estimate P? #########################

## Generate pseudo-observations
set.seed(271)
tc <- tCopula(0.6, dim = 3, df = 10) # define a t_10 copula
U <- pobs(rCopula(1000, copula = tc)) # sample pseudo-observations from it

## Maximum pseudo-likelihood estimation of P and the d.o.f.
fit.t <- fitCopula(tCopula(, dim = 3, dispstr = "un"), data = U)
fit.t@loglik # log-likelihood

## Estimation of P via Kendall's tau (and of the d.o.f. via the likelihood)
fit.t. <- fitCopula(tCopula(, dim = 3, dispstr = "un"), data = U, method = "itau.mpl")
fit.t.@loglik # log-likelihood


### 3 How well does a Gauss copula fit t copula data? ##########################

fit.N <- fitCopula(normalCopula(, dim = 3, dispstr = "un"), data = U)
fit.N@loglik


### 4 Shall we estimate correlations with the sample correlation or by inverting
###   Kendall's tau? ###########################################################

## Generate data
df <- 3
tc <- tCopula(0.5, df = df) # define a t_3 copula
set.seed(271)
U <- lapply(1:3000, function(i) rCopula(90, copula = tc)) # (90, 2, 3000)-array
X <- lapply(U, function(x) qt(x, df = df)) # map to multivariate t

## Compute the correlation coefficient and Kendall's tau
cor.pearson <- sapply(X, function(x) cor(x)[2,1])
cor.tau <- iTau(tc, tau = sapply(X, function(x) cor(x, method = "kendall")[2,1]))

## Plot
ran <- range(cor.pearson, cor.tau)
plot(cor.pearson, ylim = ran, xlab = "Sample",
     ylab = expression("Estimates of"~rho))
points(cor.tau, col = "royalblue3")
legend("bottomleft", pch = c(1,1), bty = "n", col = c("black", "royalblue"),
       legend = c("Via sample correlation coefficient", "Via inverting Kendall's tau"))

## Variance
var(cor.pearson)
var(cor.tau) # smaller
var(sapply(U, function(x) cor(x)[2,1])) # equally small (= Spearman's rho; based on U)


### 5 Approximation of Spearman’s rho for t copulas by Pearson's rho ###########

## Note: If it exists, Pearson's rho equals the correlation parameter of the
##       t copula

## Random angle for computing Spearman's rho for t copulas
Theta <- function(nu, n.MC = 1e5)
{
    W1 <- nu/rchisq(n.MC, df=nu)
    W2 <- nu/rchisq(n.MC, df=nu)
    W3 <- nu/rchisq(n.MC, df=nu)
    2 * asin(sqrt(W1 * W2 / ((W1 + W3) * (W2 + W3))))
}

## Spearman's rho for t copulas
## (see MFE (2015, Prop. 7.44 and Eq. (7.42)))
rho_t <- function(rho, nu, n.MC = 1e5)
{
    stopifnot(-1 <= rho, rho <= 1, nu > 0, n.MC > 0)
    (6 / pi) * colMeans(asin(outer(sin(Theta(nu, n.MC = n.MC) / 2), rho, FUN = "*"))) # length(rho)
}

## Generate data
set.seed(271)
nu <- 10^seq(-1, 1, length.out = 3)
rho <- seq(-1, 1, by = 0.01)
rho.t <- sapply(nu, function(nu.) rho_t(rho, nu = nu.)) # (length(rho), length(nu)-matrix

## Plot of Spearman's rho vs Pearson's rho for t copulas
cols <- c("black", "maroon3", "darkorange2", "royalblue3")
plot(rho, rho.t[,1], type = "l", col = cols[2],
     xlab = expression("Pearson's"~~rho), ylab = expression("Spearman's"~~rho[S]))
lines(rho, rho.t[,2], col = cols[3])
lines(rho, rho.t[,3], col = cols[4])
abline(0, 1)
legend("bottomright", bty = "n", lty = rep(1, 4), col = cols,
       legend = c(expression(rho[S] == rho),
                  substitute(nu==nu., list(nu.=nu[1])),
                  substitute(nu==nu., list(nu.=nu[2])),
                  substitute(nu==nu., list(nu.=nu[3]))))
