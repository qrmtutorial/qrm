## By Marius Hofert

## Visually explaining Sklar's Theorem and the invariance principle


library(mvtnorm)

set.seed(271)

## Sample from a t copula
n <- 1000 # sample size
d <- 2 # dimension
rho <- 0.7 # off-diagonal entry in the correlation matrix P
P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
diag(P) <- 1
nu <- 3.5 # degrees of freedom
set.seed(271)
X <- rmvt(n, sigma = P, df = nu) # n multiv. t observations
U <- pt(X, df = nu) # n ind. realizations from the corresponding t copula
Y <- qexp(U) # transform U (t copula) to Exp(1) margins

## Plot setup
ind <- c(A = 119, B = 516, C = 53) # use 'plot(X); identify(X)' to identify these points
cols <- c("royalblue3", "maroon3", "darkorange2")
par(pty = "s")

## Scatter plot of a bivariate t distribution
plot(X, xlab = quote(X[1]), ylab = quote(X[2]))
points(X[ind["A"],1], X[ind["A"],2], pch = 21, bg = cols[1], col = cols[1])
points(X[ind["B"],1], X[ind["B"],2], pch = 21, bg = cols[2], col = cols[2])
points(X[ind["C"],1], X[ind["C"],2], pch = 21, bg = cols[3], col = cols[3])
text(X[ind["A"],1], X[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1])
text(X[ind["B"],1], X[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2])
text(X[ind["C"],1], X[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3])

## Scatter plot of the corresponding t copula (1st part of Sklar's Theorem)
plot(U, xlab = quote(U[1]), ylab = quote(U[2]))
points(U[ind["A"],1], U[ind["A"],2], pch = 21, bg = cols[1], col = cols[1])
points(U[ind["B"],1], U[ind["B"],2], pch = 21, bg = cols[2], col = cols[2])
points(U[ind["C"],1], U[ind["C"],2], pch = 21, bg = cols[3], col = cols[3])
text(U[ind["A"],1], U[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1])
text(U[ind["B"],1], U[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2])
text(U[ind["C"],1], U[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3])

## Scatter plot of the meta-t distribution with Exp(1) margins (2nd part of Sklar's Theorem)
plot(Y, xlab = quote(Y[1]), ylab = quote(Y[2]))
points(Y[ind["A"],1], Y[ind["A"],2], pch = 21, bg = cols[1], col = cols[1])
points(Y[ind["B"],1], Y[ind["B"],2], pch = 21, bg = cols[2], col = cols[2])
points(Y[ind["C"],1], Y[ind["C"],2], pch = 21, bg = cols[3], col = cols[3])
text(Y[ind["A"],1], Y[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1])
text(Y[ind["B"],1], Y[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2])
text(Y[ind["C"],1], Y[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3])
## We see that the *relative* locations of all points remains the same.
## We thus only change the marginal distributions, but not the dependence
## (by applying the probability and quantile transformations). This also
## visually confirms the invariance principle.