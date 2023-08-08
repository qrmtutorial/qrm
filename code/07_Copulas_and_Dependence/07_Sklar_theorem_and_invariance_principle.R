## By Marius Hofert

## Visually explaining Sklar's Theorem and the invariance principle

## Setup
library(nvmix) # for rNorm() and rStudent()
set.seed(271)

n <- 1000 # sample size
d <- 2 # dimension
rho <- 0.7 # off-diagonal entry in the correlation matrix P
P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
diag(P) <- 1
nu <- 4 # degrees of freedom for Student's t distribution


### 1 Bivariate normal, bivariate normal copula, and with Exp(1) margins #######

set.seed(271) # for reproducibility
X <- rNorm(n, scale = P) # n multivariate normal observations
U <- pnorm(X) # n independent realizations from the corresponding normal copula
Y <- qexp(U) # transform U (normal copula) to Exp(1) margins

## Plot setup
ind <- c(A = 719, B = 490, C = 83) # use 'plot(X); identify(X)' to identify these points
cols <- c("royalblue3", "maroon3", "darkorange2")
opar <- par(pty = "s", pch = 20) # set plot parameters

## Scatter plot of a bivariate normal
plot(X, xlab = expression(X[1]), ylab = expression(X[2]))
points(X[ind["A"],1], X[ind["A"],2], bg = cols[1], col = cols[1])
points(X[ind["B"],1], X[ind["B"],2], bg = cols[2], col = cols[2])
points(X[ind["C"],1], X[ind["C"],2], bg = cols[3], col = cols[3])
text(X[ind["A"],1], X[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1], font = 2)
text(X[ind["B"],1], X[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2], font = 2)
text(X[ind["C"],1], X[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3], font = 2)

## Scatter plot of the corresponding normal copula sample (1st part of Sklar's Theorem)
plot(U, xlab = expression(U[1]==F(X[1])), ylab = expression(U[2]==F(X[2])))
points(U[ind["A"],1], U[ind["A"],2], bg = cols[1], col = cols[1])
points(U[ind["B"],1], U[ind["B"],2], bg = cols[2], col = cols[2])
points(U[ind["C"],1], U[ind["C"],2], bg = cols[3], col = cols[3])
text(U[ind["A"],1], U[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1], font = 2)
text(U[ind["B"],1], U[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2], font = 2)
text(U[ind["C"],1], U[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3], font = 2)

## Scatter plot of the meta-normal distribution with Exp(1) margins (2nd part of Sklar's Theorem)
plot(Y, xlab = expression(Y[1]=={F[1]^{-1}}(U[1])), ylab = expression(Y[2]=={F[2]^{-1}}(U[2])))
points(Y[ind["A"],1], Y[ind["A"],2], bg = cols[1], col = cols[1])
points(Y[ind["B"],1], Y[ind["B"],2], bg = cols[2], col = cols[2])
points(Y[ind["C"],1], Y[ind["C"],2], bg = cols[3], col = cols[3])
text(Y[ind["A"],1], Y[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1], font = 2)
text(Y[ind["B"],1], Y[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2], font = 2)
text(Y[ind["C"],1], Y[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3], font = 2)
## We see that the *relative* locations of all points remains the same.
## By applying the probability and quantile transformations, we thus only
## change the marginal distributions, not the dependence structure. This also
## visually confirms the invariance principle.


### 2 Bivariate t_4, bivariate t_4 copula, and with Exp(1) margins #############

set.seed(314) # for reproducibility
X <- rStudent(n, df = nu, scale = P) # n multivariate t observations
U <- pt(X, df = nu) # n independent realizations from the corresponding t copula
Y <- qexp(U) # transform U (t copula) to Exp(1) margins

## Plot setup
ind <- c(A = 96, B = 696, C = 112) # use 'plot(X); identify(X)' to identify these points

## Scatter plot of a bivariate t distribution
plot(X, xlab = expression(X[1]), ylab = expression(X[2]))
points(X[ind["A"],1], X[ind["A"],2], bg = cols[1], col = cols[1])
points(X[ind["B"],1], X[ind["B"],2], bg = cols[2], col = cols[2])
points(X[ind["C"],1], X[ind["C"],2], bg = cols[3], col = cols[3])
text(X[ind["A"],1], X[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1], font = 2)
text(X[ind["B"],1], X[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2], font = 2)
text(X[ind["C"],1], X[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3], font = 2)

## Scatter plot of the corresponding t copula (1st part of Sklar's Theorem)
plot(U, xlab = expression(U[1]==F(X[1])), ylab = expression(U[2]==F(X[2])))
points(U[ind["A"],1], U[ind["A"],2], bg = cols[1], col = cols[1])
points(U[ind["B"],1], U[ind["B"],2], bg = cols[2], col = cols[2])
points(U[ind["C"],1], U[ind["C"],2], bg = cols[3], col = cols[3])
text(U[ind["A"],1], U[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1], font = 2)
text(U[ind["B"],1], U[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2], font = 2)
text(U[ind["C"],1], U[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3], font = 2)

## Scatter plot of the meta-t distribution with Exp(1) margins (2nd part of Sklar's Theorem)
plot(Y, xlab = expression(Y[1]=={F[1]^{-1}}(U[1])), ylab = expression(Y[2]=={F[2]^{-1}}(U[2])))
points(Y[ind["A"],1], Y[ind["A"],2], bg = cols[1], col = cols[1])
points(Y[ind["B"],1], Y[ind["B"],2], bg = cols[2], col = cols[2])
points(Y[ind["C"],1], Y[ind["C"],2], bg = cols[3], col = cols[3])
text(Y[ind["A"],1], Y[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1], font = 2)
text(Y[ind["B"],1], Y[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2], font = 2)
text(Y[ind["C"],1], Y[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3], font = 2)

par(opar) # reset plot parameters