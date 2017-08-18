## By Marius Hofert

## Selected copulas and copula families


### Setup ######################################################################

library(copula)

n <- 1000 # sample size
d <- 5 # max. considered dimension
tau <- 0.5 # Kendall's tau
tblack <- function(alpha.f) adjustcolor("black", alpha.f = alpha.f) # color (transparent black)

set.seed(271) # for reproducibility


### 1 Fundamental copulas ######################################################

### 1.1 Independence copula ####################################################

## Define the independence copula object
ic <- indepCopula()

## Copula (wire frame and level curves)
wireframe2(ic, FUN = pCopula)
contourplot2(ic, FUN = pCopula)

## Copula density (wire frame)
wireframe2(ic, FUN = dCopula, delta = 0.001, zlim = 0:1)

## Scatter plots
U <- matrix(runif(n*d), ncol = 5)
plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
cloud2(U[,1:3], # d = 3
       xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
pairs2(U) # d = 5


### 1.2 Frechet-Hoeffding bounds W and M #######################################

## Copulas (wire frame and level curves)
n.grid <- 26 # number of grid points
u <- seq(0, 1, length.out = n.grid) # subdivison points in each dimension
grid <- expand.grid("u[1]" = u, "u[2]" = u) # build a grid
W <- pmax(grid[,1] + grid[,2] - 1, 0) # values of W on grid
M <- pmin(grid[,1], grid[,2]) # values of M on grid
val.W <- cbind(grid, "W(u[1],u[2])" = W) # append grid
val.M <- cbind(grid, "M(u[1],u[2])" = M) # append grid
wireframe2(val.W) # wire frame plot of W

contourplot2(val.W, xlim = 0:1, ylim = 0:1) # level curves of W
wireframe2(val.M) # wire frame plot of M
contourplot2(val.M, xlim = 0:1, ylim = 0:1) # level curves of M

## Scatter plots
U <- runif(n)
plot(cbind(U, 1-U), xlab = expression(U[1]), ylab = expression(U[2]), col = tblack(0.5)) # sample of W for d = 2
plot(cbind(U, U), xlab = expression(U[1]), ylab = expression(U[2]), col = tblack(0.5)) # sample of M for d = 2
cloud2(do.call(cbind, rep(list(U), 3)), col = tblack(0.05), # sample of M for d = 3
       xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
pairs2(do.call(cbind, rep(list(U), d)), col = tblack(0.05)) # sample of M for d = 5


### 2 Implicit copulas #########################################################

### 2.1 Normal (or: Gauss) copula ##############################################

## Define the normal copula object
th <- iTau(normalCopula(), tau = tau)
nc <- normalCopula(th)

## Copula (wire frame and level curves)
wireframe2(nc, FUN = pCopula)
contourplot2(nc, FUN = pCopula)

## Copula density (wire frame and level curves)
wireframe2(nc, FUN = dCopula, delta = 0.02)
contourplot2(nc, FUN = dCopula)

## Scatter plots
nc. <- normalCopula(th, dim = 5)
U <- rCopula(n, copula = nc.)
plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
cloud2(U[,1:3], # d = 3
       xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
pairs2(U, cex = 0.4, col = tblack(0.5)) # d = 5


### 2.2 t copula ###############################################################

## Define the t copula object
nu <- 3 # degrees of freedom
th <- iTau(tCopula(, df = nu), tau = tau) # correlation parameter
tc <- tCopula(th, df = nu)

## Copula (wire frame and level curves)
## Note: pCopula() is only available for integer degrees of freedom
wireframe2(tc, FUN = pCopula)
contourplot2(tc, FUN = pCopula)

## Copula density (wire frame and level curves)
wireframe2(tc, FUN = dCopula, delta = 0.02)
contourplot2(tc, FUN = dCopula)

## Scatter plots
tc. <- tCopula(th, dim = 5, df = nu)
U <- rCopula(n, copula = tc.)
plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
cloud2(U[,1:3], # d = 3
       xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
pairs2(U, cex = 0.4, col = tblack(0.5)) # d = 5

## A nonexchangeable example
## Build a block correlation matrix
rho <- c(0.3, 0.6, 0.9)
P <- matrix(rho[1], nrow = d, ncol = d)
P[1:2, 1:2] <- rho[2]
P[3:5, 3:5] <- rho[3]
diag(P) <- 1
## Define, sample and plot a t copula
tc.. <- tCopula(P2p(P), dim = d, dispstr = "un", df = 3.5)
U. <- rCopula(n, copula = tc..)
pairs2(U., cex = 0.4, col = tblack(0.5)) # d = 5

## A model more flexible than the multivariate t (more degrees of freedom)
X <- sapply(1:ncol(U.), function(j) qt(U.[,j], df = j))
pairs2(X, labels.null.lab = "X", cex = 0.4, col = tblack(0.5))


### 3 Explicit copulas #########################################################

### 3.1 Clayton copula #########################################################

## Define the Clayton copula object
th <- iTau(claytonCopula(), tau = tau)
cc <- claytonCopula(th)

## Copula (wire frame and level curves)
wireframe2(cc, FUN = pCopula)
contourplot2(cc, FUN = pCopula)

## Copula density (wire frame and level curves)
wireframe2(cc, FUN = dCopula, delta = 0.02)
contourplot2(cc, FUN = dCopula)

## Scatter plots
cc. <- claytonCopula(th, dim = 5)
U <- rCopula(n, copula = cc.)
plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
cloud2(U[,1:3], # d = 3
       xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
pairs2(U, cex = 0.4, col = tblack(0.5)) # d = 5


### 3.2 Gumbel copula ##########################################################

## Define the Gumbel copula object
th <- iTau(gumbelCopula(), tau = tau)
gc <- gumbelCopula(th)

## Copula (wire frame and level curves)
wireframe2(gc, FUN = pCopula)
contourplot2(gc, FUN = pCopula)

## Copula density (wire frame and level curves)
wireframe2(gc, FUN = dCopula, delta = 0.02)
contourplot2(gc, FUN = dCopula)

## Scatter plots
gc. <- gumbelCopula(th, dim = 5)
U <- rCopula(n, copula = gc.)
plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
cloud2(U[,1:3], # d = 3
       xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
pairs2(U, cex = 0.4, col = tblack(0.5)) # d = 5

## Scatter plot of a survival Gumbel copula
pairs2(1-U, cex = 0.4, col = tblack(0.5))

## Only the first component flipped
pairs2(U, cex = 0.4, col = tblack(0.5)) # original
pairs2(cbind(1-U[,1], U[,2:5]), cex = 0.4, col = tblack(0.5)) # first flipped


### 3.3 Outer power Clayton copula #############################################

## Note: Outer power Clayton copulas can have both lower and upper tail dependence

## Define the outer power Clayton copula
tauC <- 0.3 # Kendall's tau for the underlying Clayton copula
thC <- copClayton@iTau(tauC) # choose Clayton's generator s.t. Kendall's tau is tauC
opC <- opower(copClayton, thC) # define an outer power Clayton copula (its parameter theta is not specified yet)
th <- opC@iTau(tau) # define the outer power Clayton copula s.t. has Kendall's tau is tau
opcc <- onacopulaL(opC, list(th, 1:2)) # define the outer power Clayton copula

## Scatter plot
U <- rCopula(n, copula = opcc)
plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2


### 3.4 Marshall--Olkin copulas ################################################

## Note: Marshall--Olkin copulas have a singular component, i.e., a set of
##       Lebesgue measure 0 with positive probability mass assigned.

## Define the MO copula
C <- function(u, alpha) {
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot(0 <= alpha, alpha <= 1, length(alpha) == 2)
    pmin(u[,1]*u[,2]^(1-alpha[2]), u[,1]^(1-alpha[1])*u[,2])
}

## Sampling
rMO <- function(n, alpha) {
    stopifnot(n >= 1, 0 <= alpha, alpha <= 1, length(alpha) == 2)
    U. <- matrix(runif(n*3), ncol = 3)
    U <- cbind(pmax(U.[,1]^(1/(1-alpha[1])), U.[,3]^(1/alpha[1])),
               pmax(U.[,2]^(1/(1-alpha[2])), U.[,3]^(1/alpha[2])))
}

## Define the singular component
S.C <- function(u, alpha) {
    stopifnot(0 <= u, u <= 1,
              0 <= alpha, alpha <= 1, length(alpha) == 2)
    tau <- alpha[1] * alpha[2] / (alpha[1] + alpha[2] - alpha[1] * alpha[2])
    tau * pmin(u[,1]^alpha[1], u[,2]^alpha[2])^(1/tau)
}

## Define the absolutely continuous component
A.C <- function(u, alpha) C(u, alpha) - S.C(u, alpha)


## Scatter plot
alpha <- c(0.2, 0.8)
U <- rMO(n, alpha = alpha)
plot(U, xlab = expression(U[1]), ylab = expression(U[2]))
## Interpretation: Given u_1 = 0.6, for example, u_2 lies with a non-zero
## probability p on the curve and with the remaining probability 1-p
## anywhere else (*not*: uniform) along the line u_1 = 0.6.

## Check the margins
plot(U[,1], ylab = expression(U[1]))
plot(U[,2], ylab = expression(U[2]))

## Evaluate the copula (and its singular and absolutely continuous components)
grid <- expand.grid("u[1]" = u, "u[2]" = u) # build a grid
val.C <- cbind(grid, "C(u[1],u[2])" = C(grid, alpha = alpha)) # append C
val.S <- cbind(grid, "S[C](u[1],u[2])" = S.C(grid, alpha = alpha)) # append S.C
val.A <- cbind(grid, "A[C](u[1],u[2])" = A.C(grid, alpha = alpha)) # append S.C

## Copula (wire frame and level curves)
wireframe2(val.C) # wire frame plot
contourplot2(val.C, xlim = 0:1, ylim = 0:1) # level curves
## The copula has a kink, that is, a so-called singular component.
## A singular component is a set of Lebesgue measure 0 where C puts
## mass at (here: a curve). This is better visible from the scatter plots below.

## Singular and absolutely continuous component (wire frame)
wireframe2(val.S) # singular component
wireframe2(val.A) # absolutely continuous component
