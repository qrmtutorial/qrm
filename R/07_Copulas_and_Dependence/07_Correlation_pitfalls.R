## By Marius Hofert

## Selected correlation pitfalls


library(copula)


### Fallacy 1 (Ex. 1): F_1, F_2 and correlation rho uniquely determine F #######

## Simple example with same margins, rho = 0, but (obviously) different models
n <- 1000
set.seed(271)
Z <- rnorm(n)
U <- runif(n)
V <- rep(1, n)
V[U < 1/2] <- -1 # => V in {-1,1}, each with probability 1/2
X <- cbind(Z, Z*V) # convex combination of W and M (each with prob. 1/2 => rho = 0)
cor(X)[2,1]
Y <- matrix(rnorm(n * 2), ncol = 2) # independent N(0,1) (=> rho = 0)
cor(Y)[2,1]

## Plots
plot(X, xlab = expression(X[1]), ylab = expression(X[2]))
plot(Y, xlab = expression(Y[1]), ylab = expression(Y[2]))


### Fallacy 1 (Ex. 2): F_1, F_2 and correlation rho uniquely determine F #######

## We construct a parametric copula family here which *always* produces
## correlation 0. To this end, let C(u1,u2) = u1 * u2 + f1(u1) * f2(u2),
## for f1, f2 continuously differentiable s.t. f1(0) = f2(0) = f1(1) = f2(1) = 0
## and s.t. 1 + f1'(u1) * f2'(u2) >= 0 for all u1, u2 in [0,1].
## By computing C's density, C can be seen to be a copula. Let (U_1,U_2) ~ C.
## As the margins are U(0,1), Var(U_1) = Var(U_2) = 1/12 and thus Pearson's
## correlation coefficient rho equals rho = 12 Cov(U_1, U_2). By Hoeffding's
## identity, rho = 12 \int_0^1 \int_0^1 (C(u1,u2) - u1 * u2) du1 du2. Plugging
## in our form of C, we obtain rho = 12 \int_0^1 f1(u1) du1 * \int_0^1 f2(u2) du2.
## So if one of the f's is point symmetric about 1/2, rho equals 0. For example,
## take f1(x) = 2x(x-1/2)(x-1) and f2(x) = theta * x * (1-x). Then, for any
## theta in [-1,1], f1 and f2 fulfill all conditions and will generate a
## proper copula with Pearson's correlation coefficient equal to 0.
## Let's do that!

## Define auxiliary functions
f1 <- function(x) 2*x*(x-1/2)*(x-1) # f1
f1. <- function(x) 6*x^2-6*x+1 # derivative of f1
f2 <- function(x, th) th * x * (1-x) # f2
f2. <- function(x, th) th * (-2 * x + 1) # derivative of f2

## Define the conditional copula C(u2 | u_1)
cCop <- function(u, th) {
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot(-1 <= th, th <= 1)
    u[,2] + f1.(u[,1]) * f2(u[,2], th = th)
}

## Density
dCop <- function(u, th) {
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot(-1 <= th, th <= 1)
    1 + f1.(u[,1]) * f2.(u[,2], th = th)
}

## Conditional distribution method for this copula
## Note; We numerically invert C(u2 | u_1) to get C^{-1}(u2 | u_1)
CDM <- function(u, th) {
    if(!is.matrix(u)) u <- rbind(u)
    cbind(u[,1], apply(u, 1, function(u.)
                     uniroot(function(u2) cCop(c(u.[1], u2), th = th) - u.[2],
                             interval = 0:1, maxiter = 1000)$root))
}

## Sample
th <- 1
n <- 10000
set.seed(1)
U <- CDM(matrix(runif(2*n), ncol = 2), th = th)
plot(U, xlab = expression(U[1]), ylab = expression(U[2]))

## 'Show' that rho is (close to) 0
stopifnot(all.equal(cor(U)[2,1], 0.0016, tol = 1e-2))

## Visualize the copula density
n.grid <- 26 # number of grid points
u <- seq(0, 1, length.out = n.grid) # subdivison points in each dimension
grid <- expand.grid("u[1]" = u, "u[2]" = u) # build a grid
val.c <- cbind(grid, "c(u[1],u[2])" = dCop(grid, th = th)) # evaluate density on the grid
wireframe2(val.c) # wire frame plot of the density
contourplot2(val.c, xlim = 0:1, ylim = 0:1) # level curves


### Fallacy 2: Given F_1, F_2, any rho in [-1,1] is attainable #################

## Function to compute the correlation bounds for LN(0, sigma_.^2) margins
## Note: The derivation can be done with the moment-generating function of
##       the standard normal distribution
cor_bound_LN <- function(s, method = c("max", "min")) {
    ## s = (sigma_1, sigma_2)
    if(!is.matrix(s)) s <- rbind(s)
    method <- match.arg(method)
    if(method == "min") s[,2] <- -s[,2]
    (exp((s[,1]+s[,2])^2/2)-exp((s[,1]^2+s[,2]^2)/2)) /
        sqrt(expm1(s[,1]^2)*exp(s[,1]^2)*expm1(s[,2]^2)*exp(s[,2]^2))
}

## Evaluate correlation bounds on a grid
n.grid <- 26 # number of grid points in each dimension
s <- seq(0.01, 5, length.out = n.grid) # subdivision points in each dimension
grid <- expand.grid("sigma[1]" = s, "sigma[2]" = s) # build a grid
val.min <- cbind(grid, "underline(Cor)(sigma[1],sigma[2])" =
                 cor_bound_LN(grid, method = "min"))
val.max <- cbind(grid, "bar(Cor)(sigma[1],sigma[2])" =
                 cor_bound_LN(grid))

## Plots
wireframe2(val.min)
wireframe2(val.max)
