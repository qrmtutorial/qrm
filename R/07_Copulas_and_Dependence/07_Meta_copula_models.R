### By Marius Hofert and Alexander J. McNeil


### Copula illustrations #######################################################

require(copula) # copula functionality
require(QRM) # for BiDensPlot

## Define copulas
n <- 2000 # sample size
d <- 2 # dimension
Ga.theta <- 0.7 # Gauss copula parameter
G.theta <- 2 # Gumbel copula parameter
C.theta <- 2.2 # Clayton copula parameter
t4.theta <- 0.71 # t_4 copula parameter
Ga.cop <- ellipCopula("normal",  param=Ga.theta, dim=d) # define Ga copula
G.cop  <- archmCopula("Gumbel",  param=G.theta,  dim=d) # define Gumbel copula
C.cop  <- archmCopula("Clayton", param=C.theta,  dim=d) # define Clayton copula
t4.cop <- ellipCopula("t",       param=t4.theta, dim=d, df=4) # define t_4 copula


### Copula scatter plots #######################################################

## Generate data
set.seed(271)
U.Ga <- rCopula(n, copula=Ga.cop)
U.G  <- rCopula(n, copula=G.cop)
U.C  <- rCopula(n, copula=C.cop)
U.t4 <- rCopula(n, copula=t4.cop)

## Scatter plots (similar to Figure 7.3)
par(mfrow=c(2, 2))
plot(U.Ga, xlab=expression(U[1]), ylab=expression(U[2]), main="Gauss")
plot(U.G,  xlab=expression(U[1]), ylab=expression(U[2]), main="Gumbel")
plot(U.C,  xlab=expression(U[1]), ylab=expression(U[2]), main="Clayton")
plot(U.t4, xlab=expression(U[1]), ylab=expression(U[2]),
     main=expression(bold(t[4])))
par(mfrow=c(1, 1))


### Copula density plots #######################################################

## Density plots (similar to Figure 7.5, 7.6)
lim <- c(0, 1)
par(mfrow=c(2, 2))
BiDensPlot(function(u) dCopula(u, copula=Ga.cop), xpts=lim, ypts=lim)
mtext("Gauss copula density")
BiDensPlot(function(u) dCopula(u, copula=G.cop), xpts=lim, ypts=lim)
mtext("Gumbel copula density")
BiDensPlot(function(u) dCopula(u, copula=C.cop), xpts=lim, ypts=lim)
mtext("Clayton copula density")
BiDensPlot(function(u) dCopula(u, copula=t4.cop), xpts=lim, ypts=lim)
mtext(t[4]~"copula density")
par(mfrow=c(1, 1))


### Meta-copula scatter plots with N(0,1) margins ##############################

## Generate meta-data
U.Ga.N01  <- qnorm(U.Ga)
U.G.N01   <- qnorm(U.G)
U.C.N01   <- qnorm(U.C)
U.t4.N01  <- qnorm(U.t4)

## Scatter plots (similar to Figure 7.4)
par(mfrow=c(2, 2))
plot(U.Ga.N01, xlab=expression(U[1]), ylab=expression(U[2]), main="Meta-Gauss")
plot(U.G.N01,  xlab=expression(U[1]), ylab=expression(U[2]), main="Meta-Gumbel")
plot(U.C.N01,  xlab=expression(U[1]), ylab=expression(U[2]), main="Meta-Clayton")
plot(U.t4.N01, xlab=expression(U[1]), ylab=expression(U[2]),
     main=expression(bold("Meta-"*t[4])))
par(mfrow=c(1, 1))


### Meta-copula plots with N(0,1) margins ######################################

## Meta-copula densities with N(0,1) margins
dmeta_G_N01 <- function(x, theta)
    exp(dCopula(pnorm(x), copula=G.cop, log=TRUE) + rowSums(dnorm(x, log=TRUE)))
dmeta_C_N01 <- function(x, theta)
    exp(dCopula(pnorm(x), copula=C.cop, log=TRUE) + rowSums(dnorm(x, log=TRUE)))
dmeta_t4_N01 <- function(x, theta)
    exp(dCopula(pnorm(x), copula=t4.cop, log=TRUE) + rowSums(dnorm(x, log=TRUE)))

## Meta-density plots with N(0,1) margins
lim <- c(-3, 3)
par(mfrow=c(2, 2))
BiDensPlot(dmnorm, xpts=lim, ypts=lim, mu=c(0,0), Sigma=equicorr(2, Ga.theta))
mtext("Meta-Gauss density with N(0,1) margins")
BiDensPlot(dmeta_G_N01, xpts=lim, ypts=lim, npts=64, theta=G.theta)
mtext("Meta-Gumbel density with N(0,1) margins")
BiDensPlot(dmeta_C_N01, xpts=lim, ypts=lim, npts=64, theta=C.theta)
mtext("Meta-Clayton density with N(0,1) margins")
BiDensPlot(dmeta_t4_N01, xpts=lim, ypts=lim, npts=64, theta=t4.theta)
mtext("Meta-"*t[4]~"density with N(0,1) margins")
par(mfrow=c(1, 1))

## Meta-contour plots with N(0,1) margins
par(mfrow=c(2, 2))
BiDensPlot(dmnorm, type="contour", xpts=lim, ypts=lim, mu=c(0,0),
           Sigma=equicorr(2, Ga.theta))
mtext("Meta-Gauss contours with N(0,1) margins", line=1)
BiDensPlot(dmeta_G_N01, type="contour", xpts=lim, ypts=lim, npts=64,
           theta=G.theta)
mtext("Meta-Gumbel contours with N(0,1) margins", line=1)
BiDensPlot(dmeta_C_N01, type="contour", xpts=lim, ypts=lim, npts=64,
           theta=C.theta)
mtext("Meta-Clayton contours with N(0,1) margins", line=1)
BiDensPlot(dmeta_t4_N01, type="contour", xpts=lim, ypts=lim, npts=64,
           theta=t4.theta)
mtext("Meta-"*t[4]~"contours with N(0,1) margins", line=1)
par(mfrow=c(1, 1))
