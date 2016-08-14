## By Marius Hofert and Alexander J. McNeil

## Illustrating meta-C models


### 0 Setup ####################################################################

library(copula)
library(gridExtra) # for grid.arrange()

n <- 1000 # sample size
set.seed(271)


### 1 Generate the data ########################################################

## Define the copula parameters
## Note: They are chosen such that Pearson's correlation coefficient,
##       when computed with N(0,1) margins, is 0.7 (see below)
th.n <- 0.7 # Gauss copula parameter
th.g <- 2 # Gumbel copula parameter
th.c <- 2.2 # Clayton copula parameter
th.t <- 0.71 # t_4 copula parameter

## Define the copulas
nc <- normalCopula(th.n) # Gauss copula
gc <- gumbelCopula(th.g) # Gumbel copula
cc <- claytonCopula(th.c) # Clayton copula
tc <- tCopula(th.t) # t_4 copula

## Generate copula data
U.nc <- rCopula(n, copula = nc)
U.gc <- rCopula(n, copula = gc)
U.cc <- rCopula(n, copula = cc)
U.tc <- rCopula(n, copula = tc)

## Map to N(0,1) margins (meta-copula data)
X.nc <- qnorm(U.nc)
X.gc <- qnorm(U.gc)
X.cc <- qnorm(U.cc)
X.tc <- qnorm(U.tc)

## Correlations
cors <- vapply(list(X.nc, X.gc, X.cc, X.tc), function(x) cor(x)[1,2], NA_real_)
stopifnot(all.equal(cors, rep(0.7, 4), tol = 0.015))

## Define the corresponding (meta-C) densities (via Sklar's Theorem)
## Note: Density f(x_1, x_2) = c(F_1(x_1), F_2(x_2)) * f_1(x_1) * f_2(x_2)
##                           = exp( log(c(F_1(x_1), F_2(x_2))) + log(f_1(x_1)) + log(f_2(x_2)) )
dMetaCopulaN01 <- function(x, copula)
    exp(dCopula(pnorm(x), copula = copula, log = TRUE) + rowSums(dnorm(x, log = TRUE)))
## Alternatively, we could work with dMvdc() here


### 2 Scatter plots ############################################################

## Copula samples (see Figure 7.3)
opar <- par(pty = "s", mar = c(5.1, 4.1, 4.1, 2.1) - 1)
lay <- matrix(1:4, ncol = 2, byrow = TRUE) # layout matrix
layout(lay) # layout
plot(U.nc, xlab = expression(U[1]), ylab = expression(U[2]), # Gauss copula
     cex = 0.4, main = "Gauss copula sample")
plot(U.gc, xlab = expression(U[1]), ylab = expression(U[2]), # Gumbel copula
     cex = 0.4, main = "Gumbel copula sample")
plot(U.cc, xlab = expression(U[1]), ylab = expression(U[2]), # Clayton copula
     cex = 0.4, main = "Clayton copula sample")
plot(U.tc, xlab = expression(U[1]), ylab = expression(U[2]), # t_4 copula
     cex = 0.4, main = expression(bold(italic(t)[4]~"copula sample")))

## Meta-copula samples with N(0,1) margins (see Figure 7.4)
m <- max(abs(X.nc), abs(X.gc), abs(X.cc), abs(X.tc))
lim <- c(-m, m)
plot(X.nc, xlim = lim, ylim = lim, # meta-Gauss
     xlab = expression(X[1]), ylab = expression(X[2]),
     cex = 0.4, main = "Meta-Gauss sample")
mtext("N(0,1) margins", side = 4, line = 0.25, adj = 0, cex = 0.7)
plot(X.gc, xlim = lim, ylim = lim,
     xlab = expression(X[1]), ylab = expression(X[2]), # meta-Gumbel
     cex = 0.4, main = "Meta-Gumbel sample")
mtext("N(0,1) margins", side = 4, line = 0.25, adj = 0, cex = 0.7)
plot(X.cc, xlim = lim, ylim = lim, # meta-Clayton
     xlab = expression(X[1]), ylab = expression(X[2]),
     cex = 0.4, main = "Meta-Clayton sample")
mtext("N(0,1) margins", side = 4, line = 0.25, adj = 0, cex = 0.7)
plot(X.tc, xlim = lim, ylim = lim, # meta-t_4
     xlab = expression(X[1]), ylab = expression(X[2]),
     cex = 0.4, main = expression(bold("Meta-"*italic(t)[4]~"sample")))
mtext("N(0,1) margins", side = 4, line = 0.25, adj = 0, cex = 0.7)
par(opar) # restore graphical parameters


### 3 Density plots ############################################################

## Wire frame plots
n.grid <- 26 # number of grid points
lim <- c(-2.5, 2.5)
s <- seq(lim[1], lim[2], length.out = n.grid)
grid <- as.matrix(expand.grid("x[1]" = s, "x[2]" = s))
val.n <- cbind(grid, "f(x[1],x[2])" = dMetaCopulaN01(grid, copula = nc))
val.g <- cbind(grid, "f(x[1],x[2])" = dMetaCopulaN01(grid, copula = gc))
val.c <- cbind(grid, "f(x[1],x[2])" = dMetaCopulaN01(grid, copula = cc))
val.t <- cbind(grid, "f(x[1],x[2])" = dMetaCopulaN01(grid, copula = tc))
zlim <- c(0, max(val.n[,3], val.g[,3], val.c[,3], val.t[,3]))

## Density plots
zm <- 1.1 # zoom in a bit more
scs <- list(arrows = FALSE, col = "black", cex = 0.6) # scale down ticks
tcx <- 0.95 # scale back titles
w.n <- wireframe2(val.n, xlim = lim, ylim = lim, zlim = zlim, zoom = zm, scales = scs,
                  main = list(label = "Meta-Gauss density with N(0,1) margins", cex = tcx))
w.g <- wireframe2(val.g, xlim = lim, ylim = lim, zlim = zlim, zoom = zm, scales = scs,
                  main = list(label = "Meta-Gumbel density with N(0,1) margins", cex = tcx))
w.c <- wireframe2(val.c, xlim = lim, ylim = lim, zlim = zlim, zoom = zm, scales = scs,
                  main = list(label = "Meta-Clayton density with N(0,1) margins", cex = tcx))
w.t <- wireframe2(val.t, xlim = lim, ylim = lim, zlim = zlim, zoom = zm, scales = scs,
                  main = list(label = expression(bold("Meta-"*italic(t)[4]~"density with N(0,1) margins")),
                              cex = tcx))
grid.arrange(w.n, w.g, w.c, w.t, ncol = 2)

## Contour plots
tcx <- 0.9
c.n <- contourplot2(val.n, xlim = lim, ylim = lim,
                    main = list(label = "Meta-Gauss contours with N(0,1) margins", cex = tcx))
c.g <- contourplot2(val.g, xlim = lim, ylim = lim,
                    main = list(label = "Meta-Gumbel contours with N(0,1) margins", cex = tcx))
c.c <- contourplot2(val.c, xlim = lim, ylim = lim,
                    main = list(label = "Meta-Clayton contours with N(0,1) margins", cex = tcx))
c.t <- contourplot2(val.t, xlim = lim, ylim = lim,
                    main = list(label = expression(bold("Meta-"*italic(t)[4]~"contours with N(0,1) margins")),
                                cex = tcx))
grid.arrange(c.n, c.g, c.c, c.t, ncol = 2)
