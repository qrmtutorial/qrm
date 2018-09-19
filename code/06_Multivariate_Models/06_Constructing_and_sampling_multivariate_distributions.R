## By Marius Hofert

## Understanding the construction of multivariate normal distributions,
## normal variance mixtures, elliptical distributions etc.
## Note: This is for educational purposes only; for several of the tasks
##       presented below, functions already exist in R packages.


### 1 Basic 'ingredients' ######################################################

## Building a covariance/correlation matrix

## Covariance matrix
A <- matrix(c( 4, 0,
              -1, 1), ncol = 2, byrow = TRUE) # Cholesky factor of the ...
(Sigma <- A %*% t(A)) # ... symmetric, positive definite (covariance) matrix Sigma

## Corresponding correlation matrix
P <- outer(1:2, 1:2, Vectorize(function(r, c)
           Sigma[r,c]/(sqrt(Sigma[r,r])*sqrt(Sigma[c,c])))) # construct the corresponding correlation matrix
(P.  <- cov2cor(Sigma)) # a more elegant solution, see the source of cov2cor()
stopifnot(all.equal(P., P))
## Another option would be as.matrix(Matrix::nearPD(Sigma, corr = TRUE, maxit = 1000)$mat)
## which works differently, though (by finding a correlation matrix close to the
## given matrix in the Frobenius norm) and thus gives a different answer.


## Decomposing a covariance/correlation matrix

## We frequently need to decompose a covariance matrix Sigma (or correlation
## matrix P) as Sigma = AA^T. To this end there are several possibilities.

## The Cholesky decomposition
A. <- t(chol(Sigma)) # the Cholesky factor (lower triangular with real entries > 0)
stopifnot(all.equal(A. %*% t(A.) , Sigma), # checking decomposition
          all.equal(A., A)) # checking uniqueness of the Cholesky decomposition

## Other decompositions of Sigma than A %*% t(A) are possible, too
if(FALSE) {
    ## Eigendecomposition (or spectral decomposition)
    eig <- eigen(Sigma) # eigenvalues and eigenvectors
    V <- eig$vectors # matrix of eigenvectors
    Lambda <- diag(pmax(eig$values, 0))
    stopifnot(all.equal(Sigma, V %*% Lambda %*% t(V))) # for real, symmetric matrices
    A.eig <- V %*% sqrt(Lambda) %*% t(V)
    Sigma.eig <- A.eig %*% t(A.eig)
    stopifnot(all.equal(Sigma.eig, Sigma))
    ## ... but A.. (non-triangular) and A (triangular) are different

    ## Singular-value decomposition
    sv <- svd(Sigma) # singular values, U, V (left/right singular vectors) such that Sigma = U diag(<singular values>) V^T
    A.sv <- sv$u %*% sqrt(diag(pmax(sv$d, 0))) %*% t(sv$v)
    Sigma.sv <- A.sv %*% t(A.sv)
    stopifnot(all.equal(Sigma.sv, Sigma))
}


### 2 Sampling from the multivariate normal distribution #######################

## Sampling Z (iid N(0,1))
n <- 1000 # sample size
d <- 2 # dimension
set.seed(271) # set a seed (for reproducibility)
Z <- matrix(rnorm(n * d), ncol = d) # sample iid N(0,1) random variates
mabs <- max(abs(Z))
lim <- c(-mabs, mabs)
plot(Z, xlim = lim, ylim = lim,
     xlab = expression(Z[1]), ylab = expression(Z[2]))

## Change the covariance matrix from identity to Sigma
X. <- t(A %*% t(Z)) # ~ N_d(0, Sigma)
xlim <- c(-16, 16)
ylim <- c(-8, 8)
plot(X., xlim = xlim, ylim = ylim,
     xlab = expression(X[1]), ylab = expression(X[2]))

## Use a shift to get X ~ N_d(mu, Sigma)
mu <- c(1, 3)
X.norm <- rep(mu, each = n) + X.
plot(X.norm, xlim = xlim, ylim = ylim,
     xlab = expression(X[1]), ylab = expression(X[2]))
## Think about why the rep() and why the t()'s above!
## See also Hofert (2013, "On Sampling from the Multivariate t Distribution")

## Plot (even smaller correlation (larger absolute correlation))
P. <- matrix(c(    1, -0.95,
               -0.95,  1), ncol = d, byrow = TRUE)
Sigma. <- outer(1:d, 1:d, Vectorize(function(r, c)
                P.[r,c] * (sqrt(Sigma[r,r])*sqrt(Sigma[c,c])))) # construct the corresponding Sigma
## Note: When manually changing the off-diagonal entry of Sigma make sure that
##       |Sigma.[1,2]| <= sqrt(Sigma.[1,1]) * sqrt(Sigma.[2,2]) still holds
A. <- t(chol(Sigma.)) # new Cholesky factor
X.norm. <- rep(mu, each = n) + t(A. %*% t(Z)) # generate the sample
plot(rbind(X.norm, X.norm.), xlim = xlim, ylim = ylim,
     xlab = expression(X[1]), ylab = expression(X[2]),
     col = rep(c("black", "royalblue3"), each = n))
legend("bottomright", bty = "n", pch = c(1,1), col = c("black", "royalblue3"),
       legend = c(expression(N(mu,Sigma)),
                  expression("With larger absolute correlation")))

## Plot (switch sign of correlation)
Sigma. <- Sigma * matrix(c( 1, -1,
                           -1,  1), ncol = d, byrow = TRUE)
A. <- t(chol(Sigma.))
X.norm. <- rep(mu, each = n) + t(A. %*% t(Z))
plot(rbind(X.norm, X.norm.), xlim = xlim, ylim = ylim,
     xlab = expression(X[1]), ylab = expression(X[2]),
     col = rep(c("black", "royalblue3"), each = n))
legend("bottomright", bty = "n", pch = c(1,1), col = c("black", "royalblue3"),
       legend = c(expression(N(mu,Sigma)),
                  expression("Correlation sign switched")))

## Note: This has all been done based on the same Z (and we also recycle it below)!


### 3 Sampling from the multiv. t distribution (and other normal variance mixtures)

## We recycle the sample with positive correlation
X.norm <- X.norm.
A <- A.
Sigma <- A %*% t(A)


### 3.1 A sample from a multivariate t_nu distribution #########################

nu <- 3 # degrees of freedom
set.seed(314)
W <- 1/rgamma(n, shape = nu/2, rate = nu/2) # mixture variable for a multiv. t_nu distribution
plot(W, log = "y")
X.t <- rep(mu, each = n) + sqrt(W) * t(A %*% t(Z))
xlim <- c(-40, 80)
ylim <- c(-20, 60)
plot(rbind(X.t, X.norm), xlim = xlim, ylim = ylim,
     xlab = expression(X[1]), ylab = expression(X[2]),
     col = rep(c("maroon3", "black"), each = n))
legend("bottomright", bty = "n", pch = c(1,1), col = c("black", "maroon3"),
       legend = c(expression(N(mu,Sigma)),
                  expression(italic(t)[nu](mu,Sigma))))
## For more about sampling from the multivariate t distribution (especially
## what can go wrong!), see Hofert (2013, "On Sampling from the Multivariate t Distribution")


### 3.2 A 2-point distribution for W (= two 'overlaid' multiv. normals) ########

## Let's change the distribution of W to a 2-point distribution with
## W = 1 with prob. 0.5 and W = 2 with prob. 0.5. What we will see are two 'overlaid'
## normal distributions (as each of W = 1 and W = 2 corresponds to one).
W.binom <- 1 + rbinom(n, size = 1, prob = 0.5)
cols <- rep("royalblue3", n)
cols[W.binom == 2] <- "maroon3"
X.W.binom <- rep(mu, each = n) + sqrt(W.binom) * t(A %*% t(Z))
plot(rbind(X.t, X.W.binom), xlab = expression(X[1]), ylab = expression(X[2]),
     col = c(rep("black", n), cols))
legend("topright", bty = "n", pch = rep(1, 3), col = c("royalblue3", "maroon3", "black"),
       legend = c(expression(N(mu,W*Sigma)~"for W = 1"),
                  expression(N(mu,W*Sigma)~"for W = 2"),
                  expression(italic(t)[nu](mu,Sigma))))
## => With probability 0.5 we sample from a normal variance mixture with W = 1
##    and with probability 0.5 we sample from one with W = 2. By using a df for
##    W with infinite upper endpoint, we can reach further out in the tails than
##    with any multivariate normal distribution by overlaying normals
##    with different (unbounded) covariance matrices (and this is what is
##    creating heavier tails for the t distribution than the normal).


### 3.3 An *uncorrelated* multivariate t_nu distribution vs an independent one

A. <- diag(sqrt(diag(Sigma))) # Cholesky factor of the diagonal matrix containing the variances
X.t.uncor <- rep(mu, each = n) + sqrt(W) * t(A. %*% t(Z)) # uncorrelated t samples
X.t.ind   <- rep(mu, each = n) + matrix(rt(n * d, df = nu), ncol = d) # independent t samples
plot(rbind(X.t.uncor, # uncorrelated t; dependence only introduced by sqrt(W) => shifted spherical distribution
           X.t, # full t; dependence introduced by sqrt(W) *and* correlation
           X.t.ind), # independence
     xlim = xlim, ylim = ylim,
     xlab = expression(X[1]), ylab = expression(X[2]),
     col = rep(c("royalblue3", "maroon3", "black"), each = n))
legend("topright", bty = "n", pch = rep(1, 3), col = c("black", "royalblue3", "maroon3"),
       legend = c(expression("Independent"~italic(t)[nu]),
                  expression("Uncorrelated"~italic(t)[nu]~"(shifted spherical)"),
                  expression(italic(t)[nu](mu,Sigma))))
## => The uncorrelated (but not independent) samples show more mass in the joint
##    tails than the independent samples. The (neither uncorrelated nor independent)
##    t samples show the most mass in the joint tails.


### 3.4 A normal *mean*-variance mixture #######################################

## Now let's also 'mix' the mean (replace mu by mu(W)), to get a normal
## mean-variance mixture; here: choosing between two different locations mu
## depending on W.
mu. <- t(sapply(W.binom, function(w) mu + if(w == 1) 0 else c(15, 0)))
X.mean.var <- mu. + sqrt(W.binom) * t(A %*% t(Z))
plot(rbind(X.W.binom, X.mean.var), xlab = expression(X[1]), ylab = expression(X[2]),
     col = c(rep("black", n), cols))
legend("topleft", bty = "n", pch = rep(1, 3), col = c("black", "royalblue3", "maroon3"),
       legend = c(expression(N(mu,W*Sigma)~"for W"%in%"{1,2}"),
                  "Mean-var. mix for W = 1",
                  "Mean-var. mix for W = 2"))
## => Not only do the red points show more variation, they *also* have a different
##    location; clearly, this sample (blue + red points together) is not
##    elliptically distributed anymore.


### 4 Sampling elliptical distributions ########################################

### 4.1 Sample from a multivariate t_nu via elliptical distributions ###########

## Sample from the uniform distribution on the unit sphere (spherical distribution)
set.seed(271)
Z <- matrix(rnorm(n * d), ncol = d)
S <- Z/sqrt(rowSums(Z^2)) # Z/||Z||
lim <- c(-1.3, 1.3)
par(pty = "s")
plot(S, xlim = lim, ylim = lim,
     xlab = expression(S[1]), ylab = expression(S[2]))

## Add scale (elliptical distribution) by computing X = RAS (= ARS)
rho <- 0.8 # correlation
P <- matrix(c(  1, rho,
              rho, 1), ncol = 2) # correlation matrix
A <- chol(P) # Cholesky factor (here: upper triangular matrix for efficiency reasons)
X. <- t(A %*% t(S)) # (n, d) matrix
par(pty = "s")
plot(X., xlim = lim, ylim = lim,
     xlab = expression(X[1]), ylab = expression(X[2]))

## Add radial part of a t_nu distribution (spherical distribution)
nu <- 5
R <- sqrt(d * rf(n, df1 = d, df2 = nu))
X.. <- R * X.
ran <- c(min(X.., -max(X..)), max(X.., -min(X..)))
par(pty = "s")
plot(X.., xlim = ran, ylim = ran, xlab = expression(X[1]), ylab = expression(X[2]))

## Add location (elliptical distribution) by computing X = mu + RAS
mu <- c(ran[2]-max(X..[,1]), ran[2]-max(X..[,2])) # push sample towards the upper right corner
X <- rep(mu, each = n) + X..
par(pty = "s")
plot(X, xlim = ran, ylim = ran, xlab = expression(X[1]), ylab = expression(X[2]))


### 4.2 A 2-point distribution for R ###########################################

R <- 1 + rbinom(n, size = 1, prob = 2/3) # P(R = 1) = 1/3, P(R = 2) = 2/3
cols <- rep("royalblue3", n)
cols[R == 2] <- "maroon3" # prob. 2/3
X <- rep(mu, each = n) + t(A %*% t(R*S)) # same mu, A, S
par(pty = "s")
plot(X, xlab = expression(X[1]), ylab = expression(X[2]), col = cols)
legend("bottomright", bty = "n", pch = c(1,1), col = c("royalblue3", "maroon3"),
       legend = c("Elliptical for R = 1", "Elliptical for R = 2"))
## ... in other words, the radial part R scales each point on the ellipse AS
## to lie on another ellipse. For continuously distributed R, these ellipses
## are all overlaid (and thus not visible anymore), but for discrete R, we
## may see them (here: two ellipses).
