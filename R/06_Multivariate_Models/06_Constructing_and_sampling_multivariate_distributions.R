## By Marius Hofert

## Playground for multivariate distributions, in particular, the multivariate
## normal distribution, normal variance mixtures, elliptical distributions etc.
## Note that this is for education purposes (for several of the tasks presented
## here, functions already exist in R (packages)).


### 1 Basic 'ingredients' ######################################################

## Building a covariance/correlation matrix

## Covariance matrix
L <- matrix(c( 4, 0,
              -1, 1), ncol=2, byrow=TRUE) # Cholesky factor of the ...
(Sigma <- L %*% t(L)) # ... symmetric, positive definite (covariance) matrix Sigma

## Corresponding correlation matrix
P. <- outer(1:2, 1:2, Vectorize(function(r, c)
                          Sigma[r,c]/(sqrt(Sigma[r,r])*sqrt(Sigma[c,c])))) # construct the corresponding correlation matrix
(P  <- cov2cor(Sigma)) # a more elegant solution, see the source of cov2cor()
stopifnot(all.equal(P., P))
## Another option would be as.matrix(Matrix::nearPD(Sigma, corr=TRUE, maxit=1000)$mat)
## which works differently though (by finding a correlation matrix close to the
## given matrix in the Frobenius norm) and thus gives a different answer.


## Decomposing a covariance/correlation matrix

## We frequently need to decompose a covariance matrix Sigma (or correlation
## matrix P) as Sigma = AA^T. To this end there are several possibilities.

## The Cholesky decomposition
A <- t(chol(Sigma)) # the Cholesky factor (lower triangular with real entries > 0)
stopifnot(all.equal(A, L), all.equal(A %*% t(A) , Sigma))

## The eigendecomposition (or spectral decomposition)
eig <- eigen(Sigma) # eigenvalues and eigenvectors
V <- eig$vectors # matrix of eigenvectors
Lambda <- diag(eig$values)
Gamma <- eig$vectors # diagonal matrix with eigenvalues
stopifnot(all.equal(Sigma, V %*% Lambda %*% t(V))) # for real, symmetric matrices
A. <- V %*% sqrt(Lambda) %*% t(V)
Sigma. <- A. %*% t(A.)
stopifnot(all.equal(Sigma, Sigma.))
## However, not that A and A. are different (A. is not a lower triangular matrix!)


### 2 Sampling from the multivariate normal distribution #######################

## Sampling Z (iid N(0,1))
n <- 1000 # sample size
d <- 2 # dimension
set.seed(271) # set a seed (for reproducibility)
Z <- matrix(rnorm(n*d), ncol=d) # sample iid N(0,1) random variates

## Use a linear transformation mu+AZ to get X from N_d(mu, Sigma)
mu <- c(1, 3)
X <- rep(mu, each=n) + t(A %*% t(Z))
## Think about why the rep() and why the t()'s!
## See also Hofert (2013, "On Sampling from the Multivariate t Distribution")

## Plot
plot(X, xlab=expression(X[1]), ylab=expression(X[2]))

## Plot (different mu)
X. <- rep(-mu, each=n) + X # centers the cloud
plot(X., xlab=expression(X[1]), ylab=expression(X[2]))

## Plot (even smaller correlation (larger absolute correlation))
P. <- matrix(c(   1, -0.9,
               -0.9, 1), ncol=d, byrow=TRUE)
Sigma. <- outer(1:d, 1:d, Vectorize(function(r, c)
                P.[r,c] * (sqrt(Sigma[r,r])*sqrt(Sigma[c,c])))) # construct the corresponding Sigma
## Alternatively, change the off-diagonal entry of Sigma while making sure that
## |Sigma.[1,2]| <= sqrt(Sigma.[1,1]) * sqrt(Sigma.[2,2])
A. <- t(chol(Sigma.)) # new Cholesky factor
X. <- rep(mu, each=n) + t(A. %*% t(Z))
plot(X,  xlab=expression(X[1]), ylab=expression(X[2])) # for comparison
plot(X., xlab=expression(X[1]), ylab=expression(X[2])) # larger absolute correlation

## Plot (negative correlation)
Sigma.. <- Sigma * matrix(c( 1, -1,
                            -1,  1), ncol=d, byrow=TRUE)
A.. <- t(chol(Sigma..))
X.. <- rep(mu, each=n) + t(A.. %*% t(Z))
plot(X,   xlab=expression(X[1]), ylab=expression(X[2])) # for comparison
plot(X.., xlab=expression(X[1]), ylab=expression(X[2])) # negative correlation

## Note: This has all been done based on the same Z (and we use the same below)!


### 3 Sampling from the multiv. t distribution (and other normal variance mixtures)

## We use the sample with positive correlation as starting point
X <- X..
A <- A..
Sigma <- A %*% t(A)


### 3.1 A sample from a multivariate t_nu distribution #########################

nu <- 3 # degrees of freedom
W <- 1/rgamma(n, shape=nu/2, rate=nu/2) # mixture variable for a multiv. t_nu distribution
X. <- rep(mu, each=n) + sqrt(W) * t(A %*% t(Z))
xran <- range(X[,1], X.[,1])
yran <- range(X[,2], X.[,2])
plot(X,  xlim=xran, ylim=yran,
     xlab=expression(X[1]), ylab=expression(X[2])) # for comparison
plot(X., xlim=xran, ylim=yran,
     xlab=expression(X[1]), ylab=expression(X[2]))
## For more about sampling from the multivariate t distribution (especially
## what can go wrong!), see Hofert (2013, "On Sampling from the Multivariate t Distribution")


### 3.2 An *uncorrelated* multivariate t_nu distribution vs an independent one

A. <- sqrt(diag(diag(Sigma))) # Cholesky factor of the diagonal matrix containing the variances
X.. <- rep(mu, each=n) + sqrt(W) * t(A. %*% t(Z))
## As a comparison, we also include a sample with *independent* t_3 margins
X... <- rep(mu, each=n) + matrix(rt(n*d, df=nu), ncol=d)
xran <- range(X.[,1], X..[,1], X...[,1])
yran <- range(X.[,2], X..[,2], X...[,2])
plot(X.,  xlim=xran, ylim=yran,
     xlab=expression(X[1]), ylab=expression(X[2])) # for comparison
plot(X.., xlim=xran, ylim=yran,
     xlab=expression(X[1]), ylab=expression(X[2]))
points(X..., col=adjustcolor("royalblue3", alpha.f=0.5)) # => much less mass in the *joint* tails


### 3.3 A 2-point distribution for W (= two 'overlaid' multiv. normals) ########

## Let's change the distribution of W to a 2-point distribution with
## W=1 with prob. 0.7 and W=2 with prob. 0.3. What we will see are two 'overlaid'
## normal distributions (as each of W=1 and W=2 corresponds to one, the former
## even to the original sample, so we colorize the points W=2).
W <- 1+rbinom(n, size=1, prob=0.3)
cols <- rep("black", n)
cols[W==2] <- "royalblue3"
X. <- rep(mu, each=n) + sqrt(W) * t(A %*% t(Z))
xran <- range(X[,1], X.[,1])
yran <- range(X[,2], X.[,2])
plot(X, xlim=xran, ylim=yran,
     xlab=expression(X[1]), ylab=expression(X[2])) # for comparison
plot(X., xlim=xran, ylim=yran, col=cols,
     xlab=expression(X[1]), ylab=expression(X[2]))


### 3.4 A normal *mean*-variance mixture #######################################

## Now let's also 'mix' the mean (replace mu by mu(W)), to get a normal
## mean-variance mixture. The way we do that here is such that we 'split' the
## two clouds of 'overlaid' normal distributions.
mu. <- t(sapply(W, function(w) if(w==1) mu else mu + c(20, 0)))
X.. <- mu. + sqrt(W) * t(A %*% t(Z))
xran <- range(X.[,1], X..[,1])
yran <- range(X.[,2], X..[,2])
plot(X., xlim=xran, ylim=yran, col=cols,
     xlab=expression(X[1]), ylab=expression(X[2]))
plot(X.., xlim=xran, ylim=yran, col=cols,
     xlab=expression(X[1]), ylab=expression(X[2]))
## But, clearly, this is not an elliptical distribution anymore


### 4 Sampling elliptical distributions ########################################

### 4.1 Sample from a multivariate t_nu via elliptical distributions ###########

## Sample from the uniform distribution on the unit sphere (spherical distribution)
set.seed(271)
Z <- matrix(rnorm(n*d), ncol=d)
S <- Z/sqrt(rowSums(Z^2)) # Z/||Z||
par(pty="s")
plot(S, xlab=expression(S[1]), ylab=expression(S[2]))

## Add radial part of a t_nu distribution (spherical distribution)
R <- sqrt(d*rf(n, df1=d, df2=nu))
Y <- R*S
ran <- c(min(Y, -max(Y)), max(Y, -min(Y)))
par(pty="s")
plot(Y, xlim=ran, ylim=ran, xlab=expression(Y[1]), ylab=expression(Y[2]))

## Add scale (elliptical distribution) by computing X = RAS (= ARS)
rho <- 0.5 # correlation
P <- matrix(c(  1, rho,
              rho, 1), ncol=2) # correlation matrix
A <- chol(P) # Cholesky factor
X <- t(A %*% t(Y)) # (n, d) matrix
ran <- c(min(X, -max(X)), max(X, -min(X)))
par(pty="s")
plot(X, xlim=ran, ylim=ran, xlab=expression(X[1]), ylab=expression(X[2]))

## Add location (elliptical distribution) by computing X = mu + RAS
mu <- c(ran[2]-max(X[,1]), ran[2]-max(X[,2])) # push sample towards the upper right corner
X. <- rep(mu, each=n) + X
par(pty="s")
plot(X., xlim=ran, ylim=ran, xlab=expression(X[1]), ylab=expression(X[2]))


### 4.2 A 2-point distribution for R ###########################################

R. <- 1+rbinom(n, size=1, prob=0.3)
X.. <- rep(mu, each=n) + t(A %*% t(R.*S)) # same mu, A, S
par(pty="s")
plot(X.., xlab=expression(X[1]), ylab=expression(X[2]))
