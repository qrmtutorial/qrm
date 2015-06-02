### By Marius Hofert


### Functions ##################################################################

## Random angle for Student t distribution
Theta <- function(nu, n.MC) {
    W1 <- nu/rchisq(n.MC, df=nu)
    W2 <- nu/rchisq(n.MC, df=nu)
    W3 <- nu/rchisq(n.MC, df=nu)
    2*asin(sqrt(W1 * W2 / ((W1+W3) * (W2+W3))))
}

## Spearman's rho for Student t distribution
rho_t <- function(rho, nu, n.MC) {
    stopifnot(-1 <= rho, rho <= 1, nu > 0, n.MC > 0)
    6*colMeans(asin(outer(sin(Theta(nu, n.MC=n.MC)/2), rho)))/pi # length(rho)
}


### Spearman's rho versus linear correlation for Student t distributions #######

set.seed(271)

## Generate data
nu <- 10^seq(-1, 1)
rho <- seq(-1, 1, by=0.01)
rho.t <- sapply(nu, function(nu.) rho_t(rho, nu=nu., n.MC=1e4)) # (length(rho), length(nu)-matrix

## Plot
plot(rho, rho.t[,1], type="l",
     xlab=expression("Pearson's"~~rho), ylab=expression("Spearman's"~~rho))
lines(rho, rho.t[,2], lty=2)
lines(rho, rho.t[,3], lty=3)
abline(0, 1, lty=4)
legend("bottomright", bty="n", lty=1:3,
       legend=as.expression(c(substitute(nu==nu., list(nu.=nu[1])),
                              substitute(nu==nu., list(nu.=nu[2])),
                              substitute(nu==nu., list(nu.=nu[3])))))


### Convergence of the density of the random angle Theta to the point mass at
### pi/3 corresponding to a bivariate normal

## Generate data
nu <- 2^seq(-1, 5, length.out=4)
Theta.nu <- sapply(nu, Theta, n.MC=1e4) # (n.MC, nu)-matrix
dens.estim <- apply(Theta.nu, 2, density) # length(nu)-list; density estimates

## Plot
ylim <- range(sapply(dens.estim, function(x) x[["y"]])) # determine the y-range
plot(dens.estim[[1]], ylim=ylim, main="", xlab="x")
for(i in seq_along(nu))
    lines(dens.estim[[i]], lty=i)
abline(v=pi/3, col="royalblue3")
legend("topright", bty="n", lty=seq_along(nu),
       legend=as.expression(c(substitute(nu==nu., list(nu.=nu[1])),
                              substitute(nu==nu., list(nu.=nu[2])),
                              substitute(nu==nu., list(nu.=nu[3])),
                              substitute(nu==nu., list(nu.=nu[4])))))


