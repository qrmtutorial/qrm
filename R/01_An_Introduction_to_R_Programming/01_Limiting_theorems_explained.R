## By Marius Hofert

## R is strong on quickly visualizing concepts. We do this here by numerically
## verifying major limiting theorems in Probability (SLLN, CLT) and
## Extreme Value Theory (Gnedenko's Theorem, Pickands--Balkema--de Haan).
## Note that we use the same generated data throughout!


### Generate the data ##########################################################

require(qrmtools)

## Data from a Par(theta) distribution
n <- 50000 # sample size = number of i.i.d. random variables
th <- 3 # parameter theta
set.seed(271) # set seed for reproducibility
X <- rPar(n, theta=th) # generate data

## Build blocks of data (for CLT, Gnedenko)
m <- 500 # number of blocks
X. <- split(X, f=rep(1:m, each=floor(n/m))) # split data into blocks


### 1 Strong Law of Large Numbers (SLLN) #######################################

## Building cumulative averages (X_1/1, (X_1+X_2)/2, (X_1+X_2+X_3)/3,...)
stopifnot(th > 1)
Xn <- cumsum(X)/(1:n)

## Plot (this one path of the stochastic process (bar{X}_n)_{n=1}^{\infty})
plot(1:n, Xn, type="l", log="x", ylab="",
     xlab=expression("Number"~italic(n)~"of i.i.d. random variables"~italic((X[i])[i==1]^n)),
     main=substitute(bold("Strong Law of Large Numbers for Par("*th.*") data"), list(th.=th)))
mu <- 1/(th-1) # mean
abline(h=mu, col="blue")
legend("bottomright", lty=c(1,1), col=c("black", "blue"), bty="n", y.intersp=1.2,
       legend=c(expression(italic((bar(X)[n])[n])),
                substitute("true mean"~mu==mu., list(mu.=mu))))


### 2 Central Limit Theorem (CLT) ##############################################

## Standardize blocked data via sqrt(n) * (bar{X}_n - mu) / sigma
stopifnot(th > 2)
mu <- 1/(th-1) # mean
sig2 <- 2/((th-2)*(th-1)^2) # variance
Z <- sapply(X., function(x) sqrt(length(x))*(mean(x)-mu)/sqrt(sig2)) # standardize

## Histogram
hist(Z, probability=TRUE, ylim=c(0, dnorm(0)),
     main=substitute(bold("Central Limit Theorem for Par("*th.*") data"),
     list(th.=th)), xlab=expression("Realizations of"~sqrt(n)*(bar(X)[n]-mu)/sigma))
curve(dnorm, from=min(Z), to=max(Z), add=TRUE, col="blue") # overlay N(0,1) density
box()


### 3 Gnedenko's Theorem #######################################################

## \bar{F}(x) = x^{-\theta} L(x) for L(x) = (1+1/x)^{-\theta} (slowly varying at
## infinity). By Gnedenko's Theorem (1943), such F is in MDA(H_{1/\theta}) for
## all \theta>0 and the normalizing sequences are c_n=F^-(1-1/n) and d_n=0 for
## all n so that (M_n-d_n)/c_n is approximately distributed as H_{1/\theta}.
## Let's check that.

## Standardize blocked data via
M <- sapply(X., function(x) (max(x) - 0) / qPar(1-1/length(x), theta=th)) # standardize

## Histogram
hist(M, probability=TRUE,
     main=substitute(bold("Gnedenko's Theorem for Par("*th.*") data"),
     list(th.=th)), xlab=expression("Realizations of"~(M[n]-d[n])/c[n]~~
                                    "(for"~c[n]=={F^{-1}}(1-1/n)*";"~~d[n]==0*")"))
curve(dGEV(x, xi=1/th), from=0, to=max(M), add=TRUE, col="blue")
box()

## Q-Q plot
M. <- sort(M)
qGEV. <- qGEV(ppoints(length(M.)), xi=1/th)
plot(qGEV., M., xlab="Theoretical quantiles", ylab="Sample quantiles",
     main=substitute(bold("Gnedenko's Theorem for Par("*th.*") data"),
          list(th.=th)))
qqline(y=M., distribution=function(p) qGEV(p, xi=1/th))


### 4 Pickands--Balkema--de Haan (1974/1975) ###################################

## For sufficiently large thresholds u, excesses over u follow a GPD(\xi,\beta)
## distribution if and only if F is in MDA(H_\xi). If F is GPD(\xi,\beta)
## one can show that the excess distribution F_u over u is GPD(\xi,\beta+\xi*u).
## Since that's the case for Pareto distributions (with \xi=1/\theta, \beta=\xi),
## excesses over u should follow a GPD(1/\theta,(1/\theta)*(1+u)) distribution.
## Let's check that.

## Determine excesses
u <- quantile(X, 0.9) # use the 90% quantile (rule of thumb)
Y <- X[X>u]-u # excesses over u

## Histogram (note: true density peaks near 0)
hist(Y, probability=TRUE,
     main=substitute(bold("Pickands-Balkema-de Haan Theorem for Par("*th.*") data"),
          list(th.=th)), xlab="Realizations of excesses Y over the threshold u (90% quantile)")
curve(dGPD(x, xi=1/th, beta=1/th), from=0, to=max(Y), add=TRUE, col="blue")
box()

## Q-Q plot
Y. <- sort(Y)
qGPD. <- qGPD(ppoints(length(Y.)), xi=1/th, beta=1/th)
plot(qGPD., Y., xlab="Theoretical quantiles", ylab="Sample quantiles",
     main=substitute(bold("Pickands-Balkema-de Haan Theorem for Par("*th.*") data"),
          list(th.=th)))
qqline(y=Y., distribution=function(p) qGPD(p, xi=1/th, beta=1/th))
