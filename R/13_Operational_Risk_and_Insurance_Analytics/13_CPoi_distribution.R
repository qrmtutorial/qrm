## By Marius Hofert

## Computing approximate tail probability, VaR and ES estimates of a
## random sum S = \sum_{j=1}^N X_j for independent N ~ Poi(lam), X_j ~ LN(mu, sig^2),
## with various different methods.


### 0 Setup ####################################################################

## Fixed parameters and quantile
lam <- 50 # Poi parameter
mu <- 5 # LN meanlog mean(log(L)) parameter; E(L) = exp(mu + (sig^2)/2)
sig <- 1 # LN sdlog = sd(log(L)) parameter; var(L) = (exp(sig^2)-1)*exp(2*mu + sig^2)
q <- 20000 # quantile
alpha <- 0.99 # VaR, ES confidence level

## Theoretically correct formulas
## Note: - If unknown, these could also be estimated.
##       - For N ~ Poi(lam), EN = var(N) = lam and X ~ F (compound Poi distribution):
##         + ES = EN * EX = lam * EX
##         + var(S) = EN * var(X) + var(N) * (EX)^2
##                  = lam * (var(X) + (EX)^2) = lam * E(X^2)
##         + E(S^2) = var(S) + (ES)^2 = lam * E(X^2) + lam^2 * (EX)^2
##         + skew(S) = E(((S-mean(S))/sd(S))^3) = mu_3 / mu_2^{3/2} (3rd central moment)
##           If k_n denotes the nth cumulant, one can use that mu_2 = k_2, mu_3 = k_3
##           and that for a compound Poi distribution k_n = lam * E(X^n), to see that
##           skew(S) = k_3/k_2^{3/2} = lam * E(X^3) / sqrt((lam * E(X^2))^3)
##                   = E(X^3) / sqrt(lam * E(X^2)^3)
##           Note: One can also work with the 2nd central moment mu_2 = lam * E(X^2)
##       - For X ~ LN(mu, sig^2), E(X^n) = exp(n * mu + (n * sig)^2 / 2)
LN_moment <- function(n, mu, sig) exp(n * mu + (n * sig)^2 / 2)
E.S <- lam * LN_moment(1, mu = mu, sig = sig) # ES
var.S <- lam * LN_moment(2, mu = mu, sig = sig) # var(S)
skew.S <- LN_moment(3, mu = mu, sig = sig) /
    sqrt(lam * (LN_moment(2, mu = mu, sig = sig))^3) # skew(S)


### 1 Normal approximation #####################################################

## Idea: Assume S follows a normal distribution with parameters given by the
##       (true or estimated; here: true) mean and variance of S, so
##       S ~ N(E.S, var.S).

## Tail probability, VaR, ES
(tprob.N <- pnorm(q, mean = E.S, sd = sqrt(var.S), lower.tail = FALSE)) # P(S > q)
(VaR.N <- E.S + sqrt(var.S) * qnorm(alpha)) # VaR_alpha(S)
(ES.N <- E.S + sqrt(var.S) * dnorm(qnorm(alpha)) / (1-alpha)) # ES_alpha(S)


### 2 Translated gamma #########################################################

## Idea: Assume S follows a translated gamma distribution with parameters chosen
##       to match the (true or estimated; here: true) mean, variance and skewness
##       of S.

## Note: S = k + G where G ~ Gamma(shape, rate). Thus
##       1) ES = k + shape/rate
##       2) var(S) = shape/rate^2
##       3) skew(S) = 2/sqrt(shape)
##       and hence
##       3) => shape = (2/skew(S))^2
##       2) => rate = sqrt(shape/var(S))
##       1) => k = ES - shape/rate
shape <- (2/skew.S)^2
rate <- sqrt(shape/var.S)
k <- E.S - shape/rate

## Tail probability, VaR, ES
## Note: Let G ~ Gamma(shape, rate). Then
##       - \bar{F}_S(q) = P(S > q) = P(S-k > q-k) = P(G > q-k) = \bar{F}_G(q-k)
##       - VaR_alpha(S) = F^{-1}_S(alpha) = k + F_G^{-1}(alpha)
##       - ES_alpha(S) = (1/(1-alpha)) int_alpha^1 VaR_u(S) du
##                     = k + (1/(1-alpha)) int_alpha^1 F_G^{-1}(u) du
##                     = k + (1/(1-alpha)) int_{F_G^{-1}(alpha)}^inf x*f_G(x) dx
##         Note that x*f_G(x) = (shape/rate) * <density of Gamma(shape + 1, rate) at x>, so
##                     = k + ((shape/rate) / (1-alpha)) * P(G' > F_G^{-1}(alpha))
##                       for G' ~ Gamma(shape + 1, rate)
(tprob.tg <- pgamma(q-k, shape = shape, rate = rate, lower.tail = FALSE)) # P(S > q) = P(G > q-k)
(VaR.tg <- k + qgamma(alpha, shape = shape, rate = rate)) # VaR_alpha(S)
(ES.tg <- k + ((shape/rate) / (1-alpha)) *
     pgamma(qgamma(alpha, shape = shape, rate = rate),
            shape = shape + 1, rate = rate, lower.tail = FALSE)) # ES_alpha(S)


### 3 Panjer recursion #########################################################

## Idea: - Assume that X_j's are discrete (iid in {0,1,2,...}) with f_k = P(X_1 = k)
##         and that N is discrete (in {0,1,2,...}) with p_k = P(N = k) satisfying
##         p_k = (a+b/k) * p_{k-1} for k >= 1 for some a + b >= 0 (Panjer class;
##         so binomial, Poi or negative binomial); p_0 is such that all p's sum
##         up to 1. Furthermore, all rvs are independent.
##       - Then the compound sum S is discrete with s_{k,n} = P(S = k | N = n).
##         The probability mass function s_{k,n+1} for S given N = n+1 can then
##         be computed recursively via
##               s_{k,n+1} = sum_{j=1}^{k-1} s_{j,n} f_{k-j}
##         (compare the probability of S = k given N = n+1 with the probability
##         of reaching j with n summands and k-j for the remaining summand).
##         Based on this identity, one can show the Panjer recursion. To this end
##         let s_k = P(S = k) and G(z) = E(z^N) be the generating function of N.
##       - Panjer recursion:
##         s_0 <- if(a == 0) G(f_0) else p_0/(1-f_0*a)^(1+b/a)
##         s_k <- (1-f_0*a) * sum_{j=1}^k (a+bj/k) * f_j * s_{k-j}, k = 1,2,...

nsteps <- 50000 # k; should be > floor(q)
stopifnot(nsteps + 1 >= floor(q))

## Discretize the severity distribution
pts <- c(0, 1:(nsteps+1)- 0.5) # => probabilities in (0, 0.5], (1/2, 1.5], (1.5, 2.5], ...
f <- diff(plnorm(pts, meanlog = mu, sdlog = sig)) # length = nsteps + 1

## Panjer recursion; note that for N ~ Poi(lam), we have
## G(z) = exp(lam * (z-1)), a = 0 and b = lam
a <- 0
b <- lam
s <- numeric(nsteps+1) # note: theoretical indices k in {0,1,2,...}; R indices in {1,2,3,...}
s[1] <- exp(lam * (f[1] - 1)) # k = 0; s_0 = G(f_0) here
fct <- 1 - f[1] * a
for(k in 1:nsteps) { # k >= 1
    j <- 1:k
    s[k+1] <- fct * sum((a + b*j/k) * f[j+1] * s[k-j+1]) # s_k
}
## Note:
## - s now contains the mass function of S. This can be used to estimate tail
##   probabilities, VaR_alpha(S) and ES_alpha(S)
## - With nsteps = 25000, the tail probability and VaR_alpha(S) estimates below
##   will be ok, but there is a truncation error for ES_alpha(S) which then
##   leads to a smaller ES_alpha(S) estimate than VaR_alpha(S).

## Tail probability, VaR, ES
S.df <- cumsum(s) # df of S at 0,1,2, ...
(tprob.pr <- 1 - S.df[floor(q)]) # P(S > q) = 1 - P(S <= q)
ii <- which(S.df >= alpha) # indices for which S.df >= alpha
stopifnot(length(ii) > 0) # otherwise S.df < alpha everywhere
(VaR.pr <- min(ii)) # VaR_alpha(S)
(ES.pr <- sum((ii-1) * s[ii]) / (1-alpha)) # ES_alpha(S) = (1/(1-alpha)) * int_{VaR_alpha}^1 x * s(x) dx


### 4 Fast Fourier transform (FFT) #############################################

nsteps <- 2^16

## Discretize the severity distribution
pts <- c(0, 1:(nsteps+1)- 0.5) # => probabilities in (0, 0.5], (1/2, 1.5], (1.5, 2.5], ...
f <- diff(plnorm(pts, meanlog = mu, sdlog = sig)) # length = nsteps + 1

## Fast Fourier transform
f.hat <- fft(f) # Fast Fourier transform
gf.hat <- exp(lam * (f.hat - 1)) # apply the Poisson generating function to the transformed mass function f
s <- Re(fft(gf.hat, inverse = TRUE) / nsteps) # obtain the mass function of S
## Note: We take the real part here, due to possible round-off errors which could
##       introduce an imaginary part.

## Tail probability, VaR, ES
S.df <- cumsum(s) # df of S at 0,1,2, ...
(tprob.fft <- 1 - S.df[floor(q)]) # P(S > q) = 1 - P(S <= q)
ii <- which(S.df >= alpha) # indices for which S.df >= alpha
stopifnot(length(ii) > 0) # otherwise S.df < alpha everywhere
(VaR.fft <- min(ii)) # VaR_alpha(S)
(ES.fft <- sum((ii-1) * s[ii]) / (1-alpha)) # ES_alpha(S) = (1/(1-alpha)) * int_{VaR_alpha}^1 x * s(x) dx
