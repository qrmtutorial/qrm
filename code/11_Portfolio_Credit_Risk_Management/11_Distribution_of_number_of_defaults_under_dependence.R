## By Alexander J. McNeil and Marius Hofert

## Computing the probability mass function (pmf) of the number of defaults
## under exchangeable one-factor Bernoulli mixture models.

## Recall from MFE (2015, 11.2.1--11.2.2):
## 1) Bernoulli mixture models (in general):
##    - jth component: P(Y_j = y_j | Psi = psi) = p_j(psi)^{y_j} (1-p_j(psi))^{1-y_j}
##      where Psi is a p-dimensional factor vector, j in {1,...,m}
##    - Under conditional independence:
##      P(Y = y | Psi = psi) = prod_{j=1}^m p_j(psi)^{y_j} (1-p_j(psi))^{1-y_j},
##      where Y = (Y_1,...,Y_m) is the vector of default indicators.
## 2) Exchangeable one-factor case:
##    - Psi now a random variable ('one-factor')
##    - all functions p_j equal to p ('exchangeable')
##    - define q = p_1(psi) and Q = p_1(Psi) ~ G
##    => P(Y = y | Q = q) = prod_{j=1}^m q^{y_j} (1-q)^{1-y_j}
##                        = q^{sum_{j=1}^m y_j} (1-q)^{m-sum_{j=1}^m y_j}
## 3) Implied default probabilities:
##    - jth component: pi_1 = P(Y_j = 1) = 1 * P(Y_j = 1) + 0 * P(Y_j = 0)
##                          = E(Y_j) = E(E(Y_j | Q)) = E(Q) (as (Y_j | Q) ~ B(1, Q))
##    - k components jointly: pi_k = P(Y_{j_1} = ... = Y_{j_k} = 1)
##                                 = E(Y_{j_1} * ... * Y_{j_k})
##                                 = E(E(Y_{j_1} * ... * Y_{j_k} | Q))
##                                 = E(E(Y_{j_1} | Q) * ... * E(Y_{j_k} | Q))
##                                 = E(Q^k) since E(Y_{j_k} | Q) ~ B(1, Q) with mean Q
## 4) Implied default correlation:
##    rho_Y = Cor(Y_i, Y_j) = (E(Y_i * Y_j) - E(Y_i) E(Y_j)) /
##                            (sqrt{E(Y_i)(1-E(Y_i))} * sqrt{E(Y_j)(1-E(Y_j))})
##    by 3) => rho_Y = (pi_2 - pi_1^2) / (pi_1 (1 - pi_1)) = (E(Q^2) - E(Q)^2) / (E(Q)(1-E(Q)))
## 5) For a df G of Q, match its first two moments with the desired individual
##    default probability and the default correlation => parameters of G.
##    - pi_1 = E(Q); see 3)
##    - rho_Y = (E(Q^2) - E(Q)^2) / (E(Q)(1-E(Q))); see 4)
## 6) Now consider the number of defaults M = sum_{j=1}^m Y_j:
##    - P(M = k | Q = q) = choose(m, k) q^k (1-q)^{m-k}
##      => (M | Q = q) ~ B(m, q)
##      => P(M = k) = int_0^1 P(M = k | Q = q) dG(q)
##                  = choose(m, k) int_0^1 q^k (1-q)^{m-k} dG(q)
##    - Note: If Q = q is constant, then P(M = k) = choose(m, k) q^k (1-q)^{m-k}
##            so M ~ B(m, q).


### Setup ######################################################################

library(copula)

## Input parameters
m <- 1000 # number m of loans
pi. <- 0.01 # individual default probability pi (of each loan)
rho. <- 0.005 # default correlation rho_Y = Cor(Y_i, Y_j) (Y_. = default indicator)

## Number k of defaulted loans we consider (we plot P(M = k) for those k)
k <- 0:60


### 1 Compute default probabilities in the independent case ####################

## Here M ~ B(m, pi.); see 3) above
pmf.indep <- dbinom(k, size = m, prob = pi.) # pmf


### 2 Compute default probabilities under the beta-binomial model ##############

## We consider Q ~ B(a,b). M then follows a so-called 'beta-binomial' distribution.
## MFE (2015, Example 11.7) provides pi and rho_Y as functions
## of a, b; as in 6) above. Solving for a, b leads to:
(a <- pi. * (1/rho. - 1))
(b <- a*(1/pi. - 1))
stopifnot(a > 0, b > 0)

## Define the beta-binomial pmf of M; see MFE (2015, Example 11.7).
dbetaBinom <- function(x, size, a, b) # x = k, size = m
    exp(lchoose(size, x) + lbeta(a + k, b + size - k) - lbeta(a, b)) # see (11.16)

## Compute the pmf P(M = k) for the given k
pmf.betaBinom <- dbetaBinom(k, size = m, a = a, b = b)

## Plot P(M = k) for the given k including the independence case
ylim <- range(pmf.betaBinom, pmf.indep)
plot(k, pmf.betaBinom, type = "l", ylim = ylim, col = "maroon3",
     xlab = "Number of losses", ylab = "")
lines(k, pmf.indep) # add pmf under independence
legend("topright", bty = "n", lty = 1, col = c("maroon3", "black"),
       legend = c("Probability under beta-binomial", "Probability under independence"))

## Plot their quotient
plot(k, pmf.betaBinom/pmf.indep, type = "l", log = "y", xlab = "Number of losses", ylab = "")
abline(h = 1, lty = 2)
legend("topleft", bty = "n", lty = c(1, 2), legend = c(expression(frac("Probability under beta-binomial",
                                                                       "Probability under independence")),
                                                       "Reference line y = 1"))


### 3 Compute default probabilities under the probit-normal model ##############

## We consider Q (= p_1(Psi)) ~ Phi(mu + sig * Psi) (*) for Psi ~ N(0,1) with df
## G(q) = Phi((Phi^{-1}(q) - mu) / sig). M then follows a so-called
## 'probit-normal' distribution. To find mu and sig, we need more work here.

## First we need to find the copula parameter (here: Gauss copula parameter beta)
## such that C(pi, pi) = pi_2; see MFE (2015, (11.5)).
root <- function(x, p, p2) # x = Gauss copula parameter; p = pi, p2 = pi_2
    p2 - pCopula(rep(p, 2), copula = normalCopula(x))
pi2 <- pi.^2 + rho. * (pi. * (1-pi.)) # compute pi_2 from rho_Y and pi; see (11.3)
beta <- uniroot(root, interval = 0:1, p = pi., p2 = pi2)$root

## Compute the parameters of the probit-normal model. From (*), we have
## p_1(psi) = Phi(mu + sig * psi) and from MFE (2015, (11.21)) we
## know that p_1(psi) = Phi((Phi^{-1}(pi) + sqrt(beta) * psi) / sqrt(1-beta)).
## Matching the two leads to:
mu  <- qnorm(pi.) / sqrt(1 - beta)
sig <- sqrt(beta/(1 - beta))

## Define the probit-normal pmf of M; see MFE (2015, (11.14)).
dprobitNorm <- function(x, size, mu, sigma) # x = k, size = m
{
    ## Define integrand choose(m, k) q^k (1-q)^{m-k} phi((Phi^{-1}(q) - mu) / sigma) /
    ## (sigma * phi(Phi^{-1}(q))), that is, choose(m, k) q^k (1-q)^{m-k} times
    ## the density of G(q) as above.
    integrand <- function(q, k)
        exp(lchoose(size, k) + k * log(q) + (size-k) * log(1-q) +
            dnorm((qnorm(q)-mu)/sigma, log = TRUE) -
            dnorm(qnorm(q), log = TRUE) - log(sigma))
    ## Integrate to obtain P(M = k); see MFE (2015, (11.14))
    sapply(1:length(x), function(k)
        integrate(integrand, lower = 0, upper = 1, k = k)$value)
}

## Compute the pmf P(M = k) for the given k
pmf.probitNorm <- dprobitNorm(k, size = m, mu = mu, sigma = sig)

## Plot P(M = k) for the given k including the independence case
ylim <- range(pmf.probitNorm, pmf.indep)
plot(k, pmf.probitNorm, type = "l", ylim = ylim, col = "royalblue3",
     xlab = "Number of losses", ylab = "")
lines(k, pmf.betaBinom, col = "maroon3") # add pmf under beta-binomial
lines(k, pmf.indep) # add pmf under independence
legend("topright", bty = "n", lty = 1, col = c("maroon3", "royalblue3", "black"),
       legend = c("Probability under beta-binomial", "Probability under probit-normal",
                  "Probability under independence"))

