## By Marius Hofert

## Selected VaR superadditivity examples


### 1 Example: Skewed loss distributions #######################################

## Consider two independently defaultable zero-coupon bonds;
## see McNeil et al. (2015, Example 2.25)

## Setup
n <- 1024 # number of alphas
p <- 0.009 # default probability of each of the bonds
alpha <- seq(1e-2, 1-1e-8, length.out = n) # range of confidence levels

## Compute VaR_alpha(L_1 + L_2) and  VaR_alpha(L_1) + VaR_alpha(L_2)
qF <- function(alpha, p)
    ifelse(alpha <= 1-p, -5, 100) # valid for all alpha in (0,1]
qF_sum <- function(alpha, p) # valid for all alpha in (0,1]
    sapply(alpha, function(a) {
        if(a <= (1-p)^2) -10
        else if(a <= 1-p^2) 95
        else 200})
VaR <- 2*qF(alpha, p = p) # VaR_alpha(L_1) + VaR_alpha(L_2)
VaR.sum <- qF_sum(alpha, p = p) # VaR_alpha(L_1+L_2) for L_1, L_2 ind. Exp(1)

## Plot (Note: due to negative values, we have to shift the y values in
## order to plot in log-scale; also, log = "xy" leads to no improvement)
shift <- -min(VaR, VaR.sum) + 1 # + 1 => >= 1
VaR. <- VaR + shift # shift VaR to be > 0
VaR.sum. <- VaR.sum + shift # shift VaR.sum to be > 0
plot(1-alpha, VaR.sum., ylim = range(VaR., VaR.sum.), type = "l", log = "x",
     xlab = expression(1-alpha), ylab = substitute(VaR[alpha]+s, list(s = shift)))
lines(1-alpha, VaR., col = "royalblue3")
legend("bottomleft", lty = 1, col = c("royalblue3", "black"), bty = "n",
       legend = c(expression(VaR[alpha](L[1])+VaR[alpha](L[2])),
                  expression(VaR[alpha](L[1]+L[2]))))
mtext(expression(L[1]*","~L[2]~"two independent skewed loss distributions"),
      side = 4, line = 1, adj = 0)

## We can also detect the superadditive region by noting that
##     F^-_{L_1}(alpha) + F^-_{L_2}(alpha) < F^-_{L_1+L_2}(alpha)
## <=> F_{L_1+L_2}(F^-_{L_1}(alpha) + F^-_{L_2}(alpha)) < alpha
F_sum <- function(x, p)
    sapply(x, function(x.) {
        if(x. < -10) 0
        else if(x. < 95) (1-p)^2
        else if(x. < 200) 1-p^2
        else 1
    })
F_of_F_inverses <- function(alpha, p) F_sum(2*qF(alpha, p = p), p = p)

## Plot
y <- F_of_F_inverses(alpha, p = p)
par(mar = c(5.1, 4.1+2, 4.1, 2.1)) # more space for the y-axis label
plot(alpha, y, type = "l", ylim = range(y, alpha),
     xlab = expression(alpha),
     ylab = expression(F[sum()]*bgroup("(", sum({F[j]^{{phantom(.)%<-%phantom(.)}}}(alpha), j == 1, 2), ")")))
lines(alpha, alpha, lty = 2)
text(0.15, 0.9, label = "VaR subadditive")
text(0.82, 0.05, label = "VaR superadditive")

## With p = 0.1, one sees more:
p <- 0.1
y <- F_of_F_inverses(alpha, p = p)
par(mar = c(5.1, 4.1+2, 4.1, 2.1)) # more space for the y-axis label
plot(alpha, y, type = "l", ylim = range(y, alpha),
     xlab = expression(alpha),
     ylab = expression(F[sum()]*bgroup("(", sum({F[j]^{{phantom(.)%<-%phantom(.)}}}(alpha), j == 1, 2), ")")))
lines(alpha, alpha, lty = 2)
text(0.15, 0.9, label = "VaR subadditive")
text(0.82, 0.05, label = "VaR superadditive")
## Note: Hofert and McNeil (2014, "Subadditivity of Value-at-Risk for Bernoulli
##       random variables") showed that VaR_alpha is always subadditive for
##       alpha > 1-p.


### 2 Example: Independent, light-tailed losses and small to moderate alpha ####

## Consider two independent Exp(1) random variables and small to moderate alpha;
## see McNeil et al. (2015, Example 7.30)

## Setup
n <- 1024 # number of alphas
alpha <- seq(1e-2, 1-1e-8, length.out = n) # range of confidence levels

## Compute VaR_alpha(L_1 + L_2) and  VaR_alpha(L_1) + VaR_alpha(L_2)
VaR <- 2*qexp(alpha) # VaR_alpha(L_1) + VaR_alpha(L_2)
VaR.sum  <- qgamma(alpha, shape = 2) # VaR_alpha(L_1 + L_2) for L_1, L_2 ind. Exp(1)
alpha.super <- max(alpha[VaR.sum > VaR]) # ~ 0.7145 (max. alpha for which VaR superadd.)

## Plot (Note: no improvement by plotting in 1-alpha and with log = "xy")
plot(alpha, VaR.sum, ylim = range(VaR, VaR.sum), type = "l", log = "y",
     xlab = expression(alpha), ylab = expression(VaR[alpha]))
lines(alpha, VaR, col = "royalblue3")
abline(v = alpha.super, lty = 2)
legend("topleft", lty = 1, col = c("black", "royalblue3"), bty = "n",
       legend = c(expression(VaR[alpha](L[1]+L[2])),
                  expression(VaR[alpha](L[1])+VaR[alpha](L[2]))))
text(alpha.super-0.16, 0.025, labels = "Superadditivity\nregion")
text(alpha.super+0.16, 0.025, labels = "Subadditivity\nregion")
mtext(expression("Independent"~~L[1]*", "*L[2]%~%"Exp(1)"), side = 4, line = 1, adj = 0)

