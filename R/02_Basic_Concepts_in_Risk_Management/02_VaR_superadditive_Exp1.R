## By Marius Hofert

## VaR superadditivity for two independent Exp(1) random variables


### Setup ######################################################################

n <- 1024 # number of alphas
alpha <- seq(1e-2, 1-1e-8, length.out=n) # range of confidence levels
VaR  <- qgamma(alpha, shape=2) # VaR_alpha(L_1+L_2), L_j ind Exp(1)
VaRLL <- 2*qexp(alpha) # VaR_alpha(L_1) + VaR_alpha(L_2)
(alpha.super <- max(which(VaR > VaRLL))/n) # ~ 0.7119


### Plot #######################################################################

## Plot (note: no improvement by plotting in 1-alpha and with log="xy")
plot(alpha, VaR, ylim=range(VaR, VaRLL), type="l", log="y",
     xlab=expression(alpha), ylab=expression(VaR[alpha]))
lines(alpha, VaRLL, col="royalblue3")
abline(v=alpha.super, lty=2)
legend("topleft", inset=0.04, lty=1, col=c("black", "royalblue3"), bty="n",
       legend=c(expression(VaR[alpha](L[1]+L[2])),
                expression(VaR[alpha](L[1])+VaR[alpha](L[2]))))
text(alpha.super-0.16, 0.025, labels=expression(VaR[alpha]~"superadditive"))
text(alpha.super+0.16, 0.025, labels=expression(VaR[alpha]~"subadditive"))
mtext(expression("Independent"~~L[1]*", "*L[2]%~%"Exp(1)"), side=4, line=1, adj=0)

