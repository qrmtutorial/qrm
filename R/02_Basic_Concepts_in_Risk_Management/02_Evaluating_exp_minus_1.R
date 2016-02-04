## By Marius Hofert

## This illustrates a potential numerical issue when evaluating exp(x)-1
## appearing in the loss operator in a standard stock portfolio;
## see McNeil et al. (2015, Example 2.1).


### Setup ######################################################################

library(sfsmisc) # for eaxis()


### Investigating exp(x)-1 vs expm1(x) #########################################

x <- 10^seq(-26, 0, by=0.25) # sequence of x values
y <- exp(x)-1 # numerically critical for |x| ~= 0
y. <- expm1(x) # numerically more stable

## Plot exp(x)-1 vs expm1(x) near 0
plot(x, y., type="l", log="xy", # expm1(x)
     col="royalblue3", ann=FALSE, xaxt="n", yaxt="n")
mtext(1, text="x", line=2) # x-axis label
eaxis(1, f.smalltcl=2/5) # x-axis
eaxis(2, f.smalltcl=2/5) # y-axis
lines(x, y, col="firebrick", lwd=1.6) # exp(x)-1
legend("bottomright", legend=c("exp(x)-1", "expm1(x)"), bty="n",
       lty=c(1,1), lwd=c(1.6,1), col=c("firebrick", "royalblue3")) # legend


### Illustrating how expm1(x) works ############################################

## See the source under ./src/nmath/expm1.c

plot(x, y., type="l", lwd=2, log="xy", col="royalblue3",
     xlab="x", ylab="expm1(x)", xaxt="n", yaxt="n") # expm1(x)
eaxis(1, f.smalltcl=2/5) # x-axis
eaxis(2, f.smalltcl=2/5) # y-axis
ypos <- 1e-25 # y-axis height for labels
## Right-most part (non-critical)
text(x=3, y=ypos, "exp(x)-1", srt=90)
## Middle-right part (improvement via Newton step)
text(x=10^(-(8+0.1681302)/2), y=ypos, "y=exp(x)-1\n+ Newton step")
abline(v=0.697, lty=2)
## Middle-left part (improvement via Taylor approximation + Newton step)
text(x=10^(-(15.65356+8)/2), y=ypos, "y=(1+x/2)x\n+ Newton step")
abline(v=1e-8, lty=2)
## Left-most part (improvement via the identity, so again Taylor)
text(x=10^(-(27+15.65356)/2), y=ypos, "x")
abline(v=.Machine$double.eps, lty=2) # smallest floating-point number x > 0 s.t. 1 + x != 1
