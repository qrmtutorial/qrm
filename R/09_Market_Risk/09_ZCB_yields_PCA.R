### by Alexander J. McNeil
require(xts)
require(qrmdata)
require(termstrc)
data(ZCB_CA)



# Create a "zero-yields" object with 10 years of data now
ZCB10yr <- ZCB_CA['2002-01-02/2011-12-30']

yields <- 100*as.matrix(ZCB10yr)
dates <- time(ZCB10yr)
maturities = (1:120)*0.25
datazeroyields = zeroyields(maturities, yields, dates)
datazeroyields



# PCA Analysis of daily changes in yields
# We first difference the yields

Xyields = apply(yields,2,diff)
head(Xyields)
Xyields.series= xts(Xyields,dates[-1])
plot.zoo(Xyields.series[,1:10])

# Carry out PCA analysis
tmp = princomp(Xyields)
(results <- summary(tmp))

loads <- loadings(tmp)
loads[,1:3]

plot((1:120)/4,loads[,1],ylim=range(loads[,1:3]),type="l",xlab="Maturity",ylab="Loading")
lines((1:120)/4,loads[,2],col=2)
lines((1:120)/4,loads[,3],col=3)
abline(h=0,lty=2)
legend(x=c(15,30), y=c(0.2,0.45),legend=c("PC 1","PC 2","PC 3"),col=1:3,lty=c(1,1,1))


plot(tmp)
scores = tmp$scores
factors = scores[,1:3]

factors.series = xts(factors,dates[-1])
plotfunc <- function(x,...)
{	lines(x,...)
	abline(v=as.Date("2008-09-15"),col=4)}

plot.zoo(factors.series,yax.flip=TRUE,panel=plotfunc)

cor(factors)
acf(factors)
acf(abs(factors))





