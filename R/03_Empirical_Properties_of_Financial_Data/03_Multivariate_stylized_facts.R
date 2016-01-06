# by Alexander McNeil
require(xts)
require(qrmdata)

data("FTSE")
data("SMI")

# Index Data 
ftse100.r <- diff(log(FTSE))[-1]
smi.r <- diff(log(SMI))[-1]
FTSE.SMI <- merge(ftse100=ftse100.r,smi=smi.r,all=FALSE)
acf(FTSE.SMI)
acf(abs(FTSE.SMI))

# Multivariate stylized fact: higher correlations in stress periods
require(matrixStats)
X.2d <- FTSE.SMI

X.cor <- apply.monthly(X.2d, FUN=cor)[,2]
X.vols <- apply.monthly(X.2d, FUN=colSds)
plot.zoo(merge(X.cor,X.vols))

fisher <- function(r)
{0.5 * log((1 + r)/(1 - r))}

v1 <- as.numeric(X.vols[,1])
v2 <- as.numeric(X.vols[,2])
rho <- fisher(as.numeric(X.cor))
plot(v1,rho,xlab="Est.vol",ylab="Est.corr")
reg.model = lm(rho ~ v1)
summary(reg.model)
abline(reg.model,col=2)

plot(v2,rho,xlab="Est.vol",ylab="Est.corr")
reg.model = lm(rho ~ v2)
summary(reg.model)
abline(reg.model,col=2)

# EXERCISE: now try with these SWN data
require(mvtnorm)
X.2d <- xts(rmvt(n=dim(FTSE.SMI)[1],df=5,sigma=cor(FTSE.SMI)),time(FTSE.SMI))
X.2d <- xts(rmvnorm(n=dim(FTSE.SMI)[1],sigma=cor(FTSE.SMI)),time(FTSE.SMI))

