# by Alexander McNeil
require(QRM)
load("INDEXES-2000-2012.RData")

# Index Data 

plot(INDEXES0012)
X.INDEXES <- returns(INDEXES0012)



# Multivariate stylized fact: higher correlations in stress periods
FTSE.SMI <- X.INDEXES[,c(1,3)]
(n <- dim(FTSE.SMI)[1])

X.2d <- FTSE.SMI

X.prod <- X.2d[,1]*X.2d[,2]
by <- unique(timeLastDayInMonth(time(X.INDEXES)))
counts <- aggregate(X.2d, by, length)
counts
X.vols <- aggregate(X.2d, by, sd)
X.means <- aggregate(X.2d, by, mean)
Xprod.mean <- aggregate(X.prod, by, mean)
X.cov <- Xprod.mean - X.means[,1]*X.means[,2]
X.cor <- X.cov/(X.vols[,1]*X.vols[,2])
plot(X.cor,ylim=c(-1,1))

fisher <- function(r)
{0.5 * log((1 + r)/(1 - r))}

plot(X.vols[,1],fisher(X.cor),xlab="Est.vol",ylab="Est.corr")
reg.model = lm(fisher(X.cor)~X.vols[,1])
summary(reg.model)
abline(reg.model,col=2)

plot(X.vols[,2],fisher(X.cor),xlab="Est.vol",ylab="Est.corr")
reg.model = lm(fisher(X.cor)~X.vols[,2])
summary(reg.model)
abline(reg.model,col=2)

# EXERCISE: now try with these SWN data
X.2d <- timeSeries(rmt(n=n,df=5,Sigma=cor(FTSE.SMI)),time(FTSE.SMI))
X.2d <- timeSeries(rmnorm(n=n,Sigma=cor(FTSE.SMI)),time(FTSE.SMI))

