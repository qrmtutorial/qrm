# by Alexander McNeil
library(QRM) # for fitting functions (to change)
library(copula) # for simulation
library(qrmdata) # for data
library(xts)
#### COPULA FITTING

# Fitting copulas to data

data("SP500")
data("FTSE")
# compute returns
SP500.X <- diff(log(SP500))[-1]
FTSE.X <- diff(log(FTSE))[-1]

INDEXES.X <- merge(SP500.X,FTSE.X,all=FALSE)
INDEXES.X <- INDEXES.X['2003-01-01/2012-12-31',]
plot.zoo(INDEXES.X)
pairs(as.zoo(INDEXES.X))

# observe that there are zero values caused by market closures
exclude.rows <- (INDEXES.X[,1]==0) | (INDEXES.X[,2]==0)
sum(exclude.rows)
# remove these
INDEXES.X.d <- INDEXES.X[!exclude.rows,]
# Aggregating by week
INDEXES.X.w <- apply.weekly(INDEXES.X.d,FUN=colSums)
plot.zoo(INDEXES.X.w,type="h")
dim(INDEXES.X.w)



Xdata <- INDEXES.X.w
# Construct pseudo copula data

Udata <- apply(Xdata,2,edf,adjust=1)
plot.zoo(Udata)
plot(Udata[,1],Udata[,2])

# Compare various bivariate models

mod.gauss <- fit.gausscopula(Udata)
mod.gauss
mod.gumbel <- fit.AC(Udata,"gumbel")
mod.gumbel
mod.clayton <- fit.AC(Udata,"clayton")
mod.clayton
mod.t <- fit.tcopula(Udata, startdf =3)
mod.t

c(mod.gauss$ll.max, mod.t$ll.max, mod.gumbel$ll.max, mod.clayton$ll.max)

# Multivariate Fitting with Gauss and t: Simulation study
# We check whether we can recover true parameters

set.seed(117)
G.cop <- ellipCopula("normal", param=0.6, dim=3) # define Gauss copula
Udata <-  rCopula(n=1000, copula=G.cop)
mod.gauss <- fit.gausscopula(Udata)
mod.gauss

# How close does Spearman method get to the optimum?
Pstar <- Spearman(Udata)
sum(dcopula.gauss(Udata,Pstar,log=TRUE))
mod.gauss$ll.max

set.seed(113)
t10.cop <- ellipCopula("t", param=0.6, dim=3, df=10) # define t_10 copula
Udata <-  rCopula(n=1000, copula=t10.cop)
mod.t <- fit.tcopula(Udata)
mod.t
# How close does Kendall method get to the optimum?
mod.t2 <- fit.tcopula(Udata,method="Kendall")
mod.t2$ll.max
mod.t$ll.max
# How good a fit is Gauss copula to t copula data?
mod.gauss <- fit.gausscopula(Udata)
mod.gauss$ll.max

