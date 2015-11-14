library(QRM) # for fitting functions (to change)
library(copula) # for simulation

#### COPULA FITTING

# Fitting copulas to data

load("INDEXES-2000-2012.RData")
plot(INDEXES0012)
Xdata <- returns(INDEXES0012[,c(1,3)])
plot(Xdata)
plot(series(Xdata))
dim(Xdata)

# Construct pseudo copula data

Udata <- apply(Xdata,2,edf,adjust=1)
plot(series(Udata))
# observe that there are tied values caused by holidays
exclude.rows <- (Xdata[,1]==0) | (Xdata[,2]==0)
sum(exclude.rows)
Xdata <- Xdata[!exclude.rows,]
Udata <- apply(Xdata,2,edf,adjust=1)
plot(series(Udata))


# Compare various bivariate models

mod.gauss <- fit.gausscopula(Udata)
mod.gauss
mod.gumbel <- fit.AC(Udata,"gumbel")
mod.gumbel
mod.clayton <- fit.AC(Udata,"clayton")
mod.clayton
mod.t <- fit.tcopula(Udata)
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
Udata <- rcopula.t(1000,df=10,Sigma=P)
mod.t <- fit.tcopula(Udata)
mod.t
# How close does Kendall method get to the optimum?
mod.t2 <- fit.tcopula(Udata,method="Kendall")
mod.t2$ll.max
mod.t$ll.max
# How good a fit is Gauss copula to t copula data?
mod.gauss <- fit.gausscopula(Udata)
mod.gauss$ll.max

