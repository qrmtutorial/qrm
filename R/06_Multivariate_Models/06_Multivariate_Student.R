# Simulate multivariate Student
S <- equicorr(d = 3, rho = 0.7)
d <- 3
nu <- 10
data <- rmt(2000, Sigma = S, df=nu)
pairs(data)

# Test of marginal normality
apply(data,2,shapiro.test)

# Tests of joint normality
MardiaTest(data)
jointnormalTest(data)

# QQplots against Student t
qplotfunc <- function(data){
  qqplot(qt(ppoints(data),df=nu),data)
  qqline(data,dist=function(p){qt(p,df=nu)})
}
par(mfrow=c(1,3))
apply(data,2,qplotfunc)
par(mfrow=c(1,1))

# QQplot of Mahalanobis distances (scaled by dimension d)
# against an F(d,nu) distribution
Ddata <- mahalanobis(data,center=FALSE,cov=S)/d
qqplot(qf(ppoints(Ddata),df1=d,df2=nu),Ddata)
qqline(Ddata,dist=function(p){qf(p,df1=d,df2=nu)})

# Estimate multivariate Student
# Add library
require(QRM)

# Dow Jones Data
data(DJ)
r <- returns(DJ)
stocks <- c("AXP","EK","BA","C","KO","MSFT","HWP","INTC","JPM","DIS")
ss <- window(r[, stocks], "1993-01-01", "2000-12-31")
plot(ss)

# Fit normal
mod.norm = fit.norm(ss)
mod.norm

# Fit Student t
mod.t = fit.mst(ss,method="BFGS")
mod.t
names(mod.t)
mod.t$df

# Compare
mod.t$ll.max - mod.norm$ll.max

