# by Alexander McNeil

# This script should be run immediately after the script 05_SP_block_maxima.R
# It calculates confidence intervals by the profile likelihood method

# The objects fit.a, fit.b and the
##### for CIs
alpha <- 0.05
##############

# The following code is for profile likelihood analysis
# function to maximize profile likelihood
ML.H0 <- function(data,period,level)
{
  sigma0 <- sqrt((6 * var(data))/pi)
  xi0 <- 0.01
  theta <- c(xi0,sigma0)
  parloglik <- function(theta,maxima,k,rlevel)
  {
    implied.mu <- rlevel + theta[2]*(1-(-log(1-1/k))^(-theta[1]))/theta[1]
    -sum(dGEV(maxima,theta[1], implied.mu, abs(theta[2]),log=TRUE))
  }
  optimfit <- optim(theta, fn=parloglik, maxima=data, k=period, rlevel = level)
  par.ests <- optimfit$par
  llmax <- -parloglik(par.ests,data,period,level)
  list(par.ests=par.ests,llmax=llmax)
}


rootfunc <- function(level,data,period,global.max,alpha){
  ML.H0(data,period,level)$llmax - global.max + qchisq(1-alpha,1)/2
}

#########################

global.max.a <- fita$llmax
# CI for 10 year return level
rlevel.10year
ML.H0(M.year,10,rlevel.10year)
ll<-uniroot(rootfunc,lower=2,upper=rlevel.10year,data=M.year,period=10,global.max=global.max.a,alpha=alpha)
ll
ML.H0(M.year,10,ll$root)
uu <-uniroot(rootfunc,lower=rlevel.10year,upper=15,data=M.year,period=10,global.max=global.max.a,alpha=alpha)
uu
ML.H0(M.year,10,uu$root)
c(ll$root,uu$root)

# CI for 50 year return level
rlevel.50year
ML.H0(M.year,50,rlevel.50year)
ll<-uniroot(rootfunc,lower=4,upper=rlevel.50year,data=M.year,period=50,global.max=global.max.a,alpha=alpha)
ll
ML.H0(M.year,50,ll$root)
uu <-uniroot(rootfunc,lower=rlevel.50year,upper=25,data=M.year,period=50,global.max=global.max.a,alpha=alpha)
uu
ML.H0(M.year,50,uu$root)
c(ll$root,uu$root)

# CI for return period in years of Black Monday event
rperiodBM.annual
ML.H0(M.year,rperiodBM.annual,BlackMonday)
ll<-uniroot(rootfunc,lower=10,upper=rperiodBM.annual,data=M.year,level=BlackMonday,global.max=global.max.a,alpha=alpha)
ll
ML.H0(M.year,ll$root,BlackMonday)
# unfortunately, we can't find an upper root of the profile likelihood as following picture shows
pwr <- 1:16
profile <- rep(NA,length(pwr))
for (i in 1:length(pwr)){
  profile[i]<-rootfunc(level=BlackMonday,data=M.year,period =10^pwr[i],global.max=global.max.a,alpha=alpha)
}
plot(10^pwr,profile,type="l",log="x",xlab="return period")
abline(h=0)




global.max.b <- fitb$llmax
# CI for 20 semester return level
rlevel.20semester
ML.H0(M.halfyear,20,rlevel.20semester)
ll<-uniroot(rootfunc,lower=2,upper=rlevel.20semester,data=M.halfyear,period=20,global.max=global.max.b,alpha=alpha)
ll
ML.H0(M.halfyear,20,ll$root)
uu <-uniroot(rootfunc,lower=rlevel.20semester,upper=15,data=M.halfyear,period=20,global.max=global.max.b,alpha=alpha)
uu
ML.H0(M.halfyear,20,uu$root)
c(ll$root,uu$root)

# CI for return period in semesters of Black Monday event
rperiodBM.semester
ML.H0(M.halfyear,rperiodBM.semester,BlackMonday)
ll<-uniroot(rootfunc,lower=10,upper=rperiodBM.semester,data=M.halfyear,level=BlackMonday,global.max=global.max.b,alpha=alpha)
ll
ML.H0(M.halfyear,ll$root,BlackMonday)
uu <-uniroot(rootfunc,lower=rperiodBM.semester,upper=10^8,data=M.halfyear,level=BlackMonday,global.max=global.max.b,alpha=alpha)
uu
ML.H0(M.halfyear,uu$root,BlackMonday)
c(ll$root,uu$root)

