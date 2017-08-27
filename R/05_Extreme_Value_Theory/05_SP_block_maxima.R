## By Alexander McNeil

## Note: In comparison to 05_GEV_BMM_SP500.R, this script works with
##       sign-adjusted *percentage* *simple* returns (see 'losses')
##       rather than sign-adjusted log-returns. The results slightly
##       differ.

library(xts)
library(qrmdata)

data("SP500")
summary(SP500)
plot(SP500)

sp500.postcrash <- SP500['1987-10-19/2013-12-31']
sp500.precrash <- SP500['1960-01-01/1987-10-18']
tail(sp500.precrash)
head(sp500.postcrash)

n <- dim(sp500.precrash)[1]
levelstartofweek <- as.numeric(sp500.precrash[n-5])
levelFriday <- as.numeric(sp500.precrash[n])
# compute percentage loss over week
weekreturn <- 100*(levelstartofweek-levelFriday)/levelstartofweek
weekreturn

levelMonday <- as.numeric(sp500.postcrash[1])
# compute percentage loss on Black Monday
BlackMonday <- 100*(levelFriday-levelMonday)/levelFriday
BlackMonday

lreturns <- diff(log(sp500.precrash))[-1]
losses <- -100*(exp(lreturns)-1)
plot(losses)



# compute annual maxim
M.year <- apply.yearly(losses, FUN = max)
M.year
# remove date information
M.year <- as.numeric(M.year)

# there is no half-yearly aggregation function
# so compute quarterly maxima
M.quarterly <- apply.quarterly(losses, FUN = max)
# arrange in a matrix so that each row corresponds to half year
M.quarterly.mat <- matrix(M.quarterly, ncol = 2, byrow = TRUE)
M.quarterly.mat
# compute half-yearly maxima by taking maxima within rows
M.halfyear <- apply(M.quarterly.mat,1,max)
M.halfyear


########  Yearly analysis
library(QRM)
## fit the GEV distribution H_{xi,mu,sigma} to the block maxima
?fit.GEV
fita <- fit.GEV(M.year) # fit a GEV distribution
fita
xi.hat.a <- fita$par.ests[["xi"]]
mu.hat.a <- fita$par.ests[["mu"]]
sigma.hat.a <- fita$par.ests[["sigma"]]
rlevel.10year <- qGEV(1-1/10, xi = xi.hat.a, mu = mu.hat.a, sigma = sigma.hat.a)
rlevel.10year
rlevel.40year <- qGEV(1-1/40, xi = xi.hat.a, mu = mu.hat.a, sigma = sigma.hat.a)
rlevel.40year
rlevel.50year <- qGEV(1-1/50, xi = xi.hat.a, mu = mu.hat.a, sigma = sigma.hat.a)
rlevel.50year
rperiodBM.annual <- 1/(1-pGEV(BlackMonday, xi = xi.hat.a, mu = mu.hat.a, sigma = sigma.hat.a))
rperiodBM.annual

#### Semesterly analysis

fitb <- fit.GEV(M.halfyear) # fit a GEV distribution
fitb
xi.hat.b <- fitb$par.ests[["xi"]]
mu.hat.b <- fitb$par.ests[["mu"]]
sigma.hat.b <- fitb$par.ests[["sigma"]]
rlevel.20semester <- qGEV(1-1/20, xi = xi.hat.b, mu = mu.hat.b, sigma = sigma.hat.b)
rlevel.20semester
rperiodBM.semester <- 1/(1-pGEV(BlackMonday, xi = xi.hat.b, mu = mu.hat.b, sigma = sigma.hat.b))
rperiodBM.semester

# If you are interested in the confidence intervals for these point estimates
# please see 05_SP_rlevels_rperiods_CIs.R

