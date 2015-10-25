## by Alexander McNeil

require(QRM)


load("SP500-1960-2014.RData")
losses <- -returns(SP500[,"Adj.Close"],"simple")*100
SP500 <- cbind(SP500,losses=losses)
sp500.postcrash <- window(SP500,"1987-10-19","2013-12-31")
sp500.precrash <- window(SP500,"1960-01-01","1987-10-18")
tail(sp500.precrash)
head(sp500.postcrash)

n <- dim(sp500.precrash)[1]
levelMonday <- as.numeric(sp500.precrash[n-4,"Open"])
levelFriday <- as.numeric(sp500.precrash[n,"Adj.Close"])
weekreturn <- 100*(levelFriday-levelMonday)/levelMonday
weekreturn
BlackMonday <- as.numeric(sp500.postcrash[1,"losses"])
BlackMonday

date.split <- strsplit(as.character(time(sp500.precrash)), split="-")
head(date.split)
year <- as.numeric(unlist(lapply(date.split, `[[`, 1))) # get years
months <- as.numeric(unlist(lapply(date.split, `[[`, 2))) # get months
halfyear <- rep(1, length(months))
halfyear[months >= 7] <- 2 

blocks <- cbind(year, halfyear) 
colnames(blocks) <- c("year", "halfyear") 
sp500.precrash <- cbind(sp500.precrash,blocks)
head(sp500.precrash)


M.year <- aggregate(losses ~ year, data=sp500.precrash, FUN=max) 
colnames(M.year)[2] <- "annual.max" 
M.halfyear <- aggregate(losses ~ halfyear + year, data=sp500.precrash, FUN=max)
colnames(M.halfyear)[3] <- "semester.max"

########  Yearly analysis

## fit the GEV distribution H_{xi,mu,sigma} to the block maxima
?fit.GEV
fita <- fit.GEV(M.year[,"annual.max"]) # fit a GEV distribution
fita
xi.hat.a <- fita$par.ests[["xi"]]
mu.hat.a <- fita$par.ests[["mu"]]
sigma.hat.a <- fita$par.ests[["sigma"]]
rlevel.10year <- qGEV(1-1/10, xi=xi.hat.a, mu=mu.hat.a, sigma=sigma.hat.a)
rlevel.10year
rlevel.40year <- qGEV(1-1/40, xi=xi.hat.a, mu=mu.hat.a, sigma=sigma.hat.a)
rlevel.40year
rlevel.50year <- qGEV(1-1/50, xi=xi.hat.a, mu=mu.hat.a, sigma=sigma.hat.a)
rlevel.50year
rperiodBM.annual <- 1/(1-pGEV(BlackMonday, xi=xi.hat.a, mu=mu.hat.a, sigma=sigma.hat.a))
rperiodBM.annual

#### Semesterly analysis

fitb <- fit.GEV(M.halfyear[,"semester.max"]) # fit a GEV distribution
fitb
xi.hat.b <- fitb$par.ests[["xi"]]
mu.hat.b <- fitb$par.ests[["mu"]]
sigma.hat.b <- fitb$par.ests[["sigma"]]
rlevel.20semester <- qGEV(1-1/20, xi=xi.hat.b, mu=mu.hat.b, sigma=sigma.hat.b)
rlevel.20semester
rperiodBM.semester <- 1/(1-pGEV(BlackMonday, xi=xi.hat.b, mu=mu.hat.b, sigma=sigma.hat.b))
rperiodBM.semester

# If you are interested in the confidence intervals for these point estimates
# please see 05_SP_rlevels_rperiods_CIs.R