## By Alexander McNeil

library(QRM)
library(lme4)

data(spdata.raw)
spdata.raw
attach(spdata.raw)

# Calculate default rates

BdefaultRate <- Bdefaults/Bobligors
BBdefaultRate <- BBdefaults/BBobligors
BBBdefaultRate <- BBBdefaults/BBBobligors
AdefaultRate <- Adefaults/Aobligors
CCCdefaultRate <- CCCdefaults/CCCobligors

# Plot default rates

year <- 1981:2000
plot(year,CCCdefaultRate,xlab="Year",ylab="Rate",type="l")
lines(year,BdefaultRate,col=2)
lines(year,BBdefaultRate,col=3)
lines(year,BBBdefaultRate,col=4)
lines(year,AdefaultRate,col=5)

# simple moment estimator
pY <- momest(Bdefaults,Bobligors)
rhoY <- (pY[2]-pY[1]^2)/(pY[1]-pY[1]^2)

# Fit binomial model followed by one factor model
mod0 <- fit.binomial(Bdefaults, Bobligors)
mod1 <- fit.binomialProbitnorm(Bdefaults, Bobligors)
c(mod0$maxloglik, mod1$maxloglik)

# Advanced
# Gives table 8.8 in first edition of book

# Change format of data for glmm

defaultmatrix <- cbind(Adefaults,BBBdefaults,BBdefaults,Bdefaults,CCCdefaults)
firmmatrix <- cbind(Aobligors,BBBobligors,BBobligors,Bobligors,CCCobligors)

defaults <- as.vector(t(defaultmatrix))
firms <- as.vector(t(firmmatrix))
year <- rep(year,rep(5,20))
rating <- factor(rep(c("A","BBB","BB","B","CCC"),20),levels=c("A","BBB","BB","B","CCC"),ordered=TRUE)

tail(data.frame(defaults,firms,year,rating),n=10)

# Fit glmm

(mod <- glmer(cbind(defaults,firms-defaults) ~ -1 + rating + (1|year),family=binomial(probit)))

# Compute implied PDs and asset correlation

sigma <- mod@theta
mu <- mod@beta

beta <- sigma^2/(1+sigma^2)
# Asset correlation. Much smaller than the numbers coming out of equity analysis

PD <- pnorm(mu*sqrt(1-beta))

