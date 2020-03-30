## By Alexander McNeil and Marius Hofert

## Estimating parameters of portfolio models from default data


### Setup ######################################################################

library(QRM) # for momest(), fit.binomial(), fit.binomialProbitnorm()
library(qrmdata)
library(lme4) # for glmer()


### Data preparation ###########################################################

## Load
data(SP_defaults)
str(SP_defaults)

## Pick out defaults, firms, and ratings
defaults <- SP_defaults[,"Defaults",]
firms <- SP_defaults[,"Obligors",]
rating <- colnames(firms)

## Calculate default rates
(defRate <- apply(SP_defaults, c(1, 3), function(x) x["Defaults"]/x["Obligors"]))
nrate <- ncol(defRate)

## Plot default rates
years <- as.numeric(format(as.Date(dimnames(SP_defaults)$Time), "%Y"))
pal <- colorRampPalette(c("black", "royalblue3", "darkorange2", "maroon3"), space = "Lab")
cols <- pal(nrate)
plot(NA, xlim = range(years), ylim = range(defRate),
     xlab = "Year", ylab = "Default rate")
for(r in 1:nrate)
    lines(years, defRate[,r], type = "l", col = cols[r])
legend("topleft", bty = "n", lty = rep(1, nrate), col = cols, legend = rating)


### Fitting ####################################################################

## Simple moment estimator (for 'B')
pY <- momest(defaults[,"B"], firms[,"B"])
rhoY <- (pY[2]-pY[1]^2) / (pY[1]-pY[1]^2) # (11.3)

## Binomial model
mod0 <- fit.binomial(defaults[,"B"], firms[,"B"])

## One-factor model (p. 438)
mod1 <- fit.binomialProbitnorm(defaults[,"B"], firms[,"B"])

## Compare the latter two models
c(mod0$maxloglik, mod1$maxloglik)


## More advanced (gives MFE (2005, Table 8.8))

## Prepare data for glmer()
df <- data.frame(defaults = as.vector(t(defaults)), firms = as.vector(t(firms)),
                 year = rep(years, each = 5),
                 rating = factor(rep(rating, 20), levels = rating, ordered = TRUE))
tail(df, n = 10)

## Fit a glmm
(mod <- glmer(cbind(df$defaults, df$firms - df$defaults) ~ -1 + df$rating + (1 | df$year),
              family = binomial(probit)))

## Compute implied asset correlation beta
sigma <- mod@theta
beta <- sigma^2 / (1 + sigma^2) # see last line of Example 11.11 (solved for \beta_i)

## Compute implied default probability p (much smaller than the numbers coming out of equity analysis)
mu <- mod@beta
(PD <- pnorm(mu * sqrt(1 - beta))) # see last line of Example 11.11 (solved for p_i)
