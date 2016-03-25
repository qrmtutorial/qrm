## By Alexander McNeil

## Fitting a fundamental factor model using industry sector information


### Setup ######################################################################

library(xts) # for time series manipulation
library(qrmdata) # for SP500 (constituents) data
library(qrmtools) # for plot_NA() and plot_matrix()


### 1 Data preparation #########################################################

## Load and extract the data we work with (all available since 1995) 
## Index
data(SP500_const) # index data
SP500.const <- SP500_const['1995-01-01/',] # all since 1995

dim(SP500.const)
plot_NA(SP500.const[,1:100]) # => many stocks have missing values (look at first 100)
SP500.const <- SP500.const[, colSums(is.na(SP500.const)) <= 0.1 * nrow(SP500.const)] # omit columns with more than 10% NA
SP500.const <- na.fill(SP500.const, fill="extend") # fill the remaining NAs
dim(SP500.const)

## Caution: only execute next block if you want to see all the data plotted!!!
par(ask=TRUE)
nplots <- ceiling(dim(SP500.const)[2]/20)
for (i in 1:nplots){
  r1 <- 20*(i-1)+1
  r2 <- min(20*i,dim(SP500.const)[2])
  plot.zoo(SP500.const[,r1:r2], xlab="Time t", main="SP500 Constituents")
}
par(ask=FALSE)

## Remove certain stocks
## AIG, C, ETFC (badly affected by financial crisis)
## T, CTL, FTR, VZ (only 4 telecommunications companies so this sector is thin)
## CMCSK (moves together with CMCSA)
## FOX (moves together with FOXA)
ignore <- list("AIG","T","CTL","C","CMCSK","ETFC","FTR","FOX","VZ")
which(names(SP500.const) %in% ignore)
SP500.const <- SP500.const[,-which(names(SP500.const) %in% ignore)]
dim(SP500.const)

## Build and plot log-returns
X.const <- diff(log(SP500.const))[-1,] # compute -log-returns

## We use monthly data here as basis (and compute monthly log-returns)
X <- apply.monthly(X.const, FUN=colSums) # 252 by 367 matrix
plot.zoo(X[,1:20], type="h", xlab="Time t", main="Monthly risk-factor changes (log-returns) of SP500 constituents")

summary(SP500_const_info)

firms <- as.character(SP500_const_info$Ticker)
## make two firm names consistent
firms[firms=="BRK-B"] <- "BRK.B"
firms[firms=="BF-B"]  <- "BF.B"
sectors <- as.character(SP500_const_info$Sector)

## Construct a factor X.sectors and tabulate levels
X.sectors <- sectors[which(firms %in% colnames(X))]
X.sectors <- as.factor(X.sectors)
table(X.sectors)
levels(X.sectors) <- c("Con-Disc.","Con-Stap.","Energy",
                       "Financials","Health","Industrials","IT",
                       "Materials","Utilities")


### 2 Model fitting ############################################################

## Fit a cross-sectional regression model X_t = B*F_t + eps
## A regression is fitted at each time point
## No intercept is required (-1)

mod <- lm(t(X) ~ X.sectors -1)
B <- model.matrix(mod)

## Factors can be constructed from coefficient matrix
coef.mat <- t(coef(mod))
dimnames(coef.mat)[[2]] <- levels(X.sectors)
F <- xts(coef.mat,time(X))
plot.zoo(F,type="h")

## Calculate error variances for each stock
eps <- t(residuals(mod))
eps.variances <- diag(var(eps))
## Display them from smallest to largest
sort(eps.variances)

## Use generalized least squares to obtain better estimates
## See documentation of lm for more details
mod2 <- lm(t(X) ~ X.sectors -1, weights = 1/eps.variances)
coef.mat <- t(coef(mod2))
dimnames(coef.mat)[[2]] <- levels(X.sectors)
F <- xts(coef.mat,time(X))
plot.zoo(F,type="h")

cor(F)
eps <- t(residuals(mod))
## This matrix should be near-diagonal for perfect factor model
## Of course this won't be the case in practice
cor.eps <- cor(eps)
## Is Cor(eps) (roughly) diagonal? Check upper corner
plot_matrix(cor.eps[1:75,1:75], at=seq(-1, 1, length.out=200)) 

## Are the errors uncorrelated with the factors?
cor.eps.F <- cor(eps, F) 
summary(cor.eps.F) # => there are some large correlations

## Construct the implied covariance and correlation matrix
Ups <- cov(eps) # Upsilon (covariance matrix of epsilon)
Omega <- as.matrix(cov(F)) # Omega (covariance matrix of F; systematic risk)
Sigma <- B %*% Omega %*% t(B) + diag(diag(Ups)) # Cov(X); 
P <- cov2cor(Sigma) # Cor(X)

## Look at discrepancies between the factor model correlation matrix and the
## sample correlation matrix; check upper corner
err <- P-cor(X)
plot_matrix(err[1:75,1:75], at=seq(-1, 1, length.out=200)) 



