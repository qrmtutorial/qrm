# Library required


load("INDEXES-2000-2012.RData")
plot(INDEXES0012)
Xdata <- returns(INDEXES0012)

# sample summary statistics
colMeans(Xdata)
cor(Xdata)
cov(Xdata)
var(Xdata)
Sigma <- var(Xdata)
mu <- colMeans(Xdata)
P <- cor(Xdata)

# Cholesky decomposition
chol(Sigma)
A <- t(chol(Sigma))
A %*% t(A)

# Symmetric decomosition
eigen(Sigma)
Gamma <- eigen(Sigma)$vectors
Lambda <- eigen(Sigma)$values
Lambda
Lambda <- diag(sqrt(Lambda))
Lambda
B <- Gamma %*% Lambda %*% t(Gamma)
B %*% B
Sigma


library(QRM)
# simulation mutivariate normal data
simdata <- rmnorm(n=1000,mu=mu,Sigma=Sigma)
colMeans(simdata)
mu
var(simdata)
Sigma
cor(simdata)
P


