# by Alexander McNeil
require(mvtnorm)

# make an equicorrelation matrix P
P <- matrix(0.7,nrow=3,ncol=3)
diag(P) <- 1
P

# generate multivariate normal and Student data
data.normal <- rmvnorm(n=1000, sigma=P)
pairs(data.normal)
data.t <- rmvt(n=1000, df=4, sigma=P)
pairs(data.t)

# use probability transform to transform to uniform
data.gausscopula <- pnorm(data.normal)
data.tcopula <- pt(data.t,df=4)
pairs(data.gausscopula)
pairs(data.tcopula)

# check for uniformity
par(mfrow=(c(1,3)))
apply(data.gausscopula,2,hist)
par(mfrow=c(1,1))


