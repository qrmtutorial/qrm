## by Alexander McNeil

## Importance sampling for a one-factor Gaussian threshold model

### Definition of functions #####################################################

## Compute optimal degree of tilting for conditional IS
optimtilt <- function(psi,portfolio,threshold)
{
  tiltfunc <- function(t,psi,portfolio,threshold)
  { 
    exposures <- portfolio$exposures
    probitnorm.mu <- qnorm(portfolio$pd)/sqrt(1-portfolio$beta)
    probitnorm.sigma <- sqrt(portfolio$beta/(1-portfolio$beta))
    pPsi <- pnorm(probitnorm.mu - probitnorm.sigma*psi)
    sum(exposures*exp(t*exposures)*pPsi/(exp(t*exposures)*pPsi+1-pPsi))-threshold
  }
    tmp <- uniroot(tiltfunc,interval=c(-20,20),psi,portfolio,threshold)
    tmp$root
}

## inner tail probability
tailprobinner <- function(psi,threshold,pn.mu,pn.sigma,innerIS,portfolio,n1){
  tfactor <- 0
  if (innerIS)
    tfactor <- optimtilt(psi,portfolio,threshold)
  pPsi <- pnorm(pn.mu - pn.sigma*psi)
  wts <- exp(tfactor*portfolio$exposures)
  numer <- wts*pPsi
  denom <- numer + 1- pPsi
  qPsit <- numer/denom
  ML <- exp(sum(log(denom)))
  loss <- rep(NA,n1)
  for (j in 1:n1){
    loss[j] <- sum(portfolio$exposures*rbinom(length(qPsit),1,qPsit))
  }
  ML*mean(as.numeric(loss >= threshold)*exp(-tfactor*loss))
}

## Find optimal mu
muopt <- function(threshold,pn.mu,pn.sigma,innerIS,portfolio,n1)
{
  optfunc <- function(psi,threshold,pn.mu,pn.sigma,innerIS,portfolio,n1){
    tailprobinner(psi,threshold,pn.mu,pn.sigma,innerIS,portfolio,n1)*exp(-0.5*psi^2)
  }
  tmp <- optimize(f=optfunc,interval=c(-10,10),
                  threshold=threshold,pn.mu=pn.mu,pn.sigma=pn.sigma,
                  innerIS=innerIS,portfolio=portfolio,n1=n1,maximum=TRUE)
  tmp$maximum
}

## Carry out IS
IStailprob <- function(portfolio,threshold,outerIS=TRUE,innerIS=TRUE,n=5000,n.inner=50)
{
  probitnorm.mu <- qnorm(portfolio$pd)/sqrt(1-portfolio$beta)
  probitnorm.sigma <- sqrt(portfolio$beta/(1-portfolio$beta))
  mu <- 0
  n1 <- 1
  if (outerIS) 
    mu <- muopt(threshold,probitnorm.mu,probitnorm.sigma,innerIS,portfolio,100)
  if (innerIS)
    n1 <- n.inner
  theta.inner <- rep(NA,n)
  Psi <- rnorm(n, mean=mu, sd=1)
  for (i in seq(along=Psi)){
    theta.inner[i] <- tailprobinner(Psi[i],threshold,probitnorm.mu,probitnorm.sigma,innerIS,portfolio,n1)
  }
  r.mu <- exp(-mu*Psi+0.5*mu^2)
  theta.hat <- cumsum(theta.inner*r.mu)/(1:n)
  theta.hat
}

### Setup ######################################################################

## Specify portfolio
## m companies, with default probabilities pd and factor loadings beta
m <- 100
portfolio <- list(m=m,exposures=rep(1,m),pd=rep(0.05,m),beta=rep(0.05,m))
## Specify rare event threshold
threshold <- 20

## Test optimal mu
probitnorm.mu <- qnorm(portfolio$pd)/sqrt(1-portfolio$beta)
probitnorm.sigma <- sqrt(portfolio$beta/(1-portfolio$beta))
muopt(threshold,probitnorm.mu,probitnorm.sigma,TRUE,portfolio,100)

### 1 Sampling #########################################################

set.seed(107)
# Full IS
ISresults <- IStailprob(portfolio, threshold, outerIS=TRUE, innerIS=TRUE)
# Outer only
ISresults2 <- IStailprob(portfolio, threshold, outerIS=TRUE, innerIS=FALSE)
# Inner only
ISresults3 <- IStailprob(portfolio, threshold, outerIS=FALSE, innerIS=TRUE)
# Neither
ISresults4 <- IStailprob(portfolio, threshold, outerIS=FALSE, innerIS=FALSE)

## Plot results
## true value - see later for analytical calculation
true.theta <- 0.001121173
index <- 1:length(ISresults)
plot(index,ISresults4,type="l",
     ylim = range(ISresults,ISresults2,ISresults3),col=5,ylab="theta.hat",xlab="n outer simulations")
abline(h=true.theta)
lines(index,ISresults2,col=3)
lines(index,ISresults3,col=4)
lines(index,ISresults,col=2)
legend("topright",legend=c("no IS","outer IS","inner IS","full IS"),col=c(5,3,4,2),lty=1)


### 2 Analytical check on rare event probability ###################################

## pmf of mixed binomial
dprobitnormmix <- function(k,m,mu,sigma){
  integrand <- function(q,k,m,mu,sigma){
    dbinom(k,size=m,prob=q)*(dnorm((qnorm(q) - mu)/sigma))/(dnorm(qnorm(q)) * sigma)
  }
  output <- rep(NA,length(k))
  for (i in (1:length(k))){
    tmp <- integrate(integrand,0,1,k[i],m,mu,sigma)
    output[i] <- tmp$value
  }
  output
}

## Compute and display pmf for given parameter values
allprobs <- dprobitnormmix((0:m),m,probitnorm.mu[1],probitnorm.sigma[1])
sum(allprobs)
probs <- dprobitnormmix((0:(threshold-1)),m,probitnorm.mu[1],probitnorm.sigma[1])
barplot(probs)
## rare event probability
1-sum(probs)




