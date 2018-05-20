### by Alexander McNeil

library(sde)

# We will price an interest-rate swap in the CIR model
# and see how the value and hence the counterparty credit risk
# evolves for one of the counterparties

# Assume the short rate evolves according to the CIR model
# Using Brigo-Mercurio notation

# Dynamics under Q, the risk-neutral measure
# Long term level
theta <- 0.1
# Speed of mean reversion
k <- 0.02
# Volatility
sigma <- 0.04
# Positivity check
(2*k*theta - sigma^2 >0)

# Dynamics under P, the real-world measure
lambda <- 0.5
k.bar <- k + lambda*sigma
theta.bar <- k*theta/(k+lambda*sigma)
# Positivity check
(2*k.bar*theta.bar - sigma^2 >0)

# Functions for affine representation of term structure
# parameters k, theta, sigma fixed as above throughout
Afunc <- function(t,T){
  h=sqrt(k^2 + 2*sigma^2)
  top <- 2*h*exp((k+h)*(T-t)/2)
  bottom <- 2*h+(k+h)*(exp((T-t)*h)-1)
  (top/bottom)^{2*k*theta/(sigma^2)}
}  	
Bfunc <- function(t,T){
  h=sqrt(k^2 + 2*sigma^2)
  top <- 2*(exp((T-t)*h)-1)
  bottom <- 2*h+(k+h)*(exp((T-t)*h)-1)
  top/bottom
}

# functions for bond prices and yields
price <- function(t,T,rt){
  Afunc(t,T)*exp(-Bfunc(t,T)*rt) 
}
yield <- function(t,T,rt,k){
  -log(price(t,T,rt))/(T-t) 
}

# Current short rate
r0 = 0.005


# Simulate short rate under P
# Have to hardwire parameters into drift and vol functions
(dfunc <- parse(text=paste(k.bar,"*(",theta.bar,"-x)",sep="")))
(sfunc <- parse(text=paste(sigma,"*sqrt(x)",sep="")))

# Simulate a number of short rate paths starting at current value of r0
n.paths <- 20
# Time horizon in years
Thorizon <- 0.5
# and in days
n.days <- floor(Thorizon*365)
# Simulate paths
set.seed(13)
paths <- sde.sim(t0 = 0, T=Thorizon, N=n.days, X0=r0, drift=dfunc, sigma=sfunc, M=n.paths)
head(paths)

# Show short rate paths to time horizon
times <- (0:n.days)/365
op <- par(mar=c(3,3,2,1),mgp=c(2,1,0))
plot(times,paths[,1],type="l",xlab="t",ylab=expression(r[t]),ylim=range(paths))
for (i in 1:n.paths){
  lines(times,paths[,i],col=i)
}



# We consider a prototypical (forward-start) Interest-Rate Swap
# This is product on page 13 of Brigo-Mercurio
# start date Talpha 
# end date Tbeta 
# year fraction describing payment intervals is tau 
# nominal is N 
tau <- 0.25
N <- 1000
Talpha <- 1
Tbeta <- 4
# payment times
Tk <- seq(from=Talpha+tau, to =Tbeta, by = tau)

# A receives fixed and pays floating
# A has a receiver IRS (RFS)

# B receives floating and pays fixed
# B has a payer IRS (PFS)

PFS <- function(K,Talpha,Tbeta,Tk,N,tau,t,rt){
  N*(price(t,Talpha,rt)-price(t,Tbeta,rt))  - N*tau*K*sum(price(t,Tk,rt))  
}

# determine fixed interest rate that makes contract fair at outset
tmp <- uniroot(PFS,c(0.0001,0.1),Talpha=Talpha,Tbeta=Tbeta,Tk=Tk,N=N,tau=tau,t=0,rt=r0)
K <- tmp$root

# Initial value of PFS (zero)
PFS(K,Talpha,Tbeta,Tk,N,tau,0,r0)


# Swap value for B, i.e. value of PFS
swap.value <- matrix(NA,ncol=n.paths,nrow=length(times))
for (j in 1:n.paths){
  rt <- paths[,j]
  for (i in 1:length(times)){
    swap.value[i,j] <- PFS(K,Talpha,Tbeta,Tk,N,tau,times[i],rt[i])
  }
}
head(swap.value)

plot(times,swap.value[,1],type="l",xlab="t",ylab=expression(V[t]),ylim=range(swap.value))
for (i in 1:n.paths){
  lines(times,swap.value[,i],col=i)
}
par(op)

# initial view of the term structure of interest rates for key maturities
prices0 <- price(0,Tk,r0)
yields0 <- yield(0,Tk,r0)

op <- par(mfrow=c(4,1),mar=c(3,3,2,1),mgp=c(2,1,0))
# A scenario where interest rates rise
j <- 16
plot(times,paths[,j],type="l",xlab="t",ylab=expression(r[t]))
plot(times,swap.value[,j],type="l",xlab="t",ylab=expression(V[t]))
# Show change in term structure of prices over time horizon
prices1 <- price(Thorizon,Tk,paths[length(times),j])
plot(Tk,prices0, xlab="T", ylab=expression(p(t,T)),type="b",ylim=range(prices0,prices1))
lines(Tk,prices1,type="b",col=2)
# Show change in term structure of yields over time horizon
yields1 <- yield(Thorizon,Tk,paths[length(times),j])
plot(Tk,yields0, xlab="T", ylab=expression(y(t,T)),type="b",ylim=range(yields0,yields1))
lines(Tk,yields1,type="b",col=2)

# A scenario where interest rates fall
j <- 7
plot(times,paths[,j],type="l",xlab="t",ylab=expression(r[t]))
plot(times,swap.value[,j],type="l",xlab="t",ylab=expression(V[t]))
# Show change in term structure of prices over time horizon
prices1 <- price(Thorizon,Tk,paths[length(times),j])
plot(Tk,prices0, xlab="T", ylab=expression(p(t,T)),type="b",ylim=range(prices0,prices1))
lines(Tk,prices1,type="b",col=2)
# Show change in term structure of yields over time horizon
yields1 <- yield(Thorizon,Tk,paths[length(times),j])
plot(Tk,yields0, xlab="T", ylab=expression(y(t,T)),type="b",ylim=range(yields0,yields1))
lines(Tk,yields1,type="b",col=2)

par(op)



