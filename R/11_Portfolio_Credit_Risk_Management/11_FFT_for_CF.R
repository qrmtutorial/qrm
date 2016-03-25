# by Alexander McNeil

# Script illustrating calculation of loss distribution using FFT

laplace.transform <- function(t,pd,exposure,lgd=rep(1,length(exposure)))
{
	output <- rep(NA,length(t))
	for (i in 1:length(t))
	output[i] <- exp(sum(log(1-pd*(1- exp(-exposure*lgd*t[i])))))
	output
}


# Specify portfolio characteristics
m <- 20
exposure <- c(5,5,5,5,10,10,10,10,20,20,20,20,30,30,30,30,40,40,40,40)
pd <- c(rep(0.1,10),rep(0.05,10))
# no factor for simplicity

N <- sum(exposure)+1
t <- 2*pi*(0:(N-1))/N
cfunction <- laplace.transform(-t*(1i),pd,exposure) 
plot(cfunction)
# note that the characteristic function is the LT at -i*t

fft.out <- round(Re(fft(cfunction)),digits=20)
sum(fft.out)
probs <- fft.out/N
barplot(probs[(0:20)*5+1],names.arg=paste((0:20)*5))

# check
0.9^10 * 0.95^10
probs[1]
4*0.1*0.9^9 * 0.95^10
probs[6]
6*0.1^2 * 0.9^8 * 0.95^10 + 4*0.1*0.9^9 *0.95^10
probs[11]




# simple binomial example - even easier to check
m <- 20
exposure <- rep(1,m)
pd <- rep(0.3,m)
N <- 21
t <- 2*pi*(0:(N-1))/N
cfunction <- laplace.transform(-t*(1i),pd,exposure)
plot(cfunction)
probs <- round((Re(fft(cfunction)))/N,digits=20)
sum(probs)
plot(0:20,probs[1:21])
dbinom(0:20,20,0.3)
probs



# more complicated example
set.seed(13)
m <- 1000
exposure <- sample(1:10,m,replace=TRUE)
pd <- sample(c(0.05,0.01,0.005,0.001),m,replace=TRUE)
N <- sum(exposure)+1
t <- 2*pi*(0:(N-1))/N
cfunction <- laplace.transform(-t*(1i),pd,exposure) 
plot(cfunction)
fft.out <- round(Re(fft(cfunction)),digits=20)
probs <- fft.out/N
sum(probs*(0:(N-1)))
sum(exposure*pd)
large.losses <- (0:(N-1))[cumsum(probs)>0.995]
large.losses[1]



