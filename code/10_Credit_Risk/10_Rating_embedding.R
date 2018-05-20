#### by Alexander McNeil
library(expm)
Aaa.rates = c(0,0.082,0.0063,0,0.0003,0,0,0,0)
Aa.rates = c(0.0091,0,0.0843,0.0049,0.0006,0.0002,0.0001,0,0.0002)
A.rates = c(0.0006,0.0248,0,0.0547,0.0057,0.0011,0.0003,0,0.0006)
Baa.rates = c(0.00039,0.0017,0.0411,0,0.0405,0.0755,0.0163,0.0002,0.0017)
Ba.rates = c(0.0001,0.0005,0.0035,0.0552,0,0.0722,0.0058,0.0007,0.0106)
B.rates = c(0.0001,0.0003,0.0011,0.0032,0.0458,0,0.0581,0.0059,0.0385)
Caa.rates = c(0.0001,0.0002,0.0002,0.0012,0.0038,0.0870,0,0.0372,0.1334)
C.rates = c(0,0,0,0,0.004,0.0203,0.0938,0,0.3793)
D.rates = rep(0,length(C.rates))

Lambda = rbind(Aaa.rates,Aa.rates,A.rates,Baa.rates,Ba.rates,B.rates,Caa.rates,C.rates,D.rates)
rnames = c("Aaa","Aa","A","Baa","Ba","B","Caa","C","D")
dimnames(Lambda) =list(rnames,rnames)
Lambda.d = -apply(Lambda,1,sum)
diag(Lambda) = Lambda.d
apply(Lambda,1,sum)
Lambda

P = expm(Lambda)
apply(P,1,sum)

P2 <- P[9:1,9:1]
P2
P2 <- P2[-1,]
P2


# pick C
Cprobs <- P2[1,]
Cprobs
Ccumprobs <- cumsum(Cprobs)
Ccumprobs
thresholds <- qnorm(Ccumprobs)
thresholds

plot(seq(from=-5,to=5,length=100),dnorm(seq(from=-5,to=5,length=100)),type="l",xlab="X",ylab="density")
abline(v=thresholds,col=1:length(thresholds))

# do all
cumprobs <- t(apply(P2,1,function(v){cumsum(v)}))
cumprobs
thresholds<-qnorm(cumprobs)
thresholds

opa <- par(mfrow=c(2,4))
for (j in 1:nrow(thresholds))
{
  plot(seq(from=-5,to=5,length=100),dnorm(seq(from=-5,to=5,length=100)),type="l",
       xlab="X",ylab="density",main=rownames(thresholds)[j])
  abline(v=thresholds[j,],col=1:length(thresholds[j,]))
}
par(opa)
