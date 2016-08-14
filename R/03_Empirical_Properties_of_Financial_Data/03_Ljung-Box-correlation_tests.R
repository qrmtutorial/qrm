# By Alexander McNeil

library(xts)
library(qrmdata)

data("DJ_const")
# We extract a time period and take 29 of 30 DJ stocks
# (Visa has a very short history in the index)
DJdata <- DJ_const['2006-12-29/2015-12-31',-c(27)]
plot.zoo(DJdata)

DJ.X <- diff(log(DJdata))[-1,]
DJ.X.monthly <- apply.monthly(DJ.X,FUN=colSums)

apply(DJ.X,2,Box.test,lag=10,type="Ljung-Box")
apply(abs(DJ.X),2,Box.test,lag=10,type="Ljung-Box")
apply(DJ.X.monthly,2,Box.test,lag=10,type="Ljung-Box")
apply(abs(DJ.X.monthly),2,Box.test,lag=10,type="Ljung-Box")

LBraw.d <- apply(DJ.X,2,function(v,lag,type){Box.test(v,lag,type)$p.value},lag=10,type="Ljung-Box")
LBabs.d <- apply(abs(DJ.X),2,function(v,lag,type){Box.test(v,lag,type)$p.value},lag=10,type="Ljung-Box")
LBraw.m <- apply(DJ.X.monthly,2,function(v,lag,type){Box.test(v,lag,type)$p.value},lag=10,type="Ljung-Box")
LBabs.m <- apply(abs(DJ.X.monthly),2,function(v,lag,type){Box.test(v,lag,type)$p.value},lag=10,type="Ljung-Box")
round(cbind(LBraw.d,LBabs.d,LBraw.m,LBabs.m),2)

# This section gives exactly the Ljung-Box tests in the book
# It user older DJ data in the QRM library

library(QRM)

DJold <- as.xts(DJ)
DJold.X <- diff(log(DJold))[-1,]
DJold.X <- DJold.X['1993-01-01/2000-12-31']
DJold.X.monthly <- apply.monthly(DJold.X,FUN=colSums)

apply(DJold.X,2,Box.test,lag=10,type="Ljung-Box")
apply(abs(DJold.X),2,Box.test,lag=10,type="Ljung-Box")
apply(DJold.X.monthly,2,Box.test,lag=10,type="Ljung-Box")
apply(abs(DJold.X.monthly),2,Box.test,lag=10,type="Ljung-Box")

LBraw.d <- apply(DJold.X,2,function(v,lag,type){Box.test(v,lag,type)$p.value},lag=10,type="Ljung-Box")
LBabs.d <- apply(abs(DJold.X),2,function(v,lag,type){Box.test(v,lag,type)$p.value},lag=10,type="Ljung-Box")
LBraw.m <- apply(DJold.X.monthly,2,function(v,lag,type){Box.test(v,lag,type)$p.value},lag=10,type="Ljung-Box")
LBabs.m <- apply(abs(DJold.X.monthly),2,function(v,lag,type){Box.test(v,lag,type)$p.value},lag=10,type="Ljung-Box")
round(cbind(LBraw.d,LBabs.d,LBraw.m,LBabs.m),2)
