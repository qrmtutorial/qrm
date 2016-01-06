## by Alexander McNeil
require(rugarch)
# We use the rugarch package which uses the S4 version of the S language

# Specify a GARCH(1,1) model
GARCHspec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),mean.model=list(armaOrder=c(0,0),include.mean=FALSE),distribution.model="norm",fixed.pars=list(omega=0.02,alpha1=0.15,beta1=0.8))
GARCHspec

# Generate two realizations of length 2000
path <- ugarchpath(GARCHspec,n.sim=2000,n.start=50,m.sim=2)
path
class(path)
# path is an S4 object of class "uGARCHpath". 
# We can use getClass to see more information about such objects
# We can also use the getSlots function to see the composition.
getClass("uGARCHpath")
getSlots("uGARCHpath")

# There is a special plotting function for "uGARCHpath" objects
plot(path,which=1)
plot(path,which=2)
plot(path,which=3)
plot(path,which=4)
# How to see the documentation of the plot function
showMethods("plot")
getMethod("plot",c(x="uGARCHpath", y="missing"))

# There is also an extraction function for the volatility
vol <- sigma(path)
head(vol)
plot(vol[,1],type="h")
plot(vol[,2],type="h")

series <- path@path
# series is a simple list
class(series)
names(series)
# the actual simulated data are in the matrix/vector called "seriesSim"
X <- series$seriesSim
head(X)

# Does simulated series conform to stylized facts?
X1 <- X[,1]
acf(X1)
acf(abs(X1))
qqnorm(X1)
qqline(X1,col=2)
shapiro.test(X1)

# Exercise. Try simulating paths with different innovation distributions (e.g. Student t)
