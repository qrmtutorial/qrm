# By Alexander McNeil
# This script gives exactly the Ljung-Box tests in the book

require(QRM)

DJ.X <- returns(DJ)
DJ.X <- window(DJ.X,start=timeDate("1993-01-01"),end=timeDate("2000-12-31"))

apply(DJ.X,2,Box.test,lag=10,type="Ljung-Box")
apply(abs(DJ.X),2,Box.test,lag=10,type="Ljung-Box")

by <- unique(timeLastDayInMonth(time(DJ.X)))
DJ.X.monthly <- aggregate(DJ.X,by,sum)

apply(DJ.X.monthly,2,Box.test,lag=10,type="Ljung-Box")
apply(abs(DJ.X.monthly),2,Box.test,lag=10,type="Ljung-Box")