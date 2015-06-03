require(Matrix)

# create a generator matrix for an imaginary rating system
# time unit is years
Lambda <- structure(c(-0.0817988760059293, 0.00663655274412117, 0.000600210280876451, 
0.000345909769685178, 0.000293924513908001, 0, 0.000815871122518043, 
0, 0.0749346346627744, -0.0930093347815805, 0.0202980204078218, 
0.00200627666417403, 0.000881773541724003, 0.000762650285432701, 
0, 0, 0.00514818100736618, 0.0815905602071368, -0.0869759261560967, 
0.0417859001779695, 0.00440886770862002, 0.0022879508562981, 
0.00407935561259022, 0, 0.00114404022385915, 0.00361106546371299, 
0.0615488360753306, -0.110345216529572, 0.0534942615312562, 0.00343192628444715, 
0.0057110978576263, 0, 0.000572020111929575, 0.000585578183304809, 
0.00321930968833733, 0.0588046608464803, -0.176844582534647, 
0.0531948574089309, 0.0138698090828067, 0, 0, 0.000487981819420674, 
0.0012004205617529, 0.00615719390039617, 0.108360170794083, -0.190567240072496, 
0.134618735215477, 0, 0, 9.75963638841349e-05, 0.000109129141977537, 
0.000622637585433321, 0.00666228898191469, 0.102671794676377, 
-0.828109189355814, 0, 0, 0, 0, 0.000622637585433321, 0.00274329546314134, 
0.0282180605610099, 0.669014320464795, 0), dim = c(8, 8), dimnames = list(
    c("A1", "A2", "A3", "B1", "B2", "B3", "C", "D"), c("A1", 
    "A2", "A3", "B1", "B2", "B3", "C", "D")))

Lambda
apply(Lambda,1,sum)

# compute annual rating transition probabilities
(P <- expm(Lambda))

stayrates <- -diag(Lambda)
nr <- length(stayrates)
nms <- dimnames(Lambda)[[1]]
names(stayrates) <- nms
(stayrates <- stayrates[-nr])
(nms <- nms[-nr])
(transprobs <- Lambda[-nr,]/stayrates)
apply(transprobs,1,sum)


set.seed(13)
nfirms <- 5000
fractions <- c(0.1,0.15,0.15,0.2,0.2,0.1,0.1)
initial.ratings <- sample(nms,size=nfirms,replace=TRUE,prob=fractions)
table(initial.ratings)

T <- 20
rowmax <- nfirms*5
output <- data.frame(id=numeric(rowmax),starttime=numeric(rowmax),endtime=numeric(rowmax),startrating=rep("",rowmax),endrating=rep("",rowmax),stringsAsFactors=FALSE)
count <- 0
for (i in 1:nfirms)
{
	endtime <- 0
	endrating <- initial.ratings[i]
	while ((endtime < T) & (endrating!="D"))
	{
		count <- count + 1
		starttime <- endtime
		startrating <- endrating
		endtime <- starttime + rexp(n=1,rate=stayrates[startrating])		
		if (endtime <= T){
			pvals <- transprobs[startrating,nms!=startrating]
			endrating <- sample(names(pvals),size=1,prob=pvals)
		}
		if (endtime > T){
			endtime <- T
			}
		output[count,] <- list(i,starttime,endtime,startrating,endrating)
			}
}
output <- output[1:count,]
output$startrating <- as.factor(output$startrating)
output$endrating <- as.factor(output$endrating)
output$time <- output$endtime-output$starttime
RatingEvents <- output




head(RatingEvents,n=50)
tail(RatingEvents,n=50)
summary(RatingEvents)
dim(RatingEvents)






# The continuous-time Markov chain method

(Njktable = table(RatingEvents$startrating,RatingEvents$endrating))
(RiskSet = by(RatingEvents$time,RatingEvents$startrating,sum))
(Jlevels = levels(RatingEvents$startrating))
(Klevels = levels(RatingEvents$endrating))
(Njmatrix = matrix(nrow=length(Jlevels),ncol=length(Klevels),as.numeric(RiskSet),byrow=FALSE))

(Lambda.hat = Njktable/Njmatrix)
D = rep(0,dim(Lambda.hat)[2])
(Lambda.hat = rbind(Lambda.hat,D))
diag(Lambda.hat) = D
(rowsums = apply(Lambda.hat,1,sum))
diag(Lambda.hat) = -rowsums
Lambda.hat

apply(Lambda.hat,1,sum)

# The matrix exponential
# Annual transition probabilities
(P.hat = expm(Lambda.hat))

P




