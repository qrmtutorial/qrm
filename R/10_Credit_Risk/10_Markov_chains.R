### By Alexander McNeil

require(Matrix)
# we will need the matrix exponential fuction

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
P <- as.matrix(expm(Lambda))

# number of rating categories
nr <- dim(Lambda)[1]
# sum of migration rates from each non-default rating
migrates <- -diag(Lambda)[-nr]
# names of non-default ratings
nms <- names(migrates)
# probabilities of transition (off diagonal)
transprobs <- Lambda[-nr,]/migrates
# matrix rows sum to zero
apply(transprobs,1,sum)

# generate a random portfolio composition for 5000 firms
set.seed(17)
nfirms <- 5000
fractions <- c(0.1,0.15,0.15,0.2,0.2,0.1,0.1)
initial.ratings <- sample(nms,size=nfirms,replace=TRUE,prob=fractions)
table(initial.ratings)

# generate 20 years of migration data
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
		endtime <- starttime + rexp(n=1,rate=migrates[startrating])		
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



# inspect and summarize rating events
head(RatingEvents,n=8)
summary(RatingEvents)
dim(RatingEvents)






# The continuous-time Markov chain method

Njktable <- table(RatingEvents$startrating,RatingEvents$endrating)
RiskSet <- by(RatingEvents$time,RatingEvents$startrating,sum)
Jlevels <- levels(RatingEvents$startrating)
Klevels <- levels(RatingEvents$endrating)
Njmatrix <- matrix(nrow=length(Jlevels),ncol=length(Klevels),as.numeric(RiskSet),byrow=FALSE)

# Basic form of estimator
Lambda.hat <- Njktable/Njmatrix

# Add default row and correct diagonal
D <- rep(0,dim(Lambda.hat)[2])
Lambda.hat <- rbind(Lambda.hat,D)
diag(Lambda.hat) <- D
rowsums <- apply(Lambda.hat,1,sum)
diag(Lambda.hat) <- -rowsums

# check for valid generator
apply(Lambda.hat,1,sum)

# The matrix exponential
# Annual transition probabilities
P.hat <- as.matrix(expm(Lambda.hat))

# Now compare P.hat and P




