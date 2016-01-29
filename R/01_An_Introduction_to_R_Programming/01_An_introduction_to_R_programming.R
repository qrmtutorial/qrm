## By Marius Hofert

## Introductory R script and playground for learning R

## Note:
## - Appendix A of the manual "An Introduction to R" on
##   http://cran.r-project.org/manuals.html contains another R script which
##   you can work through
## - There are several 'good practices' mentioned on
##   http://www.math.uwaterloo.ca/~mhofert/contents/guidelines.pdf


### Simple manipulations; numbers and vectors ##################################

## Simple manipulations
1+2
1/2
1/0 # in R, Inf and -Inf exist and R can often deal with them correctly
x <- 1/0 # store the result in 'x'
class(x) # the class/type of 'x'

## Vectors (data structure which contains objects of the same mode)
numeric(0) # the empty numeric vector
length(numeric(0)) # its length
x <- c(1, 2, 3, 4, 5) # numeric vector
x # print method
(y <- 1:5) # another way of creating such a vector (and *printing* the output via '()')
(z <- seq_len(5)) # and another one (see below for the 'why')
z[6] <- 6 # append to a vector (better than z <- c(z, 6))
z

## Note: We can check whether the R objects are the same
x == y # component wise
identical(x, y) # identical as objects? why not?
class(x) # => x is a *numeric* vector
class(y) # => y is an *integer* vector
all.equal(x, y) # numerical equality; see argument 'tolerance'
identical(x, as.numeric(y)) # => also fine

## Another example of this type
x <- var(1:4)
y <- sd(1:4)^2
all.equal(x, y) # numerical equality
x == y # ... but not exactly
x-y # numerically not 0
## See also https://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-doesn_0027t-R-think-these-numbers-are-equal_003f

## Watch out
n <- 0
1:n # not the empty sequence but c(1, 0); caution in 'for loops': for(i in 1:n) ...!
seq_len(n) # better: => empty sequence
seq_along(c(3,4,2)) # 1:3; helpful to 'go along' objects

## Watch out
1:3-1 # ':' has higher priority
1:(3-1)

## Vector arithmetic
x <- c(1, 2, 3, 4, 5)
(z <- 2*x - 1) # '*' is component-wise, '-1' is recycled to the length of x

## Some functions
(x <- c(3, 4, 2))
rev(x) # change order
sort(x) # sort in increasing order
sort(x, decreasing=TRUE) # sort in decreasing order
o <- order(x) # create indices that sort x
x[o] # => sorted
length(x) # length of x
log(x) # component-wise logarithms
x^2 # component-wise squares
sum(x) # sum all numbers
cumsum(x) # compute the *cumulative* sum
prod(x) # multiply all numbers
seq(1, 7, by=2) # 1, 3, 5, 7
rep(1:3, each=3, times=2) # 1 1 1 2 2 2 3 3 3  1 1 1 2 2 2 3 3 3

## Logical vectors
logical(0) # the empty logical vector
(ind <- x >= 3) # logical vector indicating whether each element of x is >= 3
x[ind] # use that vector to index x => pick out all values of x >= 3
!ind # negate the logical vector
all(ind) # check whether all indices are TRUE (whether all x >= 3)
any(ind) # check whether any indices are TRUE (whether any x >= 3)
ind | !ind # vectorized logical OR
ind & !ind # vectorized logical AND
ind || !ind # logical OR applied to all values
ind && !ind # logical AND applied to all values
y <- c(TRUE, FALSE)
3*y # TRUE is coerced to 1, FALSE to 0
class(NA) # NA = 'not available' is 'logical' as well

## Missing values (NA), not a number (NaN)
z <- 1:3; z[5] <- 4 # two statements in one line (';'-separated)
z # => 4th element 'not available' (NA)
(z <- c(z, 0/0)) # e.g., 0/0, 0*Inf, Inf-Inf lead to 'not a number' (NaN)
class(NaN) # not a number but still of mode 'numeric'
class(z) # still numeric
is.na(z) # check for NA or NaN
is.nan(z) # check for just NaN
z[(!is.na(z)) & z >= 2] # indexing: pick out all numbers >= 2
z[(!is.na(z)) && z >= 2] # watch out (indexing by 'FALSE' => empty vector)

## Character vectors
character(0) # the empty character vector
x <- "apple"
y <- "orange"
(z <- paste(x, y)) # paste together; use sep="" or paste0() to paste without space
paste(1:3, c(x, y), sep=" - ") # recycling ("apple" appears again)

## Named vectors
(x <- c("a"=3, "b"=2)) # named vector of class 'numeric'
x["b"] # indexing elements by name (useful!)
x[["b"]] # drop the name

## Other types of objects are: arrays (incl. matrices), lists, data frames,
## factors, functions


### Arrays and matrices ########################################################

## Matrices
(A  <- matrix(1:12, ncol=4)) # watch out, R operates on/fills by *columns*
(A. <- matrix(1:12, ncol=4, byrow=TRUE)) # fills matrix row-wise
(B <- rbind(1:4, 5:8, 9:12)) # row bind
(C <- cbind(1:3, 4:6, 7:9, 10:12)) # column bind
stopifnot(identical(A, C), identical(A., B)) # check whether the constructions are identical
cbind(1:3, 5) # recycling
(A <- outer(1:4, 1:5, FUN=pmin)) # build a (4, 5)-matrix with (i,j)th element being min{i, j}
## => Lower triangular matrix contains column number, upper triangular matrix contains row number

## Some functions
nrow(A) # number of rows
ncol(A) # number of columns
dim(A) # dimension
diag(A) # 1 2 3 4; diagonal of A
diag(3) # identity (3, 3)-matrix
(D <- diag(1:3)) # diagonal matrix with elements 1, 2, 3
D %*% B # matrix multiplication
B * B # Hadamard product, i.e., element-wise product

## Build a correlation matrix and invert it
L <- matrix(c(2, 0, 0,
              6, 1, 0,
             -8, 5, 3), ncol=3, byrow=TRUE) # Cholesky factor of the ...
Sigma <- L %*% t(L) # ... real, symmetric, positive definite (covariance) matrix x
standardize <- Vectorize(function(r, c) Sigma[r,c]/(sqrt(Sigma[r,r])*sqrt(Sigma[c,c])))
P <- outer(1:3, 1:3, standardize) # construct the corresponding correlation matrix
## Alternatively, this could have been done with Matrix::nearPD(Sigma, corr=TRUE)
## which works slightly differently though (by finding a correlation matrix
## close to the given matrix in the Frobenius norm) and thus gives a different answer.
P.inv <- solve(P) # compute P^{-1}; solve(A, b) solves Ax=b (if b is omitted, it defaults to I, thus leading to A^{-1})
P %*% P.inv # (numerically close to) I
P.inv %*% P # (numerically close to) I

## Build a grid and work on it
(grid <- expand.grid(1:3, LETTERS[1:2])[,2:1]) # create a grid containing each variable combination
class(grid) # a data.frame (containing objects of different mode)
as.matrix(grid) # coercion to matrix
data.matrix(grid) # numeric matrix (without quotes)
rowSums(A) # row sums
apply(A, 1, sum) # the same
colSums(A) # column sums

## Array (data structure which contains objects of the same mode)
## Special cases: vectors (1d-arrays) and matrices (2d-arrays)
arr <- array(1:24, dim = c(2,3,4),
             dimnames = list(x = c("x1", "x2"),
                             y = c("y1", "y2", "y3"),
                             z = paste("z", 1:4, sep=""))) # (2,3,4)-array with dimensions (x,y,z)
arr # => also filled in the first dimension first, then the second, then the third
str(arr) # use str() to the *str*ucture of the object arr
arr[1,2,2] # pick out a value
arr. <- aperm(arr, perm=c(3,1,2)) # permute the array to dimensions (z,x,y)
str(arr.)
(mat <- apply(arr, 1:2, FUN=sum)) # for each combination of fixed first and second variables, sum over all others (the third dimension)


### Lists and data frames ######################################################

## data.frame (data structure which contains objects of the same length but
## possibly different type)
(df <- data.frame(group=rep(LETTERS[1:3], each=2), value=1:6))
str(df) # => first column is a factor; second an integer vector

## Note: Lists are the most general data structures in R in the sense that they
##       can contain pretty much everything, e.g., lists themselves or functions
##       or both... (and of different lengths)
(L <- list(group=LETTERS[1:4], value=1:2, sublist=list(10, function(x) x+1)))
length(L) # length of the list

## Extract elements from a list
## Version 1:
L[[1]] # get first element of the list
L[[3]][[1]] # get first element of the sub-list
## Version 2: # use '$'
L$group
L$sublist[[1]]
## Version 3: use the provided names
L[["group"]]
L[["sublist"]][[1]]

## Change a name
names(L)[3] <- "sub.list"
str(L)

## Watch out
L[[1]] # the first component
L[1] # the sub-list consisting of the first component of L
class(L[[1]]) # character
class(L[1]) # list


### Control statements (just very quickly) #####################################

## R has if() else, ifelse() (a vectorized version of 'if'), for loops (avoid or
## only use if they don't take much run time), repeat and while (with 'break' to
## exit and 'next' to advance to the next loop iteration)

## ... without going into details, note that even 'if()' is a function, so
## instead of:
x <- 4
if(x < 5) y <- 1 else y <- 0 # y is the indicator whether x < 5
## ... write (the much more readable)
y <- if(x < 5) 1 else 0
## ... or even better
(y <- x < 5) # ... as a logical
as.numeric(y) # converted to 1, but even without that...
y+2 # ... y is treated as {0,1} in calculations

## Also, loops of the type...
x <- numeric(5)
for(i in 1:5) x[i] <- i*i
## ... can typically be avoided by something like
x <- sapply(1:5, function(i) i*i) # of course we know that this is simply (1:5)^2 which is even faster

## For efficient R programming, the following functions are useful:
## caution, we enter the 'geek zone'...
lapply(1:5, function(i) i*i) # returns a list
sapply(1:5, function(i) i*i) # returns a *s*implified version (here: a vector)
unlist(lapply(1:5, function(i) i*i)) # a bit faster than sapply()
vapply(1:5, function(i) i*i, NA_real_) # even faster but we have to know the return value of the function


### Using implemented distributions ############################################

## Probability distributions (d/p/q/r*)
dexp(1.4, rate=2) # density f(x) = 2*exp(-2*x)
pexp(1.4, rate=2) # distribution function F(x) = 1-exp(-2*x)
qexp(0.3, rate=2) # quantile function F^-(y) = -log(1-y)/2
rexp(4,   rate=2) # draw random variates from Exp(2)


### Working with additional packages ###########################################

library(mvtnorm) # load mvtnorm; library = directory location where *packages* reside
library(parallel) # needed for nextRNGStream() below
## Within functions, use require() (throws a warning and continues if the package
## is unavailable) instead of library() (produces an error). See also
## http://yihui.name/en/2014/07/library-vs-require/

packageDescription("mvtnorm") # get a short description of the package
citation("mvtnorm") # how to cite a package
## Generate and plot data from the multivariate t_{4.5} distribution
pairs(rmvt(2000, sigma=diag(3), df=4.5), gap=0, pch=".")


### Random number generation ###################################################

## Remove .Random.seed (just in case it was set)
rm(.Random.seed) # remove .Random.seed (not existing yet if you started a new R session)
.Random.seed # => not there anymore; that's also the case in a new R session (try it!)

## Generate from N(0,1)
(X <- rnorm(2)) # generate two N(0,1) random variates
str(.Random.seed)
## => R generated .Random.seed which encodes the random number generator kind
##    (lowest two decimals) and the random number generator kind for generating
##    N(0,1) (highest decimal) in the first integer and then contains the seed
##    (all other components). The default kind is the "Mersenne Twister"
##    (which needs an integer(624) as seed and the current position in this
##    set => 625 numbers).
RNGkind() # => Mersenne Twister, Inversion is used for generating N(0,1)
(Y <- rnorm(2)) # => another two N(0,1) random variates (different from above)

## How can we make sure to obtain the same results (for *reproducibility*?)
all.equal(X, Y) # obviously not equal (here: with probability 1)

## Set a 'seed' so that computations are reproducible
rm(.Random.seed) # remove .Random.seed again (for demonstration purposes)
set.seed(271) # with set.seed() we can set the seed
str(.Random.seed) # => again .Random.seed now exists
X <- rnorm(2) # draw two N(0,1) random variates
set.seed(271) # set the same seed again
Y <- rnorm(2) # draw another two N(0,1) random variates
all.equal(X, Y) # => TRUE
## If you just start R without calling set.seed(), a seed is constructed from
## system time and the R process number.

## A pseudo-random number generator which allows for easily advancing the
## seed is L'Ecuyer's combined multiple-recursive generator (CMRG); see
## MRG32k3a.c and MRG32k3a.h on http://simul.iro.umontreal.ca/rng (R's version
## only makes minor modifications to this). Let's see how we can call it from R.
RNGkind() # => Mersenne Twister, Inversion is used for generating N(0,1)
RNGkind("L'Ecuyer-CMRG")
RNGkind() # => L'Ecuyer's CMRG, Inversion is used for generating N(0,1)
.Random.seed # => now of length 7: first number as above + the seed
Z <- rnorm(2) # use L'Ecuyer's CMRG for generating random numbers
.Random.seed <- nextRNGStream(.Random.seed) # advance seed by 2^127; requires 'parallel'
Z. <- rnorm(2) # generate from next stream
RNGkind("Mersenne-Twister") # switch back to Mersenne-Twister
str(.Random.seed)


### Watch out for numerical issues #############################################

## How to evaluate choose(500, 200)?
choose(500, 200)
factorial(500)/(factorial(200)*factorial(300)) # too large values
n <- 170 # last n possible
x <- 1:n
y <- sapply(1:n, factorial)
plot(x, y, type="l", log="y", xlab="n", ylab="n!") # ... always look at plots
y[length(y)] # close to largest value...
str(.Machine) # => double.xmax

## A numerical trick
log(factorial(200)) # obviously, same problem here ('mathematical composition' not doable)
lc <- lfactorial(500) - (lfactorial(200) + lfactorial(300)) # work with 'proper' logs
stopifnot(all.equal(lc, lchoose(500, 200))) # there is also lchoose()
c <- exp(lc) # => doable
choose(500, 200) # => choose() uses the same trick
stopifnot(all.equal(c, choose(500, 200)))


### Misc #######################################################################

getwd() # get the current working directory; set it with setwd()

## Not discussed here:
## 1) How to read/write data from/to a file.
##    This can be done with read.table()/write.table(), for example.
##    For .csv files, there are the convenience wrappers
##    read.csv()/write.csv().
## 2) How to load/save R objects from/to a file.
##    This can be done using load()/save()

q() # quit the R session