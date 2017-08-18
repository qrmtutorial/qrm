## By Marius Hofert

## Introductory R script and playground for learning R


### Comments ###################################################################

## Q: What is R?
## A: - R is a *free* software environment for *statistical* computing and *graphics*
##    - R was created by *R*obert Gentleman and *R*oss Ihaka in 1993.
##    - Since mid-1997, R is developed by the *R Development Core Team* (and
##      contributors)

## Q: Why R?
## A: - *Packages* (both the available ones and the possibility to write
##      your own)
##    - Ability to write *readable* code (focus is on main aspects of a problem)
##    - High(er) level (optimization, run time, debugging, parallel computing etc.)
##

## Q: Where to find R or help on R?
## A: - Main website: https://www.r-project.org/
##    - CRAN: https://cran.r-project.org/
##      + Packages (check 'Published', 'Reference manual', 'Vignettes', 'Package source')
##      + Task Views (check Finance -> rugarch or Multivariate -> copula)
##      + Manuals ("An Introduction to R" (detailed basics) and
##                 "Writing R Extensions" (package development))
##      + FAQ, in particular, FAQ 2.7 on http://cran.r-project.org/doc/FAQ/R-FAQ.html#What-documentation-exists-for-R_003f
##    - Externally (outside R or CRAN):
##      + Google ('r-help') or R-related search engines;
##        see http://tolstoy.newcastle.edu.au/R/ and http://finzi.psych.upenn.edu/
##      + R Help and other mailing lists; see https://stat.ethz.ch/mailman/listinfo/r-help
##      + Stackoverflow (tag 'R'); see http://stackoverflow.com/questions/tagged/r
##        => provide a minimal working example, see https://en. wikipedia.org/wiki/Minimal_Working_Example
##      + R graph gallery; see http://www.r-graph-gallery.com/
##      + R blog; see http://www.r-bloggers.com/
##    - Internally (from within R):
##      + Browser-based help: help.start() -> Search Engine & Keywords
##      + '?' (e.g., ?uniroot) or 'help("[[")' (for specific functions)
##        => Study the examples on the help files!
##      + Study the source code (for more ``hidden'' functions, see pp. 43 in
##        http://cran.r-project.org/doc/Rnews/Rnews_2006-4.pdf)

## Q: How can I work with R?
## A: - Software interfaces:
##      + RStudio (recommended): http://www.rstudio.com/
##      + Emacs + ESS
##    - Workflow:
##      + Write an R script (.R file) containing the source code.
##      + Execute it line-by-line (paragraph-by-paragraph etc.) or the whole script
##        at once (if in batch mode, e.g., on a computer cluster).

## Q: How to install (the latest version of) a package?
## A: - From CRAN (release):
##      install.packages("mypackage")
##    - From R-Forge (snapshot):
##      install.packages("mypackage", repos = "http://R-Forge.R-project.org")

## Q: What to watch out for when programming?
## A: - Theoretical challenges (e.g., curse of dimensionality, e.g., for computing
##      P(a < X <= b) in high dimensions)
##    - Design errors (code correct, but model wrong)
##    - Programming language related issues (R vs C vs Fortran; if you need
##      to implement your own root-finding procedure, errors are more likely)
##    - Syntactic errors (code does not run; typically easy to detect; useful tools:
##      traceback(), debug() or browser())
##    - Semantic errors (code on its own correct, but does not what was intended;
##      test your code, use plots!)
##    - Numerical errors (often undetected unless code is properly tested)
##    - Warnings (are useful! For example, if the optimum in an optimization
##      procedure has not yet been reached)
##    - Measuring run time (system vs wall clock; depends
##      on architecture, programming style, compiler, current workload, etc.; use
##      benchmarks; can also be useful for detecting semantic errors)
##    - Scaling (bigger simulations; if possible, use parallel's mclappy() and
##      parLapply() for multi-core and multi-node computations)

## Note: The code below is a medley of Appendix A of the manual
##       "An Introduction to R" on http://cran.r-project.org/manuals.html


### Simple manipulations; numbers and vectors ##################################

## Simple manipulations
1/2
1/0 # in R, Inf and -Inf exist and R can often deal with them correctly
1/-0
0/0 # ... also NaN = 'not a number' is available; 0/0, 0*Inf, Inf-Inf lead to NaN
x <- 0/0 # store the result in 'x'
class(x) # the class/type of 'x'; => NaN is still of mode 'numeric'
class(Inf) # Inf is of mode 'numeric' (although mathematically not a number); helpful in optimizations

## Vectors (data structure which contains objects of the same mode)
numeric(0) # the empty numeric vector
length(numeric(0)) # its length
x <- c(1, 2, 3, 4, 5) # numeric vector
x # print method
(y <- 1:5) # another way of creating such a vector (and *printing* the output via '()')
(z <- seq_len(5)) # and another one (see below for the 'why')
z[6] <- 6 # append to a vector (better than z <- c(z, 6)); (much) more comfortable than in C/C++
z

## Note: We can check whether the R objects are the same
x == y # component wise numerically equal
identical(x, y) # identical as objects? why not?
class(x) # => x is a *numeric* vector
class(y) # => y is an *integer* vector
all.equal(x, y) # numerical equality; see argument 'tolerance'
identical(x, as.numeric(y)) # => also fine

## Numerically not exactly the same
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
seq_along(c(3, 4, 2)) # 1:3; helpful to 'go along' objects

## Watch out
1:3-1 # ':' has higher priority; note also: the '-1' is recycled to the length of 1:3
1:(3-1)

## Some functions
(x <- c(3, 4, 2))
length(x) # as seen above
rev(x) # change order
sort(x) # sort in increasing order
sort(x, decreasing = TRUE) # sort in decreasing order
ii <- order(x) # create the indices which sort x
x[ii] # => sorted
log(x) # component-wise logarithms
x^2 # component-wise squares
sum(x) # sum all numbers
cumsum(x) # compute the *cumulative* sum
prod(x) # multiply all numbers
seq(1, 7, by = 2) # 1, 3, 5, 7
rep(1:3, each = 3, times = 2) # 1 1 1 2 2 2 3 3 3  1 1 1 2 2 2 3 3 3
tail(x, n = 1) # get the last element of a vector
head(x, n = -1) # get all but the last element

## Logical vectors
logical(0) # the empty logical vector
(ii <- x >= 3) # logical vector indicating whether each element of x is >= 3
x[ii] # use that vector to index x => pick out all values of x >= 3
!ii # negate the logical vector
all(ii) # check whether all indices are TRUE (whether all x >= 3)
any(ii) # check whether any indices are TRUE (whether any x >= 3)
ii |  !ii # vectorized logical OR (is, componentwise, any entry TRUE?)
ii &  !ii # vectorized logical AND (are, componentwise, both entries TRUE?)
ii || !ii # logical OR applied to all values (is entry any TRUE?)
ii && !ii # logical AND applied to all values (are all entries TRUE?)
3 * c(TRUE, FALSE) # TRUE is coerced to 1, FALSE to 0
class(NA) # NA = 'not available' is 'logical' as well (used for missing data)
z <- 1:3; z[5] <- 4 # two statements in one line (';'-separated)
z # => 4th element 'not available' (NA)
(z <- c(z, NaN, Inf)) # append NaN and Inf
class(z) # still numeric (although is.numeric(NA) is FALSE)
is.na(z) # check for NA or NaN
is.nan(z) # check for just NaN
is.infinite(z) # check for +/-Inf
z[(!is.na(z)) &  is.finite(z) &  z >= 2] # indexing: pick out all numbers >= 2
z[(!is.na(z)) && is.finite(z) && z >= 2] # watch out (indexing by 'FALSE' => empty vector)

## Character vectors
character(0) # the empty character vector
x <- "apple"
y <- "orange"
(z <- paste(x, y)) # paste together; use sep = "" or paste0() to paste without space
paste(1:3, c(x, y), sep = " - ") # recycling ("apple" appears again)

## Named vectors
(x <- c("a" = 3, "b" = 2)) # named vector of class 'numeric'
x["b"] # indexing elements by name (useful!)
x[["b"]] # drop the name

## Other types of objects are: arrays (incl. matrices), lists, data frames,
## factors, functions


### Arrays and matrices ########################################################

## Matrices
(A  <- matrix(1:12, ncol = 4)) # watch out, R operates on/fills by *columns*
(A. <- matrix(1:12, ncol = 4, byrow = TRUE)) # fills matrix row-wise
(B <- rbind(1:4, 5:8, 9:12)) # row bind
(C <- cbind(1:3, 4:6, 7:9, 10:12)) # column bind
stopifnot(identical(A, C), identical(A., B)) # check whether the constructions are identical
cbind(1:3, 5) # recycling
(A <- outer(1:4, 1:5, FUN = pmin)) # (4,5)-matrix with (i,j)th element min{i, j}
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
             -8, 5, 3), ncol = 3, byrow = TRUE) # Cholesky factor of the ...
Sigma <- L %*% t(L) # ... real, symmetric, positive definite (covariance) matrix Sigma
standardize <- Vectorize(function(r, c) Sigma[r,c]/(sqrt(Sigma[r,r])*sqrt(Sigma[c,c])))
(P <- outer(1:3, 1:3, standardize)) # construct the corresponding correlation matrix
stopifnot(all.equal(P, cov2cor(Sigma))) # a faster way
P.inv <- solve(P) # compute P^{-1}; solve(A, b) solves Ax = b (if b is omitted, it defaults to I, thus leading to A^{-1})
P %*% P.inv # (numerically close to) I
P.inv %*% P # (numerically close to) I
## Another useful function is Matrix::nearPD(Sigma, corr = TRUE) which finds a
## correlation matrix close to the given matrix in the Frobenius norm.

## Other useful functions
rowSums(A) # row sums
apply(A, 1, sum) # the same
colSums(A) # column sums
apply(A, 2, sum) # the same

## Array (data structure which contains objects of the same mode)
## Special cases: vectors (1d-arrays) and matrices (2d-arrays)
arr <- array(1:24, dim = c(2,3,4),
             dimnames = list(x = c("x1", "x2"),
                             y = c("y1", "y2", "y3"),
                             z = paste("z", 1:4, sep = ""))) # (2,3,4)-array with dimensions (x,y,z)
arr # => also filled in the first dimension first, then the second, then the third
str(arr) # use str() to the *str*ucture of the object arr
arr[1,2,2] # pick out a value
arr. <- aperm(arr, perm = c(3,1,2)) # permute the array to dimensions (z,x,y)
str(arr.)
(mat <- apply(arr, 1:2, FUN = sum)) # for each combination of fixed first and second variables, sum over all other dimensions


### Lists (and data frames) ####################################################

## data frames (lists containing possibly different objects of the same length)
(df <- data.frame(group = rep(LETTERS[1:3], each = 2), value = 1:6))
str(df) # => first column is a factor; second an integer vector

## Creating grids, coercion of data frames to matrices
(grid <- expand.grid(1:3, LETTERS[1:2])[,2:1]) # create a grid containing each variable combination
class(grid) # a data.frame (containing objects of different mode)
as.matrix(grid) # coercion to a character matrix
data.matrix(grid) # coercion to a numeric matrix

## Data frames are just lists
is.list(df)

## Note: Lists are the most general data structures in R in the sense that they
##       can contain pretty much everything, e.g., lists themselves or functions
##       or both... (and of different lengths)
(L <- list(group = LETTERS[1:4], value = 1:2, sublist = list(10, function(x) x+1)))

## Extract elements from a list
## Version 1:
L[[1]] # get first element of the list
L[[3]][[1]] # get first element of the sub-list
## Version 2: # use '$'
L$group
L$sublist[[1]]
## Version 3 (most readable and fail-safe): use the provided names
L[["group"]]
L[["sublist"]][[1]]

## Change a name
names(L)
names(L)[3] <- "sub.list"
str(L)

## Watch out
L[[1]] # the first component
L[1] # the sub-list containing the first component of L
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
y + 2 # ... which is internally again converted to {0,1} in calculations

## Also, loops of the type...
x <- integer(5)
for(i in 1:5) x[i] <- i * i
## ... can typically be avoided by something like
x. <- sapply(1:5, function(i) i * i) # of course we know that this is simply (1:5)^2 which is even faster
stopifnot(identical(x, x.))

## For efficient R programming, the following functions are useful:
## caution, we enter the 'geek zone'...
lapply(1:5, function(i) i * i) # returns a list
sapply(1:5, function(i) i * i) # returns a *s*implified version (here: a vector)
sapply # => calls lapply()
unlist(lapply(1:5, function(i) i * i)) # a bit faster than sapply()
vapply(1:5, function(i) i * i, NA_real_) # even faster but we have to know the return value of the function


### Using implemented distributions ############################################

## Probability distributions (d/p/q/r*)
dexp(1.4, rate = 2) # density f(x) = 2*exp(-2*x)
pexp(1.4, rate = 2) # distribution function F(x) = 1-exp(-2*x)
qexp(0.3, rate = 2) # quantile function F^-(y) = -log(1-y)/2
rexp(4,   rate = 2) # draw random variates from Exp(2)


### Working with additional packages ###########################################

## see ?install.packages()
library(mvtnorm) # load mvtnorm; library = directory location where *packages* reside
library(parallel) # needed for nextRNGStream() below
## Within functions, use require() (throws a warning and continues if the package
## is unavailable) instead of library() (produces an error). See also
## http://yihui.name/en/2014/07/library-vs-require/

packageDescription("mvtnorm") # get a short description of the package
maintainer("mvtnorm") # see citation("mvtnorm") for how to cite a package

## Generate and plot data from a multivariate t distribution
set.seed(271) # for reproducibility (see below)
X <- rmvt(2000, sigma = P, df = 4.5) # generate data from a multivariate t_4.5 distribution
library(lattice) # for the cloud plot
cloud(X[,3]~X[,1]+X[,2], scales = list(col = 1, arrows = FALSE), col = 1,
      xlab = expression(italic(X[1])), ylab = expression(italic(X[2])),
      zlab = expression(italic(X[3])),
      par.settings = list(background = list(col = "#ffffff00"),
                axis.line = list(col = "transparent"), clip = list(panel = "off")))
## => not much visible; in higher dimensions even impossible ...
pairs(X, gap = 0, pch = ".") # ... but we can use a pairs plot


### Random number generation ###################################################

## Generate from N(0,1)
(X <- rnorm(2)) # generate two N(0,1) random variates
str(.Random.seed)
## Note:
## - The first integer in .Random.seed encodes the U(0,1) random number
##   generator kind (lowest two decimals) and the one for generating N(0,1)
##   (highest decimal). The remaining integers denote the actual seed.
## - The default kind is the "Mersenne Twister" (which needs an integer(624)
##   as seed and the current position in this sequence, so 625 numbers).
RNGkind() # => Mersenne Twister, Inversion is used for generating N(0,1)
(Y <- rnorm(2)) # => another two N(0,1) random variates (obviously different)

## How can we make sure to obtain the same results (for *reproducibility*?)
all.equal(X, Y) # obviously not equal (here: with probability 1)

## Set a 'seed' so that computations are reproducible
set.seed(271) # with set.seed() we can set the seed
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
Z. <- rnorm(2) # generate from next stream => will be 'sufficiently apart' from Z
RNGkind("Mersenne-Twister") # switch back to Mersenne-Twister
RNGkind()


### Writing a function #########################################################

##' @title Nonparametric Expected Shortfall Estimator
##' @param x The vector of losses
##' @param alpha The confidence level
##' @param method
##' @param ... Additional arguments passed to quantile(); if 'type = 1',
##'        the quantile of the empirical distribution function is taken.
##' @return Nonparametric ES_alpha estimate (derived under the assumption of continuity)
##' @author Marius Hofert
##' @note - Vectorized in x and alpha
##'       - ">" : Mathematically correct for discrete dfs, but
##'               produces NaN for alpha > (n-1)/n (=> F^-(alpha) = x_{(n)} but
##'               there is no loss strictly beyond x_{(n)})
##'         ">=": mean() will always include the largest loss (so no NaN appears),
##'               but might be computed just based on this one loss.
ES <- function(x, alpha, method = c(">", ">="), ...)
{
    stopifnot(0 < alpha, alpha < 1)
    method <- match.arg(method)
    VaR <- quantile(x, probs = alpha, names = FALSE, ...) # VaR estimate(s); vectorized in x and alpha
    vapply(VaR, function(v)  # v = VaR value for one alpha
        mean(x[if(method == ">") x > v else x >= v]), # mean over all losses >(=) VaR
    NA_real_)
}

## Generate some losses and compute ES_alpha
set.seed(271)
L <- rlnorm(1000, meanlog = -1, sdlog = 2) # L ~ LN(mu, sig^2)
## Note: - meanlog = mean(log(L)) = mu, sdlog = sd(log(L)) = sig
##       - E(L) = exp(mu + (sig^2)/2), var(L) = (exp(sig^2)-1)*exp(2*mu + sig^2)
##         To obtain a sample with E(L) = a and var(L) = b, use:
##         mu = log(a)-log(1+b/a^2)/2 and sig = sqrt(log(1+b/a^2))
ES(L, alpha = 0.99)


### Misc #######################################################################

getwd() # get the current working directory; set it with setwd()

## Not discussed here:
## - How to read/write data from/to a file.
##   This can be done with read.table()/write.table(), for example.
##   For .csv files, there are the convenience wrappers
##   read.csv()/write.csv().
## - How to load/save R objects from/to a file.
##   This can be done using load()/save()
## - How to retrieve an R plot as a file (for printing, for example). For example:
##   doPDF <- TRUE
##   if(doPDF) pdf(file = (file <- "myfile.pdf"), width = 6, height = 6) # open plot device
##   plot(1:10, 10:1) # actual plot command(s)
##   if(doPDF) dev.off() # close plot device; or use crop's dev.off.crop(file)

q() # quit the R session
