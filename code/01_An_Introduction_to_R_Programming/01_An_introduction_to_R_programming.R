## By Marius Hofert

## Introductory R script and playground for learning R; see also Appendix A of
## the manual "An Introduction to R" on http://cran.r-project.org/manuals.html


### Comments ###################################################################

## Q: What are key concepts I should be aware of when programming?
## A: 1) An *algorithm* is a well-defined (unambiguous) finite set of
##       instructions (steps) for solving a problem.
##       Algorithms are typically formulated in *pseudo-code* (a high-level
##       description of an algorithm, ideally next to you when implementing it),
##       e.g. the Sieve of Eratosthenes for finding all primes <= n.
##       An *implementation* of an algorithm allows one to see how a problem is
##       solved (as opposed to a vague formulation that pseudo-code typically is).
##    2) A *minimal working example (MWE)* is
##       ... a *working example* in the sense that it allows someone else to
##           reproduce a problem or result (sufficiency)
##       ... *minimal* in the sense that it is as small as possible, without
##           non-relevant code, data or dependencies (necessity).
##       MWEs are indispensable for finding bugs or posting problems to forums;
##       see https://en.wikipedia.org/wiki/Minimal_Working_Example
##    3) Programming is learned by doing it. There are two approaches:
##       - bottom-up: Pick a problem you like and try to solve it by writing a
##                    program.
##       - top-down: Start with code you found, alter it and understand the
##                   consequences.
##    4) Software development (e.g. a package) is significantly more complicated
##       than writing standalone R scripts (maintenance, no hardcoded values,
##       initial values, unknown users (ab)use functions etc.)

## Q: What is R?
## A: - R is a *free* software environment for *statistical computing* and *graphics*
##    - R was created by *R*obert Gentleman and *R*oss Ihaka in 1993.
##    - Since mid-1997, R is developed by the *R Development Core Team*
##      (and contributors)
##    - The sources or R (.tar.gz) are on https://cran.r-project.org/ and
##      consist of base and recommended packages (involving C and Fortran code,
##      too); see installed.packages()[,"Priority"]
##    - Download and unpack it, browse /src/library/base/R (e.g., constants.R or
##      mean.R) or /src/library/parallel/R (e.g. detectCores.R)

## Q: Why R?
## A: - *Statistical* software (e.g., time series modeling, goodness-of-fit tests,...)
##    - *Contributed packages*; see
##      https://cran.r-project.org/web/packages/available_packages_by_name.html
##    - The possibility to write your own *package*; e.g.
##      https://cran.r-project.org/web/packages/qrmtools/index.html
##    - Ability to write *readable* code (focus is on main aspects of a problem)
##    - *High-level* programming language (plotting, optimization, run time
##      measurement, debugging, parallel computing etc.). Helps when exploring a
##      problem to gain understanding and insight.

## Q: How to install R (for daily use)?
## A: - Install *R* from https://cran.r-project.org/
##    - Install an *integrated development environment (IDE)* (e.g., editor, build tools):
##      + RStudio: http://www.rstudio.com/
##      + Emacs + ESS (Emacs Speaks Statistics)

## Q: How to install R (contributed) packages?
## A: - *Release* versions:
##      + Comprehensive R Archive Network (CRAN): install.packages("mypkg")
##      + Bioconductor: devtools::install_bioc("mypackage")
##                      or BiocManager::install("mypackage")
##    - *Development* versions:
##      + R-Forge: install.packages("mypkg", repos = "https://R-Forge.R-project.org")
##      + GitHub: devtools::install_github("maintainer/mypkg")
##    - From *source* (.tar.gz):
##      install.packages("mypkg.tar.gz", repos = NULL, lib = "/mydirectory")
##      (if lib is not provided, R takes the first element in .libPaths() which
##      can be set via R_LIBS_USER=/mydirectory in ∼/.Renviron).

## Q: How to work with R?
## A: - Create a *script* (a file, e.g. myscript.R) containing the R source code.
##    - Run the script interactively by executing it *line-by-line* (Ctrl + RET)
##    - For bigger simulations, run the script in *batch mode*:
##      + 'R CMD BATCH myscript.R' (output contains whole R session;
##        better in projects, e.g., when run on a cluster)
##      + 'Rscript myscript.R > myscript.Rout' (output only contains actual
##        printed output; better for shell-type R scripts with first
##        line #!/usr/bin/env Rscript)

## Q: How to find help on R?
## A: - From within an R session, manual pages of functions: '?' (e.g., ?uniroot) or
##      'help("[[")' (for specific functions). Study the examples on the help files.
##    - Online:
##      + Google ('r-help')
##      + Stackoverflow (tag 'R')
##      + mailing lists, e.g. https://stat.ethz.ch/mailman/listinfo/r-help
##      => Provide a MWE
##    - CRAN: https://cran.r-project.org/
##      + Task Views (select packages on certain topics)
##      + Manuals ("An Introduction to R" (detailed basics) and
##                 "Writing R Extensions" (package development))
##      + FAQ
##      + Packages (check 'Published' (date), 'Reference manual' (all help files),
##        'Vignettes' (PDF or HTML explanations of features), 'Package source' (.tar.gz))
##    - Study the source code (see pp. 43 in http://cran.r-project.org/doc/Rnews/Rnews_2006-4.pdf)
##      + Within an R session, execute the function and go on from there:
##        - Exported functions,
##          e.g. ncol, mvtnorm::rmvt, optimize
##        - mypkg:::myfun (for non-exported functions, i.e. functions not intended
##          to be available to the user directly),
##          e.g. qrmtools:::seq2 (only seq2 would fail even after library(qrmtools))
##        - methods(myfun) (for finding S3 and S4 methods),
##          e.g. mean -> methods(mean) -> mean.default
##        - getAnywhere(myfun) and showMethods(myfun, includeDefs = TRUE) for
##          showing S3 and S4 methods, respectively,
##          e.g. mean.zoo -> getAnywhere("mean.zoo") -> zoo::mean.zoo -> zoo:::mean.zoo
##          or library(copula) -> showMethods("rCopula") ->
##             showMethods("rCopula", classes = "claytonCopula", includeDefs = TRUE)
##      + Download the source mypkg.tar.gz and look/search inside. This must
##        be used for .Primitive(), .Internal(), .Call() etc.

## Q: What to watch out for when implementing a model?
## A: - Using the *wrong model* or *assumptions* (=> model risk),
##      e.g. predicting nuclear accidents based on too few data.
##    - Using the *wrong software* (too low-level, no standard tools like optimizers
##      available) or *wrongly using software* (e.g. recycling, sample size 1, calling
##      functions the wrong way because of misinterpreting the meaning of arguments).
##    - Theoretical hurdles, e.g. when computing P(a < X <= b) analytically in
##      high dimensions
##    - Software related errors:
##      + *Syntax errors*: source code does not compile/interpreter error because
##        of a violation of the syntax of the programming language. Functions like
##        traceback(), debugonce() or browser() are helpful to sort out syntax
##        errors in R.
##      + *Run-time errors*: correct code/syntax but errors appear when the program
##        is run (e.g. division by a numerical 0). Often numerical errors (i.e.
##        errors due to the floating-point representation of numbers).
##      + *Semantic errors*: correct code/syntax and runs correctly, but the
##        program does not compute what is intended (e.g. through a switch in a
##        logical statement). Test your code, use plots! Measure run time can
##        also be useful.
##      + *User errors*: errors caused by users of the software (e.g. wrong inputs).
##    - *Numerical* issues (sometimes kick in very slowly, only in some
##      iterations, appearing randomly etc.)
##    - *Reproducibility* (not using a seed leads to lack of reproducibility)
##    - Programs can be run *sequentially* (i.e., one computation at a time,
##      typically on one CPU) or in *parallel* (i.e. the work is shared among
##      workers (multiple cores in a CPU or multiple CPUs in a cluster),
##      coordinated by one CPU, the manager). In the latter case, don't pass
##      the same seed to every worker.
##    - *Warnings* are useful, don't suppress them. E.g. if the optimum in an
##      optimization has not yet been reached
##    - Measuring *run time* (system vs wall clock) is difficult to do properly.
##      It depends on architecture, programming style, compiler, current
##      workload etc.; use benchmarks to compare against.


### Simple manipulations; numbers and vectors ##################################

## Simple manipulations
1/2
1/0 # in R, Inf and -Inf exist and R can often deal with them correctly
1/-0
0/0 # ... also NaN = 'not a number' is available; 0/0, 0*Inf, Inf-Inf lead to NaN
x <- 0/0 # store the result in 'x'
class(x) # the class/type of 'x'; => NaN is still of mode 'numeric'
class(Inf) # Inf is of mode 'numeric' (although mathematically not a number); helpful in optimizations
class(NULL) # the R NULL object (a reserved object often used as default argument of functions or to represent special cases or undefined return values)
is.null(NULL)

## Vectors (data structure which contains objects of the same mode)
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

## Understanding all.equal()
all.equal # only see arguments 'target', 'current', no code (S3 method)
methods(all.equal) # => all.equal.numeric applies here
str(all.equal.numeric) # => 'tolerance' argument; str() prints the *str*ucture of an object
sqrt(.Machine$double.eps) # default tolerance 1.490116e-08 > 1e-8
all.equal(1e-7, 0) # reports for arguments (target, current) the error mean(abs(current - target)) / mean(abs(target)); here: relative error |0-1e-7|/1e-7 = 1
all.equal(0, 1e-7) # relative error "|1e-7-0|/0" => absolute error is used instead; see ?all.equal

## Numerically not exactly the same
x <- var(1:4)
y <- sd(1:4)^2
all.equal(x, y) # numerical equality
x == y # ... but not exactly
x - y # numerically not 0
## See also https://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-doesn_0027t-R-think-these-numbers-are-equal_003f

## Floating point numbers
1.7e+308 # = 1.7 * 10^308 (scientific notation)
1.8e+308 # => there is a largest positive (floating point, not real) number
2.48e-324 # near 0
2.47e-324 # truncation to 0 => there is a smallest positive (floating point) number
1 + 2.22e-16 # near 1
1 + 2.22e-16 == 1 # ... but actually isn't (yet)
1 + 1.11e-16 == 1 # indistinguishable from 1 (= 1 + 2.22e-16/2)
## Note: The grid near 0 is much finer than near 1

## Remark: These phenomena are better understood with the IEEE 754 standard
##         for floating point arithmetic.
str(.Machine) # lists important specifications of floating point numbers
.Machine$double.eps # smallest positive number x s.t. 1 + x != 1
.Machine$double.xmin # smallest normalized number > 0
.Machine$double.xmax # largest normalized number

## Watch out
n <- 0
1:n # not the empty sequence but c(1, 0); caution in 'for loops': for(i in 1:n) ...!
seq_len(n) # better: => empty sequence
seq_along(c(3, 4, 2)) # 1:3; helpful to 'go along' objects

## Watch out
1:3-1 # ':' has higher priority; note also: the '-1' is recycled to the length of 1:3
1:(3-1)

## Some functions (if functions exist, use them!)
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
z[(!is.na(z)) &  is.finite(z) &  z >= 2] # pick out all finite numbers >= 2
z[(!is.na(z)) && is.finite(z) && z >= 2] # watch out; used to fail; R >= 3.6.0: z[TRUE] => z

## Matching in indices or names
match(1:4, table = 3:5) # positions of elements of first in second vector (or NA)
1:4 %in% 3:5 # logical vector
which(1:4 %in% 3:5) # positions of TRUE in logical vector
which(3:5 %in% 1:4) # close to match() but without NAs
na.omit(match(1:4, table = 3:5)) # same (apart from attributes)
## Note: na.fill() from package 'zoo' is helpful in filling NAs in time series

## Character vectors
x <- "apple"
y <- "orange"
(z <- paste(x, y)) # paste together; use sep = "" or paste0() to paste without space
paste(1:3, c(x, y), sep = " - ") # recycling ("apple" appears again)

## Named vectors
(x <- c("a" = 3, "b" = 2)) # named vector of class 'numeric'
x["b"] # indexing elements by name (useful!)
x[["b"]] # drop the name


### Arrays and matrices ########################################################

## Matrices
(A  <- matrix(1:12, ncol = 4)) # watch out, R operates on/fills by columns (column-major order)
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
diag(A) # diagonal of A
diag(3) # identity (3, 3)-matrix
(D <- diag(1:3)) # diagonal matrix with elements 1, 2, 3
D %*% B # matrix multiplication
B * B # Hadamard product, i.e., element-wise product

## Build a covariance matrix, its correlation matrix and inverse
L <- matrix(c(2, 0, 0,
              6, 1, 0,
             -8, 5, 3), ncol = 3, byrow = TRUE) # Cholesky factor of the ...
Sigma <- L %*% t(L) # ... real, symmetric, positive definite (covariance) matrix Sigma
standardize <- Vectorize(function(r, c) Sigma[r,c]/(sqrt(Sigma[r,r])*sqrt(Sigma[c,c])))
(P <- outer(1:3, 1:3, standardize)) # construct the corresponding correlation matrix
stopifnot(all.equal(P, cov2cor(Sigma))) # a faster way
P.inv <- solve(P) # compute P^{-1}; solve(A, b) solves Ax = b (system of linear equations); if b is omitted, it defaults to I, thus leading to A^{-1}
P %*% P.inv # (numerically close to) I
P.inv %*% P # (numerically close to) I
## Another useful function is Matrix::nearPD(Sigma, corr = TRUE) which finds a
## correlation matrix close to the given matrix in the Frobenius norm.

## Other useful functions
rowSums(A) # row sums
apply(A, 1, sum) # the same
colSums(A) # column sums
apply(A, 2, sum) # the same
## Note that there are also faster functions .rowSums(), .colSums()

## Matrices are only vectors with attributes to determine when to 'wrap around'
matrix(, nrow = 1e6, ncol = 1e6) # fails to allocate a *vector* of length 1e12

## Array (data structure which contains objects of the same mode)
## Special cases: vectors (1d-arrays) and matrices (2d-arrays)
arr <- array(1:24, dim = c(2,3,4),
             dimnames = list(x = c("x1", "x2"),
                             y = c("y1", "y2", "y3"),
                             z = paste("z", 1:4, sep = ""))) # (2,3,4)-array with dimensions (x,y,z)
arr # => also filled in the first dimension first, then the second, then the third
str(arr) # use str() to the *str*ucture of the object arr
arr[1,2,2] # pick out a value
(mat <- apply(arr, 1:2, FUN = sum)) # for each combination of fixed first and second variables, sum over all other dimensions


### Lists (including data frames) ##############################################

## Data frames are rectangular objects containing objects of possibly different
## type of the same length
(df <- data.frame(Year = as.factor(c(2000, 2000, 2000, 2001, 2003, 2003, 2003)), # loss year
                  Line = c("A", "A", "B", "A", "B", "B", "B"), # business line
                  Loss = c(1.2, 1.1, 0.6, 0.8, 0.4, 0.2, 0.3))) # loss in M USD, say
str(df) # => first two columns are factors
is.matrix(df) # => indeed no matrix
as.matrix(df) # coercion to a character matrix
data.matrix(df) # coercion to a numeric matrix (factor are replaced according to their level index)

## Computing maximal losses per group for two different groupings
## Version 1 (leads a table structure with all combinations of selected variables):
tapply(df[,"Loss"], df[,"Year"], max) # maximal loss per Year
tapply(df[,"Loss"], df[,c("Year", "Line")], max) # maximal loss per Year-Line combination
## Version 2 (omits NAs):
aggregate(Loss ~ Year,        data = df, FUN = max) # 'aggregate' Loss per Year with max()
aggregate(Loss ~ Year * Line, data = df, FUN = max)

## Playing with the data frame
(fctr <- factor(paste0(df$Year,".",df$Line))) # build a 'Year.Line' factor
(grouped.Loss <- split(df$Loss, f = fctr)) # group the losses according to fctr
sapply(grouped.Loss, FUN = max) # maximal loss per group
sapply(grouped.Loss, FUN = length) # number of losses per group; see more on *apply() later
(ID <- unlist(sapply(grouped.Loss, FUN = function(x) seq_len(length(x))))) # unique ID per loss group
(df. <- cbind(df, ID = ID)) # paste unique ID per loss group to 'df'
(df.wide <- reshape(df., # reshape from 'long' to 'wide' format
                    v.names = "Loss", # variables to be displayed as entries in 2nd dimension
                    idvar = c("Year", "Line"), # variables to be kept in long format
                    timevar = "ID", # unique ID => wide columns have headings of the form <Loss.ID>
                    direction = "wide"))


## Lists
is.list(df) # => data frames are indeed just lists

## Lists are the most general data structures in R in the sense that they
## can contain pretty much everything, e.g., lists themselves or functions
## or both... (and of different lengths)
(L <- list(group = LETTERS[1:4], value = 1:2, sublist = list(10, function(x) x+1)))

## Extract elements from a list
## Version 1:
L[[1]] # get first element of the list
L[[3]][[1]] # get first element of the sub-list
## Version 2: use '$'
L$group
L$sublist[[1]]
## Version 3 (most readable and fail-safe): use the provided names
L[["group"]]
L[["sublist"]][[1]]

## Change a name
names(L)
names(L)[names(L) == "sublist"] <- "sub.list"
str(L)

## Watch out
L[[1]] # the first component
L[1] # the sub-list containing the first component of L
class(L[[1]]) # character
class(L[1]) # list


### Statistical functionality ##################################################

## Probability distributions (d/p/q/r*)
dexp(1.4, rate = 2) # density f(x) = 2*exp(-2*x)
pexp(1.4, rate = 2) # distribution function F(x) = 1-exp(-2*x)
qexp(0.3, rate = 2) # quantile function F^-(y) = -log(1-y)/2
rexp(4,   rate = 2) # draw random variates from Exp(2)

## Two-sample t-test
## Ex.: Is the average loss in two different business lines equal? H0: means are equal
plot(Loss ~ Line, data = df) # boxplot
(res <- t.test(Loss ~ Line, data = df)) # Two-sided Welch's test
str(res) # => class = "htest"
stopifnot(res$p.value < 0.05) # H0 rejected at 5%
t.test(subset(df, Line == "A", select = Loss), # different call
       subset(df, Line == "B", select = Loss))

## Binomial test
## Ex.: You have estimated the (1-p)-quantile of 1000 losses and recorded 18
##      exceedances over this estimate. Is the exceedance (= success) probability
##      really p? H0: p = 0.01
(res <- binom.test(18, n = 1000, p = 0.01)) # Two-sided binomial test
stopifnot(res$p.value < 0.05) # H0 rejected at 5%


### Random number generation ###################################################

.Random.seed # in a new R session, this object does not exist (until RNs are drawn)

## Generate from N(0,1)
(X <- rnorm(2)) # generate two N(0,1) random variates
str(.Random.seed) # encodes types of random number generators (RNGs) and the seed
## Note (see ?.Random.seed):
## - The first integer in .Random.seed consists of...
##   1) '03' = U(0,1) RNG (Mersenne Twister)
##   2) '4'  = N(0,1) RNG (inversion)
##   3) '10' = U{1,...,n} RNG (rejection)
##   Technical note:
##   + 3) was newly added in R 3.6.0
##   + 3) is used in sample() -> sample.int() -> C functions sample and sample2
##     -> ./src/main/names.c: do_sample and do_sample2 -> ./src/main/random.c
##     (for do_sample => uses the alias method in the non-uniform case and
##     R_unif_index in the uniform case) and ./src/main/unique.c (for do_sample2
##     => uses R_unif_index -> ./src/main/RNG.c which generates (a vector of)
##     random bits until a number < n appears (rejection algorithm); 'dn' most
##     likely stands for 'n as a double'.
## - The remaining integers denote the actual seed.
## - The default kind is the "Mersenne Twister" (which needs an integer(624)
##   as seed and the current position in this sequence, so 625 numbers).
## - If no random numbers were generated yet in an R session, .Random.seed
##   will not exist. If you start drawing random numbers without calling
##   set.seed(), .Random.seed is constructed from system time and the R process
##   number. You can also do rm(.Random.seed) and call runif(1) thereafter to
##   convince yourself that .Random.seed is newly generated.
RNGkind() # => Mersenne Twister, with inversion for N(0,1) and rejection for U{1,..,n}

## How can we make sure to obtain the same results (for *reproducibility*?)
(Y <- rnorm(2)) # => another two N(0,1) random variates
all.equal(X, Y) # obviously not equal (here: with probability 1)

## Set a 'seed' so that computations are reproducible
set.seed(271) # with set.seed() we can set the seed
X <- rnorm(2) # draw two N(0,1) random variates
set.seed(271) # set the same seed again
Y <- rnorm(2) # draw another two N(0,1) random variates
all.equal(X, Y) # => TRUE
set.seed(271)
Y <- rnorm(3)
all.equal(X, Y[1:2])

## A pseudo-random number generator which allows for easily advancing the
## seed is L'Ecuyer's combined multiple-recursive generator (CMRG); see
## MRG32k3a.c and MRG32k3a.h on http://simul.iro.umontreal.ca/rng (R's version
## only makes minor modifications to this). Let's see how we can call it from R.
RNGkind() # => Mersenne Twister, inversion is used for generating N(0,1)
RNGkind("L'Ecuyer-CMRG")
RNGkind() # => L'Ecuyer's CMRG, inversion is used for generating N(0,1)
.Random.seed # => now of length 7 (first number similarly as above and seed of length 6)
Z <- rnorm(2) # use L'Ecuyer's CMRG for generating random numbers
library(parallel) # for nextRNGStream() for advancing the seed
.Random.seed <- nextRNGStream(.Random.seed) # advance seed by 2^127
Z. <- rnorm(2) # generate from next stream => will be 'sufficiently apart' from Z
RNGkind("Mersenne-Twister") # switch back to Mersenne-Twister
RNGkind()


### Control statements #########################################################

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


### Working with additional packages ###########################################

packageDescription("qrmtools") # description of an installed package
maintainer("qrmtools") # the maintainer
citation("qrmtools") # how to cite a package

## Generate and plot data from a multivariate t distribution
set.seed(271) # for reproducibility (see below)
library(nvmix) # for rStudent()
X <- rStudent(2000, df = 4.5, scale = P) # generate data from a multivariate t_4.5 distribution
library(lattice) # for cloud()
cloud(X[,3]~X[,1]+X[,2], scales = list(col = 1, arrows = FALSE), col = 1,
      xlab = expression(italic(X[1])), ylab = expression(italic(X[2])),
      zlab = expression(italic(X[3])),
      par.settings = list(background = list(col = "#ffffff00"),
                axis.line = list(col = "transparent"), clip = list(panel = "off")))
## => not much visible; in higher dimensions even impossible ...

## Pairs plot (or scatter plot matrix)
pairs(X, gap = 0, pch = ".") # ... but we can use a pairs plot

## ... or dynamically:
library(rgl) # for plot3d()
plot3d(X, xlab = expression(italic(X)[1]), ylab = expression(italic(X)[2]),
       zlab = expression(italic(X)[3]))


### Writing functions ##########################################################

## gap = 0 is a good default for pairs(), so we can define our own scatter plot
## matrix (also with nice default labels)

##' @title Scatter Plot Matrix with Defaults
##' @param x data matrix (see ?pairs)
##' @param gap see 'gap' of pairs()
##' @param labels see 'labels' of pairs()
##' @param ... additional arguments passed to the underlying pairs()
##' @return invisible(NULL); see pairs.default
##' @author Marius Hofert
mypairs <- function(x, gap = 0, labels = paste("Variable", 1:ncol(x)), ...)
    pairs(x, gap = gap, labels = labels, ...) # ... pass through the formal arguments
## Note:
## - 'lazy evaluation': as ncol(x) is only evaluated once needed
## - one-line functions don't need the embracing '{}'
## - the last line of the function is returned by default, no need to call return()

## Calls
mypairs(X)
mypairs(X, pch = 19)


## Let's write a function for computing the mean of an Exp(lambda) analytically or
## via Monte Carlo

##' @title Function for Computing the Mean of an Exp(lambda) Distribution
##' @param lambda parameter (vector of) lambda > 0
##' @param n Monte Carlo sample size
##' @param method character string indicating the method to be used:
##'        "analytical": analytical formula 1/lambda
##'        "MC": Monte Carlo based on the sample size 'n'
##' @return computed mean(s)
##' @author Marius Hofert
##' @note vectorized in lambda
exp_mean <- function(lambda, n, method = c("analytical", "MC"))
{
    stopifnot(lambda > 0) # input check(s)
    method <- match.arg(method) # match the provided method with the two options
    switch(method, # distinguish (switch between) the two methods
           "analytical" = {
               1/lambda # vectorized in lambda (i.e., works for a vector of lambda's)
               ## Note: We did not use 'n' in this case, so although it is a
               ##       formal argument, we do not need to provide it.
           },
           "MC" = {
               if(missing(n)) # checks the formal argument 'n'
                   stop("You need to provide the Monte Carlo sample size 'n'")
               E <- rexp(n)
               sapply(lambda, function(l) mean(E/l)) # vectorized in lambda
           },
           stop("Wrong 'method'"))
}
## Note: We could have also omitted 'n' and used '...'. We would then need
##       to check in the "MC" case whether 'n' was provided. Checking can
##       be done with hasArg(n) in this case and 'n' can be obtained with
##       list(...)$n.

## Calls
(true <- exp_mean(1:4))
set.seed(271)
(sim <- exp_mean(1:4, n = 1e7, method = "MC"))
stopifnot(all.equal(true, sim, tol = 0.001))


## A word concerning efficiency

## Some function mimicking a longer computation
f <- function() {
    ## Longer computation to get 'x'
    Sys.sleep(1) # mimics a longer computation
    x <- 1
    ## Longer computation to get 'y'
    Sys.sleep(1) # mimics a longer computation
    y <- 2
    ## Return
    list(x = x, y = y)
}

## If you need both values 'x' and 'y', then do *not* do...
x <- f()$x
y <- f()$y
## ... as this calls 'f' twice. Do this instead:
res <- f() # only one function call
x <- res$x
y <- res$y


### Misc #######################################################################

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
