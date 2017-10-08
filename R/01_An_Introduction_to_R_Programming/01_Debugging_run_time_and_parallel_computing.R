## By Marius Hofert

## In this example we learn how to do multi-node and multi-core computing in R.
## To this end we study the error when computing probabilities of the form
## P(a < U <= b) for a d-dimensional random vector U with independent U(0,1)
## distributed components. Furthermore, we quickly demonstrate how debugging
## and run time measurement can be done in R.


### In general: Watch out for numerical issues #################################

## How to evaluate choose(500, 200)?
choose(500, 200)
factorial(500)/(factorial(200)*factorial(300)) # too large values
n <- 500
x <- 1:n
y <- sapply(1:n, factorial)
plot(x, y, type = "l", log = "y", xlab = "n", ylab = "n!") # ... always look at plots
str(.Machine) # => double.xmax
summary(y) # => beyond double.xmax, R uses Inf here

## A numerical trick
log(factorial(200)) # obviously, same problem here ('mathematical composition' not doable)
lc <- lfactorial(500) - (lfactorial(200) + lfactorial(300)) # work with 'proper' logs
c <- exp(lc) # due to the '-/+', the result is of a size representable in computer arithmetic
c. <- exp(sum(log(1:500)) - sum(log(1:200)) - sum(log(1:300))) # we can mimic this trick
stopifnot(all.equal(choose(500, 200), c, c.)) # note: choose(500, 200) uses the same trick


### 1 Auxiliary functions ######################################################

##' @title Relative error of P(a < U <= b) under U ~ Pi for a randomly drawn (a,b]
##' @param d dimension
##' @return relative error when computing P(a < U <= b) based on U ~ Pi
##'         for a randomly constructed d-dimensional hypercube (a, b]
##' @author Marius Hofert
rel_err_P_under_Pi_random_vol <- function(d) # dimension
{
    ## Randomly generate the lower-left and upper-right corner of the cube
    pt <- matrix(runif(d * 2), ncol = 2) # 2 points in (0,1)^d (one in each col)
    cube <- cbind(a = pmin(pt[,1], pt[,2]),
                  b = pmax(pt[,1], pt[,2])) # choose a, b such that a <= b

    ## Compute P(a < U <= b) directly
    len <- apply(cube, 1, diff) # cube edge lengths
    prob <- prod(len) # or: prod(len) = exp(log(prod(len))) = exp(sum(log(len)))

    ## Compute P(a < U <= b) by evaluating Pi at all corners of (a,b],
    ## and using the 'checkerboard' way of determining the signs.

    ## For iterating over all 2^d corners, we iterate over the binary
    ## representation of m = 0:(2^d-1); see, for example, sfsmisc::digitsBase().
    D <- 2^d
    m <- 0:(D-1)
    II <- matrix(0, nrow = D, ncol = d)
    for (i in d:1L) {
        II[,i] <- m %% 2L + 1L
        if (i > 1) m <- m %/% 2L
    }
    Sign <- c(1,-1)[1L + (- rowSums(II)) %% 2] # (checkerboard) signs
    corners <- array(cube[cbind(c(col(II)), c(II))], dim = dim(II)) # (2^d, d)-matrix of corners
    C <- apply(corners, 1, prod) # evaluate the U(0,1)^d distribution function
    prob. <- sum(Sign * C) # add the checkerboard-sign-adjusted values together

    ## Compute relative error
    ## Note: |(x-y)/y| = eps => |x-y| = eps*|y| => (+/-)(x-y) = (+/-)*eps*y
    ##       => (x-y) = (+/-)*eps*y => x = (1 +/- eps) * y
    abs((prob. - prob) / prob)
}


### 2 Helpful features for software development in R ###########################

if(FALSE)
{
    ## Debugging
    ## Executing rel_err_P_under_Pi_random_vol() line-wise and interacting with it
    debug(rel_err_P_under_Pi_random_vol) # put rel_err_P_under_Pi_random_vol() in debug mode
    rel_err_P_under_Pi_random_vol(5) # call it
    undebug(rel_err_P_under_Pi_random_vol) # disable debug mode for rel_err_P_under_Pi_random_vol()

    ## Measuring run time
    ## 1) For smaller pieces of source code, use 'microbenchmark' (more reliable than
    ##    a simple system.time(replicate(1000), <function>) or even just
    ##    system.time(<function>)). microbenchmark() estimates the run time
    ##    in one of s, ms (10^{-3}s), us (10^{-6}s) or ns (10^{-9}s).
    system.time(replicate(1000, rel_err_P_under_Pi_random_vol(5))) # a simple version
    library(microbenchmark) # more sophisticated
    mbm <- microbenchmark(rel_err_P_under_Pi_random_vol(5), times = 1000)
    mbm # in ms; see also the argument 'unit'
    boxplot(mbm, unit = "ms", ylab = "Time in ms in log-scale") # boxplot in ms; note: the output is still a bit random
    ## 2) For larger chunks, 'profile' your code with Rprof(). This is also useful
    ##    to figure out where run time is 'lost' and to improve the run time of a
    ##    function.
    Rprof(profiling <- tempfile(), line.profiling = TRUE) # enable profiling
    err <- rel_err_P_under_Pi_random_vol(20)
    Rprof(NULL) # disable profiling
    (profile <- summaryRprof(profiling, lines = "both")) # get a summary
    str(profile)
    profile$by.self # profiling output for each relevant function (as 'self' and 'total', in absolute time and %)
}


### 3 Parallel computing #######################################################

## Setup
library(parallel) # for parLapply() (multi-node) and mclapply() (multi-core) functionality
d <- 8:16 # dimensions in which we compute the volumes
B <- 100 # number of replications (of measuring the relative error)

## Determine the number of cores/nodes to use
(num.cores <- detectCores()) # detect the number of available cores
num.workers <- num.cores

## Parallel computation
## Trade-off: Iterate over the B repetitions or the dimensions d first?
## 1) If we iterate over B first, all workers have (roughly) the same run time
##    (all have to work on the dimensions d)
## 2) If we iterate over d first, we can use the same seed on each worker
##    (guarantees reproducibility), but some workers are way faster
##    (biggest d takes longest), so we would require load balancing
## We use 1) here.
fun <- function(x, fun, method = c("mc", "sc")) {
    if(method == "mc") mclapply(x, fun,
                                mc.cores = num.workers,
                                mc.preschedule = FALSE) # with load balancing
    else lapply(x, fun)
}

## Sequentially
set.seed(271) # set seed here (sequential order guarantees reproducibility)
system.time(res. <- fun(seq_len(B), function(b) { # multi-core (with load balancing)
                vapply(d, function(d.) rel_err_P_under_Pi_random_vol(d.), NA_real_)
            }, method = "sc"))

## In parallel
set.seed(271) # set seed here (but not guaranteed that workers proceed in same order again)
system.time(res <- fun(seq_len(B), function(b) { # multi-core (with load balancing)
                vapply(d, function(d.) rel_err_P_under_Pi_random_vol(d.), NA_real_)
            }, method = "mc"))
## ... but method = "mc" is *not* num.cores-times faster! (due to overhead)

## Boxplot
res.mat <- matrix(unlist(res), nrow = B, ncol = length(d), byrow = TRUE) # (B, length(d))-matrix
df <- data.frame(rel.err = as.vector(res.mat), d = factor(rep(d, each = B))) # data.frame
sum(df$rel.err == 0, na.rm = TRUE) # number of relative errors equal to 0
df. <- df[df$rel.err > 0,] # remove 0s (to plot y-axis in log-scale)
boxplot(rel.err~d, data = df., boxwex = 0.3, outline = FALSE, log = "y",
        xlab = "d", ylab = expression("Relative error of P(a < U"<="b) under independence"))
mtext(substitute(B == B.~"replications", list(B. = B)), side = 4, line = 0.5, adj = 0)

## Remark:
## - Up to simulation error, the error gets larger in higher dimensions (more corners).
## - In high dimensions, this can lead to P(a < U <= b) being outside [0,1].
## - Especially a problem for copulas evaluated by MC
##   (like the multivariate normal or t distributions in >= 3 dimensions).
## - The run time is *not* "sequential run time" / num.workers (due to overhead)