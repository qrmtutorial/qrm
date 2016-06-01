## By Marius Hofert

## In this example we learn how to do multi-node and multi-core computing in R.
## To this end we study the error (in log-scale) when computing probabilities
## of the form P(a < U <= b) for a d-dimensional random vector U with
## independent U(0,1) distributed components. Furthermore, we quickly
## demonstrate how debugging and run time measurement can be done in R.


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

##' @title Relative Error in log-Scale when Computing the Probability
##'        P(a < U <= b) for a Randomly Chosen d-dimensional Volume (a,b] and
##'        U ~ U(0,1)^d
##' @param d dimension
##' @return Relative error computed in logarithmic scale
##' @author Marius Hofert
err_random_vol <- function(d)
{
    ## Randomly generate the lower-left and upper-right corner of the cube
    pt <- matrix(runif(d*2), ncol=2) # 2 points in (0,1)^d (one in each col)
    cube <- cbind(a=apply(pt, 1, min), b=apply(pt, 1, max)) # choose a, b such that a <= b

    ## Compute the logarithm of the probability P(a < U <= b) directly
    len <- apply(cube, 1, diff) # cube edge lengths
    log.prob <- sum(log(len)) # prod(len) = exp(log(prod(len))) = exp(sum(log(len)))

    ## Compute the logarithm of the probability P(a < U <= b) by evaluating
    ## the distribution function of U (simply the product) at all corners of
    ## (a,b], using the 'checkerboard' way of determining the signs and adding
    ## everything together. For iterating over all 2^d corners, we iterate over
    ## the binary representation of m=0:(2^d-1); see, for example, sfsmisc::digitsBase().
    D <- 2^d
    m <- 0:(D-1)
    II <- matrix(0, nrow=D, ncol=d)
    for (i in d:1L) {
        II[,i] <- m %% 2L + 1L
        if (i > 1) m <- m %/% 2L
    }
    Sign <- c(1,-1)[1L + (- rowSums(II)) %% 2] # signs
    corners <- array(cube[cbind(c(col(II)), c(II))], dim=dim(II)) # (2^d, d)-matrix containing all corners
    C <- apply(corners, 1, prod) # evaluate the U(0,1)^d distribution function
    log.prob. <- log(sum(Sign * C)) # add the checkerboard-sign-adjusted values together and compute log()

    ## Compute relative error in log scale
    ## Note: |(x-y)/y| = eps => |x-y| = eps*|y| => (+/-)(x-y) = (+/-)*eps*y
    ##       => (x-y) = (+/-)*eps*y => x = (1 +/- eps) * y
    abs((log.prob. - log.prob) / log.prob)
}


### 2 Helpful features for software development in R ###########################

if(FALSE)
{
    ## Debugging
    ## Executing err_random_vol() line-wise and interacting with it
    debug(err_random_vol) # put err_random_vol() in debug mode
    err_random_vol(5) # call it
    undebug(err_random_vol) # disable debug mode for err_random_vol()

    ## Measuring run time
    ## 1) For smaller pieces of source code, use 'microbenchmark' (more reliable than
    ##    a simple system.time(replicate(1000), <function>) or even just
    ##    system.time(<function>)). microbenchmark() estimates the run time
    ##    in one of s, ms (10^{-3}s), us (10^{-6}s) or ns (10^{-9}s).
    system.time(replicate(1000, err_random_vol(5))) # a simple version
    library(microbenchmark) # more sophisticated
    mbm <- microbenchmark(err_random_vol(5), times=1000)
    mbm # in ms; see also the argument 'unit'
    boxplot(mbm, unit="ms", ylab="Time in ms in log-scale") # boxplot in ms; note: the output is still a bit random
    ## 2) For larger chunks, 'profile' your code with Rprof(). This is also useful
    ##    to figure out where run time is 'lost' and to improve the run time of a
    ##    function.
    Rprof(profiling <- tempfile(), line.profiling=TRUE) # enable profiling
    err <- err_random_vol(20)
    Rprof(NULL) # disable profiling
    (profile <- summaryRprof(profiling, lines="both")) # get a summary
    str(profile)
    profile$by.self # profiling output for each relevant function (as 'self' and 'total', in absolute time and %)
}


### 3 Parallel computing #######################################################

## Setup
library(parallel) # for parLapply() (multi-node) and mclapply() (multi-core) functionality
d <- 8:15 # dimensions in which we compute the volumes
B <- 100 # number of replications (of measuring the relative error)

## Determine the number of cores/nodes to use
(num.cores <- detectCores()) # detect the number of available cores
(num.workers <- if(num.cores >= 16) 12 else if(num.cores >= 8) 6 else
                if(num.cores >= 4) 3 else if(num.cores >=  2) 2 else
                stop("Only one core/node was detected; cannot test parallel functionality."))

## Parallel computation
## Trade-off: Shall we iterate over the B repetitions or our chosen dimensions d?
## - If we iterate over B first, all workers have (roughly) the same run time
## - If we iterate over d first, we can use the same seed on each worker
##   (guarantees reproducibility), but some workers are significantly faster
##   done than others (mostly the one dealing with the largest d is occupied then).
## => Iterate over B first. We still set the seed (although there's no guarantee that
##    the workers do their work in the same order if the computations are repeated).
multi.core <- if(num.cores > 1) TRUE else NA # NA = serial; TRUE = multi-core; FALSE = multi-node
if(is.na(multi.core)) { # single-core case

    set.seed(271) # set seed
    system.time(res <- lapply(seq_len(B), function(b) { # single-core
        vapply(d, function(d.) err_random_vol(d.), NA_real_)
    }))

} else if(multi.core) { # multi-core case

    set.seed(271) # set seed
    system.time(res <- mclapply(seq_len(B), function(b) { # multi-core (with load balancing)
        vapply(d, function(d.) err_random_vol(d.), NA_real_)
    }, mc.cores=num.workers, mc.preschedule=FALSE)) # load balancing

} else { # multi-node case

    cluster <- makeCluster(num.workers, type="PSOCK") # create cluster (type "MPI" is preferred if available)
    clusterExport(cluster, varlist=c("d", "err_random_vol")) # export variables to cluster
    set.seed(271) # set seed
    system.time(res <- parLapplyLB(cluster, seq_len(B), function(b) { # multi-node (with load balancing)
                    vapply(d, function(d.) err_random_vol(d.), NA_real_)}))
    stopCluster(cluster) # stop cluster

}

## Extract result
str(res) # B-list with length(d) values each
res.mat <- matrix(unlist(res), nrow=B, ncol=length(d), byrow=TRUE) # (B, length(d))-matrix


### 4 Analysis #################################################################

## Boxplot
res.df <- data.frame(rel.err=as.vector(res.mat), d=factor(rep(d, each=B))) # data.frame
sum(res.df$rel.err == 0, na.rm=TRUE) # number of relative errors equal to 0
res.df. <- res.df[res.df$rel.err > 0,] # remove 0s (to plot y-axis in log-scale)
boxplot(rel.err~d, data=res.df., boxwex=0.3, outline=FALSE, log="y",
        xlab="d", ylab="Simulated relative error (in log-scale)",
        main=expression(bold("Computing probabilities of"~U(0,1)^d~"of randomly chosen cubes")))
mtext(substitute(B==B.~"replications", list(B.=B)), side=4, line=0.5, adj=0)

## Result:
## Of course, up to Monte Carlo (MC) error the error gets larger in larger
## dimensions (there are more corners). In high dimensions, this can lead to
## the fact that the resulting probability of falling in a cube can be
## outside [0,1]. This is especially a problem for distributions other than
## U(0,1)^d as their evaluation might be numerically more demanding, some can
## even only be evaluated by MC themselves (like the multivariate normal or t
## distributions in >= 3 dimensions).