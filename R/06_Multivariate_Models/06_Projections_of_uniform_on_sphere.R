## By Marius Hofert

## Orthogonal projections of the uniform distribution on the unit sphere S_{d-1}
## in d dimensions down by two dimensions are uniform in the (d-2)-ball;
## see http://www.math.univ-toulouse.fr/~letac/Archimede.pdf
## Let's check if that's the case!


### Setup ######################################################################

library(rgl)

## Sample from the uniform distribution on the unit d-sphere S_{d-1}
n <- 5e3 # sample size
d <- 4 # maximal dimension investigated
set.seed(271) # for reproducibility
Z <- matrix(rnorm(n*d), ncol = d) # generate independent N(0,1)


### 1 Projections from d = 3 to d = 2 ##########################################

S2 <- Z[,1:3]/sqrt(rowSums(Z[,1:3]^2)) # ~ U(S_2)
pairs(S2, gap = 0, pch = ".") # => does not look uniform (fine)


### 2 Projections from d = 3 to d = 1 ##########################################

plot(S2[,1], ylab = expression(S[1]), cex = 0.5) # => looks uniform in [0,1]
plot(S2[,2], ylab = expression(S[2]), cex = 0.5) # => looks uniform in [0,1]
plot(S2[,3], ylab = expression(S[3]), cex = 0.5) # => looks uniform in [0,1]


### 3 Projections from d = 4 to d = 3 ##########################################

S3 <- Z[,1:4]/sqrt(rowSums(Z[,1:4]^2)) # ~ U(S_3)
colnames(S3) <- paste0("S", 1:4)
plot3d(S3[,c(2,3,4)]) # => looks uniform in the ball in 3d but difficult to judge
S3. <- S3[0.2 < S3[,2] & S3[,2] <= 0.3, c(2,3,4)] # pick out a slice (0.2, 0.3]
plot(S3.) # slice looks uniform
plot3d(S3[,c(1,3,4)])
plot3d(S3[,c(1,2,4)])
plot3d(S3[,c(1,2,3)])

## Note: Although not covered by the above result, it *seems* that projections
##       from U(S_3) to IR^3 are also uniform in the 3-ball. More investigations
##       necessary.


### 4 Projections from d = 4 to d = 2 ##########################################

pairs(S3, gap = 0, pch = ".", # => looks (pairwise) uniform in the unit disk
      labels = as.expression(sapply(1:4, function(j) bquote(S[.(j)]))))
