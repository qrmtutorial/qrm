## By Marius Hofert

## Computing and plotting tail probabilities


### Setup ######################################################################

library(copula)

set.seed(271) # due to the Monte Carlo (MC) error for evaluating Ga/t copulas


### Example 7.40 (Financial example with five stocks) ##########################

d <- 5
cop <- ellipCopula("normal", param=0.5, dim=d)
pCopula(rep(0.01, d), copula=cop) # ~= 7.44e-05 (for our seed)


### Table 7.2 (Joint tail probabilities P(U_1 > u, U_2 > u) as a plot ##########

n <- 128
u <- seq(0.95, to=0.9999, length.out=n) # quantiles u where to evaluate P(U_1 > u, U_2 > u)
rho <- c(0.75, 0.5) # correlation parameter rho
nu <- c(3, 4, 8, Inf) # degrees of freedom
len <- length(rho)*length(nu)
col <- c("maroon3", "darkorange2", "royalblue3", "black") # colors for different nu
tail.prob <- matrix(, nrow=length(u), ncol=1+len, # tail probabilities
                    dimnames=list(rownames=seq_along(u),
                                  colnames=c("u", "rho.075.nu.3", "rho.075.nu.4",
                                             "rho.075.nu.8", "rho.075.nu.Inf",
                                             "rho.05.nu.3", "rho.05.nu.4",
                                             "rho.05.nu.8", "rho.05.nu.Inf")))
tail.prob[,1] <- u
expr <- vector("expression", len) # vector of expressions
ltys <- numeric(len) # line types
cols <- numeric(len) # colors
for(i in seq_along(rho)) { # rho
    for(j in seq_along(nu)) { # degrees of freedom
        k <- length(nu)*(i-1)+j
        ## Create the copula
        cop <- ellipCopula("t", param=rho[i], df=nu[j])
        ## Evaluate P(U_1 > u, U_2 > u) = P(U_1 <= 1-u, U_2 <= 1-u)
        tail.prob[,k+1] <- pCopula(cbind(1-u, 1-u), copula=cop)
        ## Create plot information
        expr[k] <- as.expression(
            substitute(group("(",list(rho, nu),")")==group("(", list(r., n.), ")"),
                       list(r.=rho[i], n.=if(nu[j]==Inf) bquote(infinity) else nu[j])))
        ltys[k] <- length(rho)-i+1
        cols[k] <- col[j]
    }
}

## Standardize w.r.t. Gauss case
tail.prob.fact <- tail.prob # tail probabilities as factors in comparison to Gauss case
tail.prob.fact[,2:5] <- tail.prob.fact[,2:5]/tail.prob.fact[,5]
tail.prob.fact[,6:9] <- tail.prob.fact[,6:9]/tail.prob.fact[,9]

## Plot tail probabilities
matplot(tail.prob[,1], tail.prob[,-1], type="l", lty=ltys, col=cols,
        xlab="Quantile u", ylab=expression("Tail probability"~P(U[1]>u,U[2]>u)))
legend("topright", inset=0.01, bty="n", lty=ltys, col=cols, legend=expr)

## Plot standardized tail probabilities
matplot(tail.prob.fact[,1], tail.prob.fact[,-1], log="y", type="l", lty=ltys, col=cols,
        lwd=c(1,1,1,1.6,1,1,1,1), xlab="Quantile u",
        ylab=expression("Tail probability"~P(U[1]>u,U[2]>u)~"standardized by Gauss case"))
legend("topleft", inset=0.01, bty="n", lty=ltys, col=cols, legend=expr)


### Table 7.3 (Joint tail probabilities P(U_1 > u, ..., U_d > u) as a plot #####

## Note: Due to radial symmetry, we can compute probabilities in the lower tail
d <- 2:20
u <- 0.99 # quantile u where to evaluate P(U_1 > u, ..., U_d > u)
rho <- c(0.75, 0.5) # correlation parameter rho
nu <- c(3, 4, 8, Inf) # degrees of freedom
len <- length(rho)*length(nu)
col <- c("maroon3", "darkorange2", "royalblue3", "black") # colors for different nu
tail.prob <- matrix(, nrow=length(d), ncol=1+len, # tail probabilities
                    dimnames=list(rownames=seq_along(d),
                                  colnames=c("d", "rho.075.nu.3", "rho.075.nu.4",
                                             "rho.075.nu.8", "rho.075.nu.Inf",
                                             "rho.05.nu.3", "rho.05.nu.4",
                                             "rho.05.nu.8", "rho.05.nu.Inf")))
tail.prob[,1] <- d
expr <- vector("expression", len) # vector of expressions
ltys <- numeric(len) # line types
cols <- numeric(len) # colors
for(i in seq_along(rho)) { # rho
    for(j in seq_along(nu)) { # degrees of freedom
        k <- length(nu)*(i-1)+j
        tail.prob <- cbind(tail.prob, rep(NA, length(d)))
        for(l in seq_along(d)) { # dimension
            ## Create the copula
            cop <- ellipCopula("t", param=rho[i], dim=d[l], df=nu[j])
            ## Evaluate P(U_1 > u, ..., U_d > u) = P(U_1 <= 1-u, ..., U_d <= 1-u)
            tail.prob[l,k+1] <- pCopula(rep(1-u, d[l]), copula=cop)
        }
        ## Create plot information
        expr[k] <- as.expression(
            substitute(group("(", list(rho, nu), ")")==group("(", list(r., n.), ")"),
                       list(r.=rho[i], n.=if(nu[j]==Inf) bquote(infinity) else nu[j])))
        ltys[k] <- length(rho)-i+1
        cols[k] <- col[j]
    }
}

## Standardize w.r.t. Gauss case
tail.prob.fact <- tail.prob # tail probabilities as factors in comparison to Gauss case
tail.prob.fact[,2:5] <- tail.prob.fact[,2:5]/tail.prob.fact[,5]
tail.prob.fact[,6:9] <- tail.prob.fact[,6:9]/tail.prob.fact[,9]

## Plot tail probabilities
matplot(tail.prob[,1], tail.prob[,-1], type="l", lty=ltys, col=cols,
        xlab="Dimension d", ylab=expression("Tail probability"~P(U[1]>u,...,U[d]>u)))
legend("topright", inset=0.01, bty="n", lty=ltys, col=cols, legend=expr)

## Plot standardized tail probabilities
matplot(tail.prob.fact[,1], tail.prob.fact[,-1], log="y", type="l", lty=ltys, col=cols,
        lwd=c(1,1,1,1.6,1,1,1,1), xlab="Dimension d",
        ylab=expression("Tail probability"~P(U[1]>u,..,U[d]>u)~"standardized by Gauss case"))
legend("topleft", inset=0.01, bty="n", lty=ltys, col=cols, legend=expr)

## Results:
##
## 'Table 7.2' case:
## - Obviously, P(U_1 > u, U_2 > u) is decreasing in u.
## - The larger rho or the smaller nu, the larger P(U_1 > u, U_2 > u).
## - Divided by the tail probability for the Gauss copula, P(U_1 > u, U_2 > u) is...
##   + ... increasing in u (the further we are in the tail, the more pronounced the difference);
##   + ... larger the smaller nu (the more we deviate from the Gauss copula);
##   + ... larger the smaller rho.
##
## 'Table 7.3' case:
## - Obviously, P(U_1 > u, ..., U_d > u) is decreasing in d.
## - The larger rho or the smaller nu, the larger P(U_1 > u, ..., U_d > u).
## - Divided by the tail probability for the Gauss copula, P(U_1 > u, ..., U_d > u) is...
##   + ... increasing in d (the larger the dimension, the more pronounced the difference);
##   + ... larger the smaller nu (the more we deviate from the Gauss copula);
##   + ... larger the smaller rho.
