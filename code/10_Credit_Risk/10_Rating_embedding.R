## By Alexander McNeil

library(expm)


Lambda <- matrix(c(0,       0.082,  0.0063, 0,      0.0003, 0,      0,      0,      0,
                   0.0091,  0,      0.0843, 0.0049, 0.0006, 0.0002, 0.0001, 0,      0.0002,
                   0.0006,  0.0248, 0,      0.0547, 0.0057, 0.0011, 0.0003, 0,      0.0006,
                   0.00039, 0.0017, 0.0411, 0,      0.0405, 0.0755, 0.0163, 0.0002, 0.0017,
                   0.0001,  0.0005, 0.0035, 0.0552, 0,      0.0722, 0.0058, 0.0007, 0.0106,
                   0.0001,  0.0003, 0.0011, 0.0032, 0.0458, 0,      0.0581, 0.0059, 0.0385,
                   0.0001,  0.0002, 0.0002, 0.0012, 0.0038, 0.0870, 0,      0.0372, 0.1334,
                   0,       0,      0,      0,      0.004,  0.0203, 0.0938, 0,      0.3793,
                   0,       0,      0,      0,      0,      0,      0,      0,      0),
                 ncol = 9, byrow = TRUE)
nms <- c("Aaa","Aa","A","Baa","Ba","B","Caa","C","D")
dimnames(Lambda) <- list(nms, nms)
Lambda

Lambda.d <- -rowSums(Lambda)
diag(Lambda) <- Lambda.d
rowSums(Lambda)
Lambda

P <- expm(Lambda) # matrix exponential
rowSums(P)

(P2 <- P[nrow(Lambda):1, ncol(Lambda):1]) # rows and columns sorted in inverse order
(P2 <- P2[-1, ]) # without first row

## Pick C
(Cprobs <- P2[1,])
(Ccumprobs <- cumsum(Cprobs))
(thresholds <- qnorm(Ccumprobs))

## Plot
x <- seq(-5, 5, length.out = 100)
plot(x, dnorm(x), type = "l", xlab = "X", ylab = "Density")
abline(v = thresholds, col = 1:length(thresholds))

## Do all
(cumprobs <- t(apply(P2, 1, function(v) cumsum(v))))
(thresholds <- qnorm(cumprobs))
layout(matrix(1:8, nrow = 2, byrow = TRUE))
for (j in 1:nrow(thresholds)) {
    plot(x, dnorm(x), type = "l",
         xlab = "X", ylab = "density", main = rownames(thresholds)[j])
    abline(v = thresholds[j,], col = 1:length(thresholds[j,]))
}
layout(1)
