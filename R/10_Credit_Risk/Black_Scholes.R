BlackScholes <- function(t, S, r, sigma, K, T, type=c("call","put"))
{
    d1 <- (log(S/K) + (r+sigma^2/2)*(T-t))/(sigma*sqrt(T-t))
    d2 <- d1-sigma*sqrt(T-t)
    type <- match.arg(type)
    switch(type,
           "call" = {
               S*pnorm(d1)-K*exp(-r*(T-t))*pnorm(d2)
           },
           "put" =  {
               S*(pnorm(d1)-1)+K*exp(-r*(T-t))*(1-pnorm(d2))
           },
           stop("Wrong type"))
}

BlackScholesGreeks <- function(t, S, r, sigma, K, T)
{
    d1 <- (log(S/K) + (r+sigma^2/2)*(T-t))/(sigma*sqrt(T-t))
    d2 <- d1-sigma*sqrt(T-t)
    cbind( ## First-order
        delta = pnorm(d1),
        theta = -(S*dnorm(d1)*sigma)/(2*sqrt(T-t)) - r*K*exp(-r*(T-t))*pnorm(d2),
        rho = K*(T-t)*pnorm(d2)*exp(-r*(T-t)),
        vega = S*dnorm(d1)*sqrt(T-t),
        ## Second-order
        gamma = dnorm(d1)/(S*sigma*sqrt(T-t)),
        vanna = -dnorm(d1)*d2/sigma,
        vomma = S*dnorm(d1)*sqrt(T-t)*d1*d2/sigma)
}
