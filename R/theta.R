## theta.R (2009-10-03)

##   Population Parameter THETA

## theta.h: using homozigosity
## theta.k: using expected number of alleles
## theta.s: using segregating sites in DNA sequences
## theta.tree: using a genealogy

## Copyright 2002-2009 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

theta.h <- function(x, standard.error = FALSE)
{
    HE <- H(x, variance = TRUE)
    sdH <- HE[2]
    HE <- HE[1]
    f <- function(th) HE - th * (1 + (2 * (1 + th)) / ((2 + th) * (3 + th)))
    th <- uniroot(f, interval = c(0, 1))$root
    if (standard.error) {
        SE <- (2 + th)^2 * (2 + th)^3 * sdH /
          HE^2 * (1 + th) * ((2 + th) * (3 + th) * (4 + th) + 10 * (2 + th) + 4)
        return(c(th, SE))
    }
    else return(th)
}

theta.k <- function(x, n = NULL, k = NULL)
{
    if (is.null(n)) {
        if (!is.factor(x)) {
            if (is.numeric(x)) {
                n <- sum(x)
                k <- length(x)
            }
            else x <- factor(x)
        }
        if (is.factor(x)) { # ne pas remplacer par `else'...
            n <- length(x)
            k <- nlevels(x)
        }
    }
    f <- function(th) th * sum(1 / (th + (0:(n - 1)))) - k
    th <- uniroot(f, interval = c(1e-8, 100))$root
    return(th)
}

theta.s <- function(s, n, variance = FALSE)
{
    a1 <- sum(1 / (1:(n - 1)))
    th <- s / a1
    if (variance) {
        a2 <- sum(1 / (1:(n - 1))^2)
        var.th <- (a1^2 * s + a2 * s^2) / (a1^2 * (a1^2 + a2))
        return(c(th, var.th))
    }
    else return(th)
}

theta.tree <- function(phy, theta, fixed = FALSE, log = TRUE)
{
    ## coalescent intervals from the oldest to most recent one:
    x <- rev(diff(c(0, sort(branching.times(phy)))))
    k <- 2:length(phy$tip.label)
    K <- length(k)
    tmp <- choose(k, 2)
    tmp2 <- sum(x * tmp)
    sltmp <- sum(log(tmp))
    if (fixed) {
        res <- sltmp - K * log(theta) - tmp2/theta
        if (!log) res <- exp(res)
    } else {
        minusLogLik <- function(theta) # vectorized on 'theta'
            -(sltmp - K*log(theta) - tmp2/theta)
        gr <- function(theta) K/theta - tmp2/theta^2
        out <- nlminb(theta[1], minusLogLik, gr,
                      lower = .Machine$double.eps, upper = Inf)
        res <- list(theta = out$par, logLik = -out$objective)
    }
### alternative version based on L-BFGS-B
###out <- optim(theta[1], minusLogLik, gr, method = "L-BFGS-B",
###             lower = .Machine$double.eps, upper = Inf,
###             hessian = TRUE)
###res <- list(theta = out$par, se = sqrt(1/out$hessian[, ]),
###            logLik = -out$value)
### I prefered nlminb() because it is slightly faster and in most cases
### the hessian-based estimate of SE(theta) are not needed
    res
}
