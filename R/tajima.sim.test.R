#
# Tajima's D with a simulated p-value (see tajima.test for p-value based on
# normal or beta distribution)
#


tajima.test2 <- function(n,khat,S){
    tmp <- 1:(n - 1)
    a1 <- sum(1/tmp)
    a2 <- sum(1/tmp^2)
    b1 <- (n + 1)/(3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
    c1 <- b1 - 1/a1
    c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    D <- (khat - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
    return(D)
}

#this is wrong at the moment!!
#Returns skewed distirbution of D from simulated data
tajima.sim.test <- function(x, nreps=1000, plot=TRUE){
    observed.D <- tajima.test(x)$D
    n <- if (is.list(x)) length(x) else dim(x)[1]
    S <- length(seg.sites(x))
    cmd <- paste('ms', n, nreps, '-s', S)
    sims <- system(cmd, intern=TRUE)
    #first two lines are junk, and throw the index off
    sims <- sims[3:length(sims)]
    hap.index <- seq(5, length(sims), by=n+4)
    get.seqs <- function(i){
      return(matrix(unlist(strsplit(sims[i:(i+n-1)],"")), nrow=n, ncol=S))
    }
    m <- lapply(hap.index, get.seqs)
    khats <- sapply(m, function(x) mean(dist(x,'manhattan')))
    Ds <- tajima.test2(n, khats, S)
    centered <- scale(Ds, scale=FALSE)
    pval <- mean(abs(centered) > abs(observed.D))
    if (plot){
      hist(Ds)
     abline(v=observed.D, col="blue")
    }
    return(list(D=observed.D, p.value=pval))
}
