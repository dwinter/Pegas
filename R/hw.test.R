## hw.test.R (2009-05-10)

##   Test of Hardy--Weinberg Equilibrium

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

expand.genotype <- function(n, alleles = NULL, ploidy = 2, matrix = FALSE)
{
    if (!is.null(alleles)) n <- length(alleles)
    if (ploidy == 2) {
        unan <- 1:n
        ans <- matrix(c(unan, unan), n, 2)
        ans <- rbind(ans, t(combn(unan, 2)))
    }
    if (ploidy == 4) {
        ans <- matrix(0, 0, 4)
        for (a1 in 1:n)
            for (a2 in a1:n)
                for (a3 in a2:n)
                    for (a4 in a3:n)
                        ans <- rbind(ans, c(a1, a2, a3, a4))
    }
    if (is.character(alleles))
        ans <- matrix(alleles[ans], ncol = ploidy)
    if (!matrix)
        ans <- apply(ans, 1, paste, collapse = "/")
    ans
}

hw.test <- function(x, B = 1000)
{
    multinomialCoef <- function(x) {
        n <- nrow(x)
        numerator <- factorial(ncol(x))
        coef <- numeric(n)
        for (i in 1:n) {
            y <- tabulate(x[i, ])
            y <- y[y > 0]
            coef[i] <- numerator/prod(factorial(y))
        }
        coef
    }
    test.polyploid <- function(x, ploidy) {
        nms <- names(x$allele)
        Nall <- length(nms)
        all.prop <- x$allele/sum(x$allele)
        all.geno <- expand.genotype(Nall, ploidy = ploidy, matrix = TRUE)
        Ngeno <- nrow(all.geno)
        ## expected proportions:
        E <- numeric(Ngeno)
        for (i in 1:Ngeno) E[i] <- prod(all.prop[all.geno[i, ]])
        E <- multinomialCoef(all.geno) * E
        names(E) <- expand.genotype(alleles = nms, ploidy = ploidy)
        O <- E
        O[] <- 0
        O[names(x$genotype)] <- x$genotype
        E <- sum(O) * E
        chi2 <- sum((O - E)^2/E)
        DF <- length(E) - Nall + 1
        c(chi2, DF, 1 - pchisq(chi2, DF))
    }
    test.diploid <- function(x) {
        nms <- names(x$allele)
        all.prop <- x$allele/sum(x$allele)
        E <- all.prop %o% all.prop
        ## expected proportions (homozygotes first):
        E <- c(diag(E), 2 * E[col(E) < row(E)])
        tmp <- combn(nms, 2)
        names(E) <- c(paste(nms, nms, sep="/"),
                      paste(tmp[1, ], tmp[2, ], sep="/"))
        O <- E
        O[] <- 0
        O[names(x$genotype)] <- x$genotype
        E <- sum(O) * E
        chi2 <- sum((O - E)^2/E)
        DF <- length(E) - length(nms) + 1
        c(chi2, DF, 1 - pchisq(chi2, DF))
    }
    y <- summary.loci(x)
    plo <- getPloidy(x)[1]
    if (plo == 2) {
        ans <- lapply(y, test.diploid)
    } else {
        ans <- lapply(y, test.polyploid, ploidy = plo)
        if (B) {
            warning("no Monte Carlo test available for polyploids")
            B <- 0
        }
    }
    ans <- matrix(unlist(ans), ncol = 3, byrow = TRUE)
    dimnames(ans) <- list(names(y), c("chi^2", "df", "Pr(chi^2 >)"))
    if (B) {
        LN2 <- log(2)
        test.mc <- function(x) {
            ## accept only diploids for the moment
            n <- sum(x$genotype) # Nb of individuals
            Nall <- length(x$allele) # Nb of alleles
            all.geno <- expand.genotype(Nall)
            Ngeno <- length(all.geno)
            p <- numeric(B)
            ## the constant term below is dropped:
            ## lfactorial(n) - lfactorial(2*n) + sum(lfactorial(x$allele))
            ma <- unlist(mapply(rep, 1:Nall, x$allele))
            for (k in 1:B) {
                m <- sample(ma)
                dim(m) <- c(n, 2)
                ## check the order of both alleles:
                o <- m[, 1] > m[, 2]
                m[o, 1:2] <- m[o, 2:1]
                s <- factor(paste(m[, 1], m[, 2], sep = "/"), levels = all.geno)
                f <- tabulate(s, Ngeno)
                names(f) <- all.geno
                ## we know the homozygotes are the 1st Nall genotypes:
                #h <- unlist(lapply(strsplit(names(f), "/"),
                #                   function(x) length(unique(x)))) == 1
                p[k] <- -sum(lfactorial(f)) + LN2*sum(f[-(1:Nall)])
            }
            ## find the homozygotes:
            h <- unlist(lapply(strsplit(names(x$genotype), "/"),
                               function(x) length(unique(x)))) == 1
            p0 <- -sum(lfactorial(x$genotype)) + LN2*sum(x$genotype[!h])
            sum(p <= p0)/B
        }
        ans <- cbind(ans, "Pr.exact" = unlist(lapply(y, test.mc)))
    }
    ans
}

## version avec rhyper(): (Huber et al. 2006, Biometrics)
## (doesn't work)
#for (k in 1:B) {
#    N <- n
#    f <- x$allele
#    ff <- matrix(0, m, m)
#    for (i in 1:m) {
#        a <- rhyper(1, N, N, f[i]) #HG(2 * N, N, f[i])
#        y <- rhyper(1, a, N - a, f[i] - a)
#        ff[i, i] <- y #HG(N, a, f[i] - a)
#        N <- N - (f[i] - y)
#        f[i] <- f[i] - 2*y
#        b <- 2*N - f[i]
#        if (i != m) {
#            for (j in (i + 1):m) {
#                ff[i, j] <- rhyper(1, f[i], b - f[i], f[j]) #HG(b, f[i], f[j])
#                b <- b - f[j]
#                f[j] <- f[j] - ff[i, j]
#                f[i] <- f[i] - ff[i, j]
#            }
#        }
#    }
#    p[k] <- tmp - sum(lfactorial(ff)) + sum(diag(ff)) * LN2
#}
