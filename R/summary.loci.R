## summary.loci.R (2009-11-30)

##   Print and Summaries of Loci Objects

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

getPloidy <- function(x)
    unlist(lapply(x[, attr(x, "locicol"), drop = FALSE],
                  function(x) sum(charToRaw(levels(x)[1]) == 47) + 1))

getAlleles <- function(x)
    lapply(x[, attr(x, "locicol"), drop = FALSE],
           function(x) unique(unlist(strsplit(levels(x), "/"))))

getGenotypes <- function(x)
    lapply(x[, attr(x, "locicol"), drop = FALSE], levels)

print.loci <- function(x, details = FALSE, ...)
{
    if (details) print.data.frame(x) else {
        n <- dim(x)
        nloci <- length(attr(x, "locicol"))
        cat("Allelic data frame:", n[1])
        if (n[1] == 1) cat(" individual\n") else cat(" individuals\n")
        cat("                   ", nloci)
        if (nloci == 1) cat(" locus\n") else cat(" loci\n")
        nav <- n[2] - nloci
        if (nav) {
            cat("                   ", nav)
            if (nav == 1) cat(" additional variable\n")
            else cat(" additional variables\n")
        }
    }
}

summary.loci <- function(object, ...)
{
    ## modified tabulate to allow NA's
    tabul <- function(x, n)
        .C("R_tabulate", x, length(x), n, integer(n),
           NAOK = TRUE, DUP = FALSE, PACKAGE = "base")[[4]]

    L <- attr(object, "locicol")
    ans <- vector("list", length(L))
    names(ans) <- names(object[L])
    ii <- 1L
    for (i in L) {
        geno <- levels(object[, i])
        alle <- strsplit(geno, "/")
        unialle <- sort(unique(unlist(alle)))
        l <- tabul(object[, i], length(geno))
        names(l) <- geno
        tab <- matrix(0, length(unialle), length(geno),
                      dimnames = list(unialle, geno))
        for (j in seq_along(alle))
            for (k in alle[[j]]) tab[k, j] <- tab[k, j] + 1
        ans[[ii]] <- list(genotype = l, allele = drop(tab %*% l))
        ii <- ii + 1L
    }
    class(ans) <- "summary.loci"
    ans
}

print.summary.loci <- function(x, ...)
{
    nms <- names(x)
    for (i in 1:length(x)) {
        cat("Locus", nms[i], ":\n")
        cat("-- Genotype frequencies:\n")
        print(x[[i]][[1]])
        cat("-- Allele frequencies:\n")
        print(x[[i]][[2]])
        cat("\n")
    }
}

"[.loci" <- function(x, i, j, drop = TRUE)
{
    oc <- oldClass(x)
    loci.nms <- names(x)[attr(x, "locicol")]
    class(x) <- "data.frame"
    x <- NextMethod("[")
    ## restore the class and the "locicol" attribute only if there
    ## is at least 2 col *and* at least one loci returned:
    if (length(x) > 1) {
        locicol <- match(loci.nms, names(x))
        locicol <- locicol[!is.na(locicol)]
        if (length(locicol)) {
            attr(x, "locicol") <- locicol
            class(x) <- oc
        }
    }
    x
}

rbind.loci <- function(...) NextMethod("rbind")
## No need to drop the class cause it calls rbind.data.frame.
## The successive bindings eventually reorders the columns to
## agree with the 1st data frame, AND its "locicol" attribute
## is kept, so no need to change anything.

cbind.loci <- function(...)
{
    x <- list(...)
    n <- length(x)
    if (n == 1) return(x[[1]])
    NC <- unlist(lapply(x, ncol))
    LOCICOL <- lapply(x, attr, "locicol")
    offset.col <- cumsum(NC)
    for (i in 2:n) LOCICOL[[i]] <- LOCICOL[[i]] + offset.col[i - 1]
    for (i in 2:n) x[[1]] <- cbind.data.frame(x[[1]], x[[i]])
    x <- x[[1]]
    ## need to restore the class, so can't use NextMethod()
    class(x) <- c("loci", "data.frame")
    attr(x, "locicol") <- unlist(LOCICOL)
    x
}
