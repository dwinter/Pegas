## conversion.R (2010-11-10)

##   Conversion Among Allelic Data Classes

## Copyright 2009-2010 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

## as.genind <- function(x)
loci2genind <- function(x)
{
    ipop <- which(names(x) == "population")
    pop <- if (length(ipop)) x[, ipop] else NULL
    df2genind(as.matrix(x[, attr(x, "locicol")]), sep = "/", pop = pop)
}

as.loci <- function(x, ...) UseMethod("as.loci")

as.loci.genind <- function(x, ...)
{
    obj <- genind2df(x, sep = "/")
    icol <- 1:ncol(obj)
    pop <- which(names(obj) == "pop")
    if (length(pop)) {
        names(obj)[pop] <- "population"
        icol <- icol[-pop]
    }
    for (i in icol) obj[, i] <- factor(obj[, i] )
    class(obj) <- c("loci", "data.frame")
    attr(obj, "locicol") <- icol
    obj
}

genind2loci <- function(x) as.loci.genind(x)

.check.order.alleles <- function(x)
{
    for (k in attr(x, "locicol")) {
        y <- x[, k]
        if (is.numeric(y)) { # haploid with alleles coded with numerics
            x[, k] <- factor(y)
            next
        }
        lv <- levels(y)
        if (!length(grep("/", lv))) next # if haploid
        n <- length(lv)
        if (n == 1) next
        ## works with all levels of ploidy > 1:
        a <- matrix(unlist(strsplit(lv, "/")), nrow = n, byrow = TRUE)
        a <- t(apply(a, 1, sort))
        #levels(y) <- apply(a, 1, paste, collapse = "/")
        #drop <- FALSE
        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                if (all(a[i, ] == a[j, ])) {
                    y[y == lv[j]] <- lv[i]
                    #drop <- TRUE
                    #break
                }
            }
        }
        x[, k] <- factor(y)
    }
    x
}

as.loci.data.frame <-
    function(x, allele.sep = "/", col.pop = NULL, col.loci = NULL, ...)
{
    if (is.null(col.pop)) {
        ipop <- which(tolower(names(x)) == "population")
        if (length(ipop)) col.pop <- ipop
    }
    if (is.character(col.pop))
        col.pop <- which(names(x) == col.pop)
    if (is.numeric(col.pop)) {
        names(x)[col.pop] <- "population"
        x[, col.pop] <- factor(x[, col.pop])
    }
    if (is.null(col.loci)) {
        col.loci <- 1:ncol(x)
        if (is.numeric(col.pop))
            col.loci <- col.loci[-col.pop]
    }
    if (is.character(col.loci))
        col.loci <- match(col.loci, names(x))
    if (allele.sep != "/") {
        if (allele.sep == "")
            stop("alleles within a genotype must be separated")
        for (i in 1:col.loci)
            levels(x[, i]) <- gsub(allele.sep, "/", levels(x[, i]))
    }
    class(x) <- c("loci", "data.frame")
    attr(x, "locicol") <- col.loci
    .check.order.alleles(x)
}

as.loci.factor <- function(x, allele.sep = "/", ...)
    as.loci.data.frame(data.frame(x), allele.sep = allele.sep, ...)

as.loci.character <- function(x, allele.sep = "/", ...)
    as.loci.data.frame(data.frame(factor(x)), allele.sep = allele.sep, ...)
