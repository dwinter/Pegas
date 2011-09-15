## haplotype.R (2011-07-11)

##   Haplotype Extraction, Frequencies, and Networks

## Copyright 2009-2011 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

.TempletonProb <- function(j, S, b = 2, r = 1)
{
    br <- b * r
    P <- numeric(max(j))
    L_jm <- function(q, j, m) {
        jm1 <- j - 1
        qonbr <- q/br
        (2*q)^jm1 * (1 - q)^(2*m + 1) * (1 - qonbr) *
            (2 - q*(br + 1)/br)^jm1 *
                (1 - 2*q*(1 - qonbr))
    }
    for (i in seq_along(P)) {
        M <- S - i
        denom <- integrate(L_jm, 0, 1, j = i, m = M)$value
        ## eq.7 from Templeton et al. 1992:
        out <- integrate(function(q) q*L_jm(q, j = i, m = M), 0, 1)$value/denom
        P[i] <- 1 - out
    }
    cumprod(P)[j]
}

haplotype <- function(x, labels = NULL)
{
    nms.x <- deparse(substitute(x))
    if (is.list(x)) x <- as.matrix(x)
    rownames(x) <- NULL
    y <- apply(x, 1, rawToChar)
    n <- length(y)
    keep <- nhaplo <- 1L
    no <- list(1L)
    for (i in 2:n) {
        already.seen <- FALSE
        j <- 1L
        while (j <= nhaplo) {
            if (y[i] == y[keep[j]]) {
                no[[j]] <- c(no[[j]], i)
                already.seen <- TRUE
                break
            }
            j <- j + 1L
        }
        if (!already.seen) {
            keep <- c(keep, i)
            nhaplo <- nhaplo + 1L
            no[[nhaplo]] <- i
        }
    }
    obj <- x[keep, ]
    if (is.null(labels))
        labels <- as.character(as.roman(1:length(keep)))
    rownames(obj) <- labels
    class(obj) <- c("haplotype", "DNAbin")
    attr(obj, "index") <- no
    attr(obj, "from") <- nms.x
    obj
}

haploNet <- function(h, d = NULL)
{
    if (!inherits(h, "haplotype"))
        stop("'h' must be of class 'haplotype'")
    freq <- sapply(attr(h, "index"), length)
    n <- length(freq) # number of haplotypes
    link <- matrix(0, 0, 3)
    if (is.null(d)) d <- dist.dna(h, "N")
    d <- as.matrix(d)
    d[col(d) >= row(d)] <- NA # put NA's in the diag and above-diag elts
    dimnames(d) <- list(1:n, 1:n)
    step <- 1
    gr <- 1:n
    while (length(unique(gr)) > 1) {
        newLinks <- which(d == step, TRUE)
        if (length(newLinks)) {
            del <- NULL
            for (i in 1:nrow(newLinks)) {
                a <- gr[newLinks[i, 1]]
                b <- gr[newLinks[i, 2]]
                ## if both nodes are already in the
                ## same subnet, then drop this link
                if (a == b) del <- c(del, i)
                else gr[which(gr == b)] <- a
            }
            if (!is.null(del)) newLinks <- newLinks[-del, , drop=FALSE]
            link <- rbind(link, cbind(newLinks, rep(step, nrow(newLinks))))
        }
        step <- step + 1
    }
    link <- cbind(link, .TempletonProb(link[, 3], ncol(h)))
    dimnames(link) <- list(NULL, c("", "", "step", "Prob"))
    attr(link, "freq") <- freq
    attr(link, "labels") <- rownames(h)
    class(link) <- "haploNet"
    link
}

plot.haploNet <-
    function(x, size = 1, col = "black", bg = "white",
             col.link = "black", lwd = 1, lty = 1, pie = NULL,
             labels = TRUE, font = 2, cex = 1, scale.ratio = 1,
             asp = 1, legend = FALSE, fast = FALSE, ...)
{
    par(xpd = TRUE)
    link <- x[, 1:2]
    l1 <- x[, 1]
    l2 <- x[, 2]
    ld <- x[, 3] * scale.ratio
#    ld[] <- 1

    tab <- tabulate(link)
    n <- length(tab)
    xx <- yy <- angle <- theta <- r <- numeric(n)
    avlb <- !logical(length(ld))

    ## adjust 'ld' wrt the size of the symbols:
    size <- rep(size, length.out = n)
    ld <- ld + (size[l1] + size[l2])/2

    H <- vector("list", n) # the list of hierarchy of nodes...

    foo <- function(i) {
        j <- NULL # indices of the haplotypes linked to 'i'
        for (ii in 1:2) { # look at both columns
            ll <- which(link[, ii] == i & avlb)
            if (length(ll)) {
                newj <- link[ll, -ii]
                r[newj] <<- ld[ll]
                j <- c(j, newj)
                avlb[ll] <<- FALSE
            }
        }
        nlink <- length(j)
        if (nlink) {
            H[[i]] <<- j
            start <- theta[i] - angle[i]/2
            theta[j] <<-
                seq(start, by = angle[i]/nlink, length.out = nlink)
            angle[j] <<- angle[i]/nlink # *sqrt(xx[i]^2 + yy[i]^2)
            xx[j] <<- r[j] * cos(theta[j]) + xx[i]
            yy[j] <<- r[j] * sin(theta[j]) + yy[i]
            for (ii in j) foo(ii)
        }
    }

    ## start with the haplotype with the most links:
    central <- which.max(tab)
    angle[central] <- 2*pi
    foo(central)

if (!fast) {
    fCollect <- function(i) {
        ## find all nodes to move simultaneously
        ii <- H[[i]]
        if (!is.null(ii)) {
            j <<- c(j, ii)
            for (jj in ii) fCollect(jj)
        }
    }

    ## Version qui ne prend pas en compte les angles
    ## dans le calcul de l'énergie du réseau: cela semble
    ## mieux marcher mais il y a encore des line-crossings

    energy <- function(xx, yy) {
        ## check line crossings
        nlink <- length(l1)
        ## round everyhting to avoid problems:
        x0 <- round(xx[l1], 1e-6)
        y0 <- round(yy[l1], 1e-6)
        x1 <- round(xx[l2], 1e-6)
        y1 <- round(yy[l2], 1e-6)
        ## compute all the slopes and intercepts:
        beta <- (y1 - y0)/(x1 - x0)
#        browser()
        intp <- y0 - beta*x0
        for (i in 1:(nlink - 1)) {
            for (j in (i + 1):nlink) {
                ## in case they are nearly parallel:
                if (any(is.na(beta))) return(Inf)
                if (beta[i] == beta[j]) next
                ## if both lines are vertical:
                if (abs(beta[i]) == Inf && abs(beta[j]) == Inf) next
                ## find the point where both lines cross
                ## in case the 1st line is vertical...
                if (abs(beta[i]) == Inf) {
                    xi <- x0[i]
                    yi <- beta[j]*xi + intp[j]
                } else if (abs(beta[j]) == Inf) { # ... or the 2nd one
                    xi <- x0[j]
                    yi <- beta[i]*xi + intp[i]
                } else {
                    xi <- (intp[j] - intp[i])/(beta[i] - beta[j])
                    yi <- beta[i]*xi + intp[i]
                }
                ## rounding again
                xi <- round(xi, 1e-6)
                yi <- round(yi, 1e-6)

                if (x0[i] < x1[i]) {
                    if (xi <= x0[i] || xi >= x1[i]) next
                } else {
                    if (xi >= x0[i] || xi <= x1[i]) next
                }

                if (y0[i] < y1[i]) {
                    if (yi <= y0[i] || yi >= y1[i]) next
                } else {
                    if (yi >= y0[i] || yi <= y1[i]) next
                }

                ## if we reach here, the intersection point is within
                ## the 1st segment, check if it is within the 2nd one
                if (x0[j] < x1[j]) {
                    if (xi <= x0[j] || xi >= x1[j]) next
                } else {
                    if (xi >= x0[j] || xi <= x1[j]) next
                }

                if (y0[i] < y1[j]) {
                    if (yi <= y0[j] || yi >= y1[j]) next
                } else {
                    if (yi >= y0[j] || yi <= y1[j]) next
                }
#                cat("i =", i, "j =", j, "\n")
#                browser()
                return(Inf)
            }
        }
        Dang <- NULL
        for (i in NODES) {
            #cat("i =", i, "\n")
            j <- c(link[l1 == i, 2], link[l2 == i, 1])
            if (length(j) == 1) next
            alpha <- atan2(yy[j] - yy[i], xx[j] - xx[i])
            Dang <- c(Dang, diff(alpha))
            #print(Dang)
        }
        #x <- (xx[l1] + xx[l2])/2 # add the edge
        #y <- (yy[l1] + yy[l2])/2 # midpoints
        D <- dist(cbind(xx, yy))
        sum(1/c(D)^2, na.rm = TRUE)
    }

    Rotation <- function(rot, i, beta) {
        ## rot: indices of the nodes to rotate
        ## i: index of the node connected to 'rot' (= fixed rotation point)
        xx.rot <- xx[rot] - xx[i]
        yy.rot <- yy[rot] - yy[i]
        theta <- atan2(yy.rot, xx.rot) + beta
##        cat("i =", i, "rot =", rot, "theta =", theta, "beta =", beta, "\n")
        h <- sqrt(xx.rot^2 + yy.rot^2)
        new.xx[rot] <<- h*cos(theta) + xx[i]
        new.yy[rot] <<- h*sin(theta) + yy[i]
    }

    OptimizeRotation <- function(node, rot) {
        ## test the direction 1st
        inc <- pi/90
        Rotation(rot, node, inc)
        new.E <- energy(new.xx, new.yy)
        ##cat("node =", node, "rot = ", rot, "E =", E, "new.E =", new.E, "\n")
        if (new.E >= E) {
            inc <- -inc
            Rotation(rot, node, inc)
            new.E <- energy(new.xx, new.yy)
        }
        ##cat("node =", node, "rot = ", rot, "E =", E, "new.E =", new.E, "\n")
        ##plot(xx, yy, type = "n")
        ##segments(xx[l1], yy[l1], xx[l2], yy[l2], lwd = 3)
        ##text(xx, yy, as.character(as.roman(1:n)), font = 2)
        ##plot(new.xx, new.yy, type = "n")
        ##segments(new.xx[l1], new.yy[l1], new.xx[l2], new.yy[l2], lwd = 3)
        ##text(new.xx, new.yy, as.character(as.roman(1:n)), font = 2, col = "blue")

        while (new.E < E) {
            ##plot(xx, yy, type = "n")
            ##segments(xx[l1], yy[l1], xx[l2], yy[l2], lwd = 3)
            ##text(xx, yy, as.character(as.roman(1:n)), font = 2)
            ##plot(new.xx, new.yy, type = "n")
            ##segments(new.xx[l1], new.yy[l1], new.xx[l2], new.yy[l2], lwd = 3)
            ##cat("node =", node, "E =", E, "new.E =", new.E, "\n")
            xx <<- new.xx
            yy <<- new.yy
            E <<- new.E
            Rotation(rot, node, inc)
            new.E <- energy(new.xx, new.yy)
        }
    }

    NODES <- 1:n
    E <- energy(xx, yy)
    new.xx <- xx
    new.yy <- yy
###HH <<- H
    nextOnes <- NULL
    for (i in H[[central]][-1]) {
        ## collect the nodes descending from 'i':
        j <- NULL # j must be initialized before calling fCollect
        fCollect(i)
        rot <- c(i, j) # index of the nodes that will be moved
        OptimizeRotation(central, rot)
        nextOnes <- c(nextOnes, i)
    }

    while (!is.null(nextOnes)) {
        newNodes <- nextOnes
        nextOnes <- NULL
        for (i in newNodes) {
            if (is.null(H[[i]])) next
            for (j in H[[i]]) {
                fCollect(j)
                rot <- j
                OptimizeRotation(i, rot)
                nextOnes <- c(nextOnes, rot)
            }
        }
    }
}

    plot(xx, yy, type = "n", xlab = "", ylab = "",
         axes = FALSE, bty = "n", asp = asp, ...)
    segments(xx[l1], yy[l1], xx[l2], yy[l2], lwd = lwd,
             lty = lty, col = col.link)
    if (is.null(pie))
        symbols(xx, yy, circles = size/2, inches = FALSE,
                add = TRUE, fg = col, bg = bg)
    else {
        for (i in 1:n)
            floating.pie.asp(xx[i], yy[i], pie[i, ], radius = size[i]/2)
    }
    if (labels)
        text(xx, yy, attr(x, "labels"), font = font, cex = cex)
    if (legend[1]) {
        if (is.logical(legend)) {
            cat("Click where you want to draw the legend\n")
            xy <- unlist(locator(1))
        } else xy <- legend
        segments(xy[1], xy[2], xy[1] + scale.ratio, xy[2])
        text(xy[1] + scale.ratio, xy[2], " 1", adj = 0)
        if (length(unique(size)) > 1) {
            vspace <- strheight(" ")
            symbols(xy[1] + 0.5, xy[2] - 2*vspace, circles = 0.5,
                    inches = FALSE, add = TRUE)
            text(xy[1] + 0.5, xy[2] - 2*vspace, "  1", adj = 0)
        }
        if (!is.null(pie)) {
            nc <- ncol(pie)
            co <- rainbow(nc)
            TEXT <- paste(" ", colnames(pie))
            for (i in 1:nc) {
                Y <- xy[2] - 2*(i+1)*vspace
                symbols(xy[1] + 0.5, Y, circles = 0.5,
                    inches = FALSE, add = TRUE, bg = co[i])
                text(xy[1] + 0.5, Y, TEXT[i], adj = 0)
            }
        }
    }
    #browser()
}

plot.haplotype <- function(x, ...)
{
    barplot(sapply(attr(x, "index"), length), xlab = "Haplotype",
            ylab = "Number", names.arg = rownames(x), ...)
}

print.haplotype <- function(x, ...)
{
    d <- dim(x)
    DF <- sapply(attr(x, "index"), length)
    names(DF) <- rownames(x)
    cat("\nHaplotypes extracted from:", attr(x, "from"), "\n\n")
    cat("    Number of haplotypes:", d[1], "\n")
    cat("         Sequence length:", d[2], "\n\n")
    cat("Haplotype labels and frequencies:\n\n")
    print(DF)
}

as.network.haploNet <- function(x, directed = FALSE, ...)
{
    res <- network(x[, 1:2], directed = directed, ...)
    network.vertex.names(res) <- attr(x, "labels")
    res
}

as.igraph.haploNet <- function(x, directed = FALSE, use.labels = TRUE, ...)
{
    directed <- directed
    y <- x[, 1:2]
    y <-
        if (use.labels) matrix(attr(x, "labels")[y], ncol = 2)
        else y - 1L
    graph.edgelist(y, directed = directed, ...)
}
