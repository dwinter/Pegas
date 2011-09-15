\name{haplotype}
\alias{haplotype}
\alias{plot.haplotype}
\alias{print.haplotype}
\title{Haplotype Extraction and Frequencies}
\description{
  \code{haplotype} extracts the haplotypes from a set of DNA
  sequences. The result can be plotted with the appropriate function.
}
\usage{
haplotype(x, labels = NULL)
\method{plot}{haplotype}(x, ...)
\method{print}{haplotype}(x, ...)
}
\arguments{
  \item{x}{a set of DNA sequences (as an object of class
    \code{"DNAbin"}, or an object of class \code{"haplotype"}.}
  \item{labels}{a vector of character strings used as names for the rows
    of the returned object. By default, Roman numerals are given.}
  \item{\dots}{further arguments passed to
    \code{\link[graphics]{barplot}} (unused in \code{print}).}
}
\value{
  \code{haplotype} returns an object of class \code{c("haplotype",
    "DNAbin")} which is an object of class \code{"DNAbin"} with two
  additional attributes: \code{"index"} identifying the index of each
  observation that share the same haplotype, and \code{"from"} giving
  the name of the original data.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}}, \code{\link[ape]{DNAbin}} for manipulation of
  DNA sequences in R.
}
\examples{
## generate some artificial data from 'woodmouse':
data(woodmouse)
x <- woodmouse[sample(15, size = 110, replace = TRUE), ]
(h <- haplotype(x))
## the indices of the individuals belonging to the 1st haplotype:
attr(h, "index")[[1]]
}
\keyword{manip}
\keyword{hplot}
