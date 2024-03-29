\name{as.loci}
\alias{as.loci}
\alias{as.loci.genind}
\alias{as.loci.data.frame}
\alias{as.loci.factor}
\alias{as.loci.character}
\alias{genind2loci}
\alias{loci2genind}
\title{Conversion Among Allelic Data Classes}
\description{
  These functions do conversion among different allelic data classes.
}
\usage{
as.loci(x, ...)
\method{as.loci}{genind}(x, ...)
genind2loci(x)
\method{as.loci}{data.frame}(x, allele.sep = "/", col.pop = NULL, col.loci = NULL, ...)
loci2genind(x)
\method{as.loci}{factor}(x, allele.sep = "/", ...)
\method{as.loci}{character}(x, allele.sep = "/", ...)
}
\arguments{
  \item{x}{an object of class \code{"loci"} or \code{"genind"}, a data
    frame, a factor, or a vector of mode character.}
    \item{allele.sep}{the character(s) separating the alleles for each locus
    in the data file (a forward slash by default).}
  \item{col.pop}{specifies whether one of the column of the data file
    identifies the population; default \code{NULL}, otherwise an integer
    or a character giving the number or the name of the column.}
  \item{col.loci}{a vector of integers or of characters specifying the
    indices or the names of the columns that are loci. By default, all
    columns are taken as loci except the one labelled "population", if
    present or specified.}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
  The main objectives of these functions is to provide easy conversion
  between the data structures of \pkg{adegenet} and \pkg{pegas}, so both
  packages can be used together smoothly. In addition, it is possible to
  create a \code{"loci"} object directly from a data frame, a vector, or
  a factor.

  \code{genind2loci(x)} and \code{as.loci(x)} are the same if \code{x}
  is of class \code{"genind"}.
}
\value{
  An object of class \code{c("loci", "data.frame")} for \code{as.loci}
  and \code{genind2loci}; an object of class \code{"genind"} for
  \code{loci2genind}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.loci}}, \code{\link[adegenet]{genind}},
  \code{\link[adegenet]{df2genind}} for converting data frames to \code{"genind"}
}
\examples{
x <- c("A-A", "A-a", "a-a")
as.loci(x, allele.sep = "-")
require(adegenet)
data(nancycats)
x <- as.loci(nancycats)
y <- loci2genind(x) # back to "genind"
identical(nancycats@tab, y@tab)
identical(nancycats@pop, y@pop)
}
\keyword{IO}
