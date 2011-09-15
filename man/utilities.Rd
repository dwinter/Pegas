\name{utilities}
\alias{getPloidy}
\alias{getAlleles}
\alias{getGenotypes}
\alias{expand.genotype}
\title{Utily Functions for pegas}
\description{
  The first three functions extract information on loci;
  \code{expand.genotype} creates a table of all possible genotypes given
  a set of alleles.
}
\usage{
getPloidy(x)
getAlleles(x)
getGenotypes(x)
expand.genotype(n, alleles = NULL, ploidy = 2, matrix = FALSE)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{n}{an integer giving how many alleles to consider (ignored if
    \code{alleles} is used).}
  \item{alleles}{the allele names as a vector of mode character.}
  \item{ploidy}{an integer giving the ploidy level (either 2 or 4 for
    the moment).}
  \item{matrix}{a logical specifying whether to return the genotypes in
    a matrix or as a character vector.}
}
\value{
  \code{getPloidy} returns the ploidy level of all loci in an object of
  class \code{"loci"} as a numeric vector.

  \code{getAlleles} and \code{getGenotpes} return the alleles and
  genotypes, respectively, observed in all loci in an object of class
  \code{"loci"} as a list.

  \code{expand.genotype} returns a character vector (the default) or a
  matrix where the rows are the genotypes and the columns are the
  alleles. The matrix is numeric by default, or character if the
  argument \code{alleles} is given.
}
\author{Emmanuel Paradis}
\examples{
require(adegenet)
data(nancycats)
X <- as.loci(nancycats)[, 2:3]
getAlleles(X)
getGenotypes(X)
expand.genotype(2)
expand.genotype(2, LETTERS[1:3])
expand.genotype(3, ploidy = 4)
}
\keyword{manip}
