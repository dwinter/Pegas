\name{theta.tree}
\Rdversion{1.1}
\alias{theta.tree}
\title{
  Population Parameter THETA Using Genealogy
}
\description{
  This function estimates the population parameter \eqn{\theta}{THETA}
  from a genealogy (coded a as phylogenetic tree) under the coalescent.
}
\usage{
theta.tree(phy, theta, fixed = FALSE, log = TRUE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{theta}{a numeric vector.}
  \item{fixed}{a logical specifying whether to estimate \code{theta}
    (the default), or to return the likelihoods for all values in
    \code{theta}.}
  \item{log}{a logical specifying whether to return the likelihoods on a
    log scale (the default); ignored if \code{fixed = FALSE}.}
}
\details{
  The tree \code{phy} is considered as a genealogy, and therefore should
  be ultrametric. By default, \eqn{\theta}{THETA} is estimated by
  maximum likelihood and the value given in \code{theta} is used as
  starting value for the minimisation function (if several values are
  given as a vector the first one is used). If \code{fixed = TRUE},
  then the [log-]likelihood values are returned corresponding to each
  value in \code{theta}.

  The present implementation does a numerical optimisation of the
  log-likelihood function (with \code{\link[stats]{nlminb}}) with the
  first partial derivative as gradient. It is possible to solve the
  latter and have a direct analytical MLE of \eqn{\theta}{THETA} (and
  its standard-error), but this does not seem to be faster.
}
\value{
  If \code{fixed = FALSE}, a list with two elements:

  \item{theta}{the maximum likelihood estimate of \eqn{\theta}{THETA};}
  \item{logLik}{the log-likelihood at its maximum.}

  If \code{fixed = TRUE}, a numeric vector with the [log-]likelihood
  values.
}
\references{
  Kingman, J. F. C. (1982) The coalescent. \emph{Stochastic Processes
  and their Applications}, \bold{13}, 235--248.

  Kingman, J. F. C. (1982) On the genealogy of large
  populations. \emph{Journal of Applied Probability}, \bold{19A},
  27--43.

  Wakeley, J. (2009) Coalescent Theory: An Introduction. Greenwood
  Village, CO: Roberts and Company Publishers.
}
\author{
  Emmanuel Paradis
}
\seealso{
  \code{\link{theta.h}}, \code{\link{theta.s}}, \code{\link{theta.k}}
}
\examples{
tr <- rcoal(50) # assumes theta = 1
theta.tree(tr, 10)
## profile log-likelihood:
THETA <- seq(0.5, 1.5, 0.01)
logLikelihood <- theta.tree(tr, THETA, fixed = TRUE)
plot(THETA, logLikelihood, type = "l")
}
\keyword{models}
