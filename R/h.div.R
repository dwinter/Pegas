#
# Calculate haplotype diversity (h) from a haplotype object.
#

h.div <- function(h){
  if (!inherits(h, "haplotype")) 
    stop("'h' must be of class 'haplotype'")
  freqs <- sapply(attr(h, "index"), length)
  n <- sum(freqs)
  little_h <- (n/(n-1)) * (1 - sum( (freqs/n)**2 ))
  return(little_h)
}
