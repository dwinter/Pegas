#
#Calculate the expected distribution of mismatches for an alignment under teh 
#assumption of stable population size
#

MMD.stable <- function(alignment){
  mismatch <- dist.dna(alignment, "N")
  theta = mean(mismatch)
  upper = ceiling(max(mismatch))
  count <- sapply(0:upper, function(i) (theta^i) / ((theta + 1)^(i+1)))
  return(count)
}

#> MMD(woodmouse, legend=FALSE)
#> lines(MMD.stable(woodmouse), col='red')
#> legend('topleft', c("Emperical", "Stable expectation"), col=c("blue", "red"), bg="white", lty=1)

