stable.dist <- function(alignment){
  mismatch <- dist.dna(alignment, "N")
  theta = mean(mismatch)
  upper = ceiling(max(mismatch))
  mm <- 0:upper
  count <- numeric()
  for (i in 0:upper){
    count[i] <- (theta^i) / ((theta + 1)^(i+1) )
  }
  return(count)
}
