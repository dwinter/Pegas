##
# diversity.indices 
#
# calculate a few molcular diversity statistics for an aligment
##

diversity.indices <- function(d){
  haplo <- haplotype(d)
  result <- data.frame(
    "n" = dim(d)[1], 
    "n haplotypes" = dim(haplo)[1],
    "h" = h.div(haplo),
    "seggregating sites" = length(seg.sites(d)),
    "nucleotide diversity" =  nuc.div(d)
  )
  row.names(result) <- as.character(bquote(d))
  return(result)
}

#
#data(woodmouse)
#> diversity.indices(woodmouse)
#   n n.haplotypes h seggregating.sites nucleotide.diversity
#d 15           15 1                 56           0.01294610
#set.seed(123)
#> to.sample <- sample(1:dim(woodmouse)[1], 20, replace=TRUE)
#> d <- woodmouse[to.sample, ]
#> diversity.indices(d)
#   n n.haplotypes         h seggregating.sites nucleotide.diversity
#d 20           11 0.9368421                 55           0.01372922

