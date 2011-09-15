##
# diversity.indices 
#
# summarise a few molcular diversity statistics for an aligment
##

diversity.indices <- function(d){
  haplo <- haplotype(d)
  result <- c(
    "n" = dim(d)[1], 
    "n haplotypes" = dim(haplo)[1],
    "h" = h.div(haplo),
    "seggregating sites" = length(seg.sites(d)),
    "nucleotide diversity" =  nuc.div(d)
  )
  #class(result) <- "diversity.index"
  return(result)
}



