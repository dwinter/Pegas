#
#Calculate summary statistics for one or more populations from a single
#alignment
#

compare.diversity <- function(d, sites){
  #TODO 
  #catch error re passing a haplotype (not DNAbin)
  n <-table(sites)
  names(sites) <- 1:dim(d)[1]
  h <- haplotype(d)
  u.sites <- unique(sites)
  h.index <- sapply(attr(h, 'index'), function(s) unique(sites[s]))
  private <- sapply(u.sites, function(s) sum(h.index==s))
  n.haplo <- sapply(u.sites, function(s) sum(sapply(h.index, function(h) s %in% h)))
  pi <- sapply(u.sites, function(s) nuc.div(d[sites==s,]))
  segregating.sites <- sapply(unique(sites), function(s) length(seg.sites(d[sites==s,])))
  h.diversity <- sapply(sapply(u.sites, function(s) haplotype(d[sites==s,])), h.div)
  return(cbind(n, n.haplo,private,h.diversity,pi, segregating.sites))
}

#> compare.diversity(woodmouse, c(rep('ni', 10), rep('si', 5)))
#    n n.haplo private h.diversity         pi segregating.sites
#ni 10      10      10           1 0.01222921                45
#si  5       5       5           1 0.01540984                30
 

