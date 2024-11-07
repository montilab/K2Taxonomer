.K2filterGenes <- function(e, m, thr) {
  eg0 <- e > 0
  pu <- unique(m)
  mk <- rowMax(do.call(cbind, lapply(pu, function(x) {
    es <- eg0[,m == x]
    rowMeans(es)
  }))) > thr
  return(e[mk,])
}