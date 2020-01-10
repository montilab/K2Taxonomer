#' Take difference of two paired GSVA scores
#'
#' Replaces GSVA results from paired up- and down- gene sets with the difference of
#' the up-regulated genes and down-regulated genes
#' @param aggList A list where each item is a vector of 3 items:
#' the new name, the name of the 'up' gene set, and the name of the 'down'
#' gene set.
#' @param K2res An object of class K2. The output of runDGE_mods().
#' @return An object of class K2.
#' @keywords clustering
#' @export
#' @import GSVA
#' @import Biobase
#' @examples
#' aggregateGSVAscores(aggList, K2res)

aggregateGSVAscores <- function(aggList, K2res) {
  
  # Only perform this method for gsva runs
  if(K2meta(K2res)$ssGSEAalg != "gsva") {
    stop("Enrichment score aggregation only supported when ssGSEAalg == 'gsva'")
  }
  
  # Get gSet
  gSet <- K2gSet(K2res)
  
  # For each item in aggList subtract negative signature score from positive
  gAgg <- do.call(rbind, lapply(aggList, function(agg) {
    gSubUp <- exprs(gSet)[agg[2],]
    gSubDown <- exprs(gSet)[agg[3],]
    gSubAgg <- gSubUp - gSubDown
  }))
  rownames(gAgg) <- unlist(lapply(aggList, function(x) x[1]))
  
  # Remove original scores
  origNames <- unlist(lapply(aggList, function(x) x[-1]))
  gSet <- gSet[ !rownames(gSet) %in% origNames,]
  
  # Add new scores
  gExprs <- rbind(exprs(gSet), gAgg)
  gNew <- ExpressionSet(assayData = gExprs)
  pData(gNew) <- pData(gSet)
  
  # Replace gSet
  K2gSet(K2res) <- gNew
  
  return(K2res)
  
}