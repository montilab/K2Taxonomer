#' Extract gene set hits from K2 object
#'
#' Creates names list of genes found in gene sets for different results.
#' @param K2res An object of class K2. The output of runDGE_mods().
#' getDGEhits(K2res)

getDGEhits <- function(K2res) {
  
  # Get gene list subset for each hit set
  gseHitList <- unlist(lapply(names(K2results(K2res)), function(x) {
    
    K2node <- K2results(K2res)[[x]]
    gseList <- unlist(lapply(names(K2node$gse), function(y) {
      
      gseFram <- K2node$gse[[y]]
      gseFram <- gseFram[gseFram$fdr <= K2meta(K2res)$ssGSEAminQvalue,]
      if (nrow(gseFram) > 0) {
        hits <- strsplit(gseFram$hits, ",")
        names(hits) <- paste(gseFram$category, x, y, sep = "_")
        hits <- hits[unlist(lapply(hits, length)) > 0]
      }
    }), recursive = FALSE)
    
  }), recursive = FALSE)
  
  return(gseHitList)
}