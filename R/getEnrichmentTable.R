#' Extract table of enrichment results from "K2" object
#'
#' Create table hyper- and single-sample enrichment results from "K2" object.
#' @param K2res An object of class K2 or K2results().
#' @return A data.frame object.
#' @keywords clustering
#' @export
#' @examples
#' getEnrichementTable(K2res)
  
getEnrichmentTable <- function(K2res) {
  
  # If class is "K2" extract K2results()
  if(class(K2res) == "K2") {
    K2res <- K2results(K2res)
  }
  
  # Format hyperenrichment table
  EnrTable <- do.call(rbind, lapply(names(K2res), function(x) {
    
    # Get GSE tables
    GSEtabList <- K2res[[x]]$gse
    
    # Add group informaiton and extract tables
    GSEtab <- do.call(rbind, lapply(names(GSEtabList), function(y) {
      GSEsub <- GSEtabList[[y]]
      if(nrow(GSEsub) > 0) {
        GSEsub$split <- x # Add split ID
        GSEsub$mod <- gsub("g|_up|_down", "", y) # Add group ID
        GSEsub$direction <- gsub("g1_|g2_", "", y) # Add  direction
      }
      return(GSEsub)
    }))
    
    # Make pval and fdr column names unique
    colnames(GSEtab)[colnames(GSEtab) %in% c("pval", "fdr")] <- 
      paste(colnames(GSEtab)[colnames(GSEtab) %in% c("pval", "fdr")], "hyper", sep = "_")
    
    return(GSEtab)
  }))
  
  if(!is.null(K2res[[1]]$dss)) {
    # Format single-sample enrichment results
    ssEnrTable <- do.call(rbind, lapply(names(K2res), function(x) {
      
      # Get SSGSEA tables
      SSGSEAtab <- K2res[[x]]$dsse
      
      # Add group information and extract tables
      SSGSEAtab$split <- x
      
      # Add direction information
      SSGSEAtab$direction <- c("down", "up")[as.numeric(SSGSEAtab$t > 0) + 1]
      
      # Make pval and fdr column names unique
      colnames(SSGSEAtab)[colnames(SSGSEAtab) %in% c("pval", "fdr")] <- 
        paste(colnames(SSGSEAtab)[colnames(SSGSEAtab) %in% c("pval", "fdr")], "limma", sep = "_")
      
      return(SSGSEAtab)
    }))
    
    # Merge the two and sort by hyper p-value
    EnrTable <- merge(EnrTable, ssEnrTable, all.x = TRUE)
    
    
    # Sort columns
    EnrTable <- EnrTable[, c("category", "split", "mod", "direction", "pval_hyper", "fdr_hyper", "nhits", "ndrawn", "ncats", "ntot", "pval_limma", "fdr_limma", "coef", "mean", "t", "B", "hits")]
    
  }
  
  # Sort by p-value
  EnrTable <- EnrTable[order(EnrTable$pval_hyper),]
  
  rownames(EnrTable) <- NULL
  
  return(EnrTable)
}

