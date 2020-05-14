#' Extract table of differential results from 'K2' object
#'
#' Create table differential analysis results from 'K2' object.
#' @param K2res An object of class K2 or K2results().
#' @return A data.frame object.
#' @keywords clustering
#' @export
#' @examples
#' getDGETable(K2res)

getDGETable <- function(K2res) {
    
    ## If class is 'K2' extract K2results()
    if (class(K2res) == "K2") {
        K2res <- K2results(K2res)
    }
    
    ## Format single-sample enrichment results
    dgeTable <- do.call(rbind, lapply(names(K2res), function(x) {
        
        ## Get SSGSEA tables
        CompTab <- K2res[[x]]$dge
        
        ## Add group information and extract tables
        CompTab$split <- x
        
        ## Add direction information
        CompTab$direction <- c("down", "up")[as.numeric(CompTab$t > 0) + 1]
        
        return(CompTab)
    }))
    
    ## Sort by p-value
    dgeTable <- dgeTable[order(dgeTable$pval), ]
    
    rownames(dgeTable) <- NULL
    
    return(dgeTable)
}

