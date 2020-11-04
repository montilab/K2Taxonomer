#' Extract table of differential results from 'K2' object
#'
#' Create table differential analysis results from 'K2' object.
#' @param K2res An object of class K2.
#' @return A data.frame object.
#' @keywords clustering
#' @export
#' @examples
#' ## Read in ExpressionSet object
#' library(Biobase)
#' data(sample.ExpressionSet)
#' 
#' ## Pre-process and create K2 object
#' K2res <- K2preproc(sample.ExpressionSet)
#' 
#' ## Run K2 Taxonomer algorithm
#' K2res <- K2tax(K2res,
#'                stabThresh = 0.5)
#' 
#' ## Run differential analysis on each partition
#' K2res <- runDGEmods(K2res)
#' 
#' ## Extract table of differential results
#' head(getDGETable(K2res))
#' 

getDGETable <- function(K2res) {
    
    ## Run checks
    .isK2(K2res)
    
    K2resList <- K2results(K2res)
    
    ## Format single-sample enrichment results
    dgeTable <- do.call(rbind, lapply(names(K2resList), function(x, K2resList) {
        
        ## Get SSGSEA tables
        CompTab <- K2resList[[x]]$dge
        
        ## Add group information and extract tables
        CompTab$node <- x
        
        ## Add direction information
        CompTab$direction <- c("down", "up")[as.numeric(CompTab$t > 0) + 1]
        
        return(CompTab)
    }, K2resList))
    
    ## Sort by p-value
    dgeTable <- dgeTable[order(dgeTable$pval), ]
    
    rownames(dgeTable) <- NULL
    
    return(dgeTable)
}

