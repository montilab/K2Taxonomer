#' Extract table of differential results from 'K2' object
#'
#' Create table differential analysis results from 'K2' object.
#' @param K2res An object of class K2.
#' @return A data.frame object with the following columns: 
#' \itemize{
#'  \item{gene: }{The gene ids from the `ExpressionSet` object}
#'  \item{coef: }{The difference between the means of each subgroup at 
#'  a given partition}
#'  \item{mean: }{The mean across all observations at the given 
#'  partition}
#'  \item{t: }{The test statistic estimated by the `limma` R package}
#'  \item{pval: }{The p-value estimated by the `limma` R package}
#'  \item{fdr: }{The multiple hypothesis corrected fdr p-value, adjusted across all partitions}
#'  \item{B: }{The B-statistic estimated by the `limma` R package}
#'  \item{edge: }{Indication of which subgroup the gene was assigned at a given partition}
#'  \item{node: }{The identifier of the partition}
#'  \item{direction: }{The direction of coefficient for the assigned gene}
#' }
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{limma}{K2Taxonomer}
#'  \insertRef{bh}{K2Taxonomer}
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
        
        ## Get tables
        CompTab <- K2resList[[x]]$dge
        
        if(!is.null(CompTab)) {
            ## Add group information and extract tables
            CompTab$node <- x
            
            ## Add direction information
            CompTab$direction <- c("down", "up")[as.numeric(CompTab$t > 0) + 1]
        }
            
        return(CompTab)
    }, K2resList))
    
    ## Sort by p-value
    dgeTable <- dgeTable[order(dgeTable$pval), ]
    
    rownames(dgeTable) <- NULL
    
    return(dgeTable)
}

