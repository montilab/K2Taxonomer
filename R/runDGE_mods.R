#' Perform differential analysis of gene splits
#'
#' Adds limma differential analysis results to the output of K2tax().
#' @param K2res An object of class K2. The output of K2tax().
#' @return An object of class K2.
#' @keywords clustering
#' @export
#' @import limma
#' @import Biobase
#' @examples
#' runDGE_mods(K2res)

runDGE_mods <- function(K2res){
  
  K2results(K2res) <- lapply(K2results(K2res), function(x){
    
    # Create module variable
    mods <- as.factor(c(rep(1, length(x$obs[[1]])), rep(2, length(x$obs[[2]]))))
    names(mods) <- c(x$obs[[1]], x$obs[[2]])
    
    # Perform differential analysis
    x$dge <- .signature_wrapper(K2eSet(K2res), 
                                K2meta(K2res)$cohorts, 
                                mods, 
                                K2meta(K2res)$vehicle, 
                                K2meta(K2res)$covariates)
    if (!is.null(x$dge)) {
      x$dge$gene <- rownames(x$dge)
      x$dge <- x$dge[,c(ncol(x$dge), 1:(ncol(x$dge)-1))]
    }
    
    # Set x$gse to NULL, if values are here they need to be re-run
    x$gse <- NULL
    
    return(x)
  })
  
  # Make K2gSet empty
  K2gSet(K2res) <- ExpressionSet()
  
  return(K2res)
}