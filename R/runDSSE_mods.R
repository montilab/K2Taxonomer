#' Perform differential analysis of single-sample gene set enrichment
#'
#' Adds limma differential analysis results of single-sample enrichment scores to to the output of K2tax().
#' @param K2res An object of class K2. The output of K2tax().
#' @return An object of class K2.
#' @keywords clustering
#' @export
#' @import limma
#' @import Biobase
#' @examples
#' runDSSE_mods(K2res)

runDSSE_mods <- function(K2res){
  
  K2results(K2res) <- lapply(K2results(K2res), function(x){

    # Create module variable
    mods <- as.factor(c(rep(1, length(x$obs[[1]])), rep(2, length(x$obs[[2]]))))
    names(mods) <- c(x$obs[[1]], x$obs[[2]])
    
    # Perform differential analysis
    x$dsse <- .signature_wrapper(K2gSet(K2res), 
                                K2meta(K2res)$cohorts, 
                                mods, 
                                K2meta(K2res)$vehicle, 
                                K2meta(K2res)$covariates,
                                K2meta(K2res)$block)
    if (!is.null(x$dsse)) {
      x$dsse$category <- rownames(x$dsse)
      x$dsse <- x$dsse[, c(ncol(x$dsse), 1:(ncol(x$dsse)-1)) ]
    }
    
    return(x)
  })
  
  # Fix FDR values
  K2res <- .fixFDR(K2res, "dsse")
  
  return(K2res)
}
