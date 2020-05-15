#' Perform differential analysis of single-sample gene set enrichment
#'
#' Adds limma differential analysis results of single-sample enrichment scores
#' to to the output of K2tax().
#' @param K2res An object of class K2. The output of K2tax().
#' @return An object of class K2.
#' @keywords clustering
#' @export
#' @import limma
#' @import Biobase
#' @examples
#' runDSSEmods(K2res)

runDSSEmods <- function(K2res) {
    
    ## Run checks
    .isK2(K2res)
    
    ## Check K2 object
    k2Check <- .checkK2(K2res)
    
    ## K2 algorithm
    if (length(K2results(K2res)) == 0) {
        "No results found. Please run K2tax() or runK2Taxonomer().\n"
    }
    
    ## DGE
    if (is.null(K2results(K2res)[[1]]$dge)) {
        "No differential analysis results found. Please run runDGEmods().\n"
    }
    
    ## GSE
    if (is.null(K2results(K2res)[[1]]$gse)) {
        "No enrichment results found. Please run runDGEmods().\n"
    }
    
    ## GSVA
    if (ncol(K2gSet(K2res)) == 0) {
        "No ssGSEA data found. Please run runGSVAmods().\n"
    }
    
    K2results(K2res) <- lapply(K2results(K2res), function(x) {
        
        ## Create module variable
        mods <- as.factor(c(rep(1, length(x$obs[[1]])), rep(2, length(x$obs[[2]]))))
        names(mods) <- c(x$obs[[1]], x$obs[[2]])
        
        ## Perform differential analysis
        x$dsse <- .signatureWrapper(K2gSet(K2res), K2meta(K2res)$cohorts, mods, K2meta(K2res)$vehicle, 
            K2meta(K2res)$covariates, K2meta(K2res)$block)
        if (!is.null(x$dsse)) {
            x$dsse$category <- rownames(x$dsse)
            x$dsse <- x$dsse[, c(ncol(x$dsse), 1:(ncol(x$dsse) - 1))]
        }
        
        return(x)
    })
    
    ## Fix FDR values
    K2res <- .fixFDR(K2res, "dsse")
    
    return(K2res)
}
