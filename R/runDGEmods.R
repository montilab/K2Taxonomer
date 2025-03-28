#' Perform differential analysis between subgroups at each partition
#'
#' Adds differential analysis results to the output of K2tax().
#' @return An object of class K2.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{limma}{K2Taxonomer}
#'  \insertRef{bh}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2preproc
#' @export
#' @import limma
#' @import Biobase

runDGEmods <- function(K2res,
                       DGEmethod = NULL,
                       DGEexpThreshold = NULL,
                       vehicle = NULL,
                       variables = NULL,
                       block = NULL,
                       logCounts = NULL) {

    ## Run checks
    .isK2(K2res)
  
    ## Change meta data if new value is specific
    K2meta(K2res)$DGEmethod <- DGEmethod <- .checkK2(K2res, "DGEmethod",
                                                     DGEmethod)
    K2meta(K2res)$DGEexpThreshold <- DGEexpThreshold <- 
      .checkK2(K2res, "DGEexpThreshold", DGEexpThreshold)
    K2meta(K2res)$vehicle <- vehicle <- .checkK2(K2res, "vehicle",
                                                 vehicle)
    K2meta(K2res)$variables <- variables <- .checkK2(K2res, "variables",
                                                       variables)
    K2meta(K2res)$block <- block <- .checkK2(K2res, "block",
                                             block)
    K2meta(K2res)$logCounts <- logCounts <- .checkK2(K2res, "logCounts",
                                                     logCounts)
  

    ## Check K2 object
    k2Check <- .checkK2(K2res)

    ## K2 algorithm
    if (length(K2results(K2res)) == 0) {
        stop("No results found. Please run K2tax().\n")
    }
    
    modVec <- unlist(lapply(K2results(K2res), function(x) paste(x$obs[[1]], collapse = "_")))
    modL <- length(modVec)
    
    cat("Running DGE for partition:\n")
    K2results(K2res) <- lapply(K2results(K2res), function(x) {
      
        ## Create module variable
        obs1 <- x$obs[[1]]
        obs2 <- x$obs[[2]]
      
        mods <- as.factor(c(rep(1, length(obs1)), rep(2,
            length(obs2))))
        names(mods) <- c(x$obs[[1]], x$obs[[2]])
        
        ## Print progress
        cat(" ", which(modVec == paste(obs1, collapse = "_")), "/", modL, "\n")

        ## Perform differential analysis
        if(DGEmethod == "limma") {
          dgeRes <- .signatureWrapper(K2res, mods, GENE = TRUE)
        } else {
          dgeRes <- .signatureWrapper_MAST(K2res, mods, GENE = TRUE)
        }
        
        x$dge <- dgeRes$modStats
        x$dgeFormula <- dgeRes$formula
        
        if (!is.null(x$dge)) {
          colnames(x$dge)[1] <- "gene"
        }

        ## Set x$gse to NULL, if values are here they need to be
        ## re-run
        x$gse <- NULL
        x$dsse <- NULL

        return(x)
    })

    ## Fix FDR values
    K2res <- .fixFDR(K2res, "dge")

    ## Make K2gSet empty
    K2gMat(K2res) <- new("matrix")

    return(K2res)
}
