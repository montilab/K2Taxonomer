#' Perform differential analysis of enrichment scores between subgroups at each
#' partition
#'
#' Adds differential analysis results of single-sample enrichment scores
#' to to the output of K2tax().
#' @return An object of class K2.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{limma}{K2Taxonomer}
#'  \insertRef{bh}{K2Taxonomer}
#'  \insertRef{gsva}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2preproc
#' @export
#' @import limma
#' @import Biobase

runDSSEmods <- function(K2res,
                        vehicle = NULL,
                        variables = NULL,
                        block = NULL
                        ) {

    ## Run checks
    .isK2(K2res)
    
    K2meta(K2res)$vehicle <- vehicle <- .checkK2(K2res, "vehicle",
                                                 vehicle)
    K2meta(K2res)$variables <- variables <- .checkK2(K2res, "variables",
                                                     variables)
    K2meta(K2res)$block <- block <- .checkK2(K2res, "block",
                                             block)

    ## Check K2 object
    k2Check <- .checkK2(K2res)

    ## K2 algorithm
    if (length(K2results(K2res)) == 0) {
        stop("No results found. Please run K2tax() or runK2Taxonomer().\n")
    }

    ## DGE
    if (is.null(K2results(K2res)[[1]]$dge)) {
        stop("No differential analysis results found. Please run runDGEmods().\n")
    }

    ## GSE
    if (is.null(K2results(K2res)[[1]]$gse)) {
        stop("No enrichment results found. Please run runDGEmods().\n")
    }

    ## GSVA
    if (ncol(K2gMat(K2res)) == 0) {
        stop("No ssGSEA data found. Please run runGSVAmods().\n")
    }
    
    modVec <- unlist(lapply(K2results(K2res), function(x) paste(x$obs[[1]], collapse = "_")))
    modL <- length(modVec)

    
    cat("Running different enrichment score for partition:\n")
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
        dsseRes <- .signatureWrapper(K2res, mods, GENE = FALSE)
        
        x$dsse <- dsseRes$modStats
        x$dsseFormula <- dsseRes$formula
        
        if (!is.null(x$dsse)) {
          colnames(x$dsse)[1] <- "category"
        }

        return(x)
    })

    ## Fix FDR values
    K2res <- .fixFDR(K2res, "dsse")

    return(K2res)
}
