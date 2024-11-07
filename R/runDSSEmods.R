#' Perform differential analysis of single-sample gene set enrichment
#'
#' Adds limma differential analysis results of single-sample enrichment scores
#' to to the output of K2tax().
#' @param K2res An object of class K2. The output of K2tax().
#' @return An object of class K2.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{limma}{K2Taxonomer}
#'  \insertRef{bh}{K2Taxonomer}
#'  \insertRef{gsva}{K2Taxonomer}
#' @keywords clustering
#' @export
#' @import limma
#' @import Biobase
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
#'             stabThresh=0.5)
#'
#' ## Run differential analysis on each partition
#' K2res <- runDGEmods(K2res)
#'
#' ## Create dummy set of gene sets
#' DGEtable <- getDGETable(K2res)
#' genes <- unique(DGEtable$gene)
#' genesetsMadeUp <- list(
#'     GS1=genes[1:50],
#'     GS2=genes[51:100],
#'     GS3=genes[101:150])
#'
#' ## Run gene set hyperenrichment
#' K2res <- runGSEmods(K2res,
#'                 genesets=genesetsMadeUp,
#'                 qthresh=0.1)
#'
#' ## Run GSVA on genesets
#' K2res <- runGSVAmods(K2res,
#'                 ssGSEAalg='gsva',
#'                 ssGSEAcores=1,
#'                 verbose=FALSE)
#'
#' ## Run differential analysis on GSVA results
#' K2res <- runDSSEmods(K2res)
#'

runDSSEmods <- function(K2res,
                        cohorts = NULL,
                        vehicle = NULL,
                        covariates = NULL,
                        block = NULL) {

    ## Run checks
    .isK2(K2res)

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
