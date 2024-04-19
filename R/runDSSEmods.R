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

runDSSEmods <- function(K2res) {

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
    if (ncol(K2gSet(K2res)) == 0) {
        stop("No ssGSEA data found. Please run runGSVAmods().\n")
    }

    K2results(K2res) <- lapply(K2results(K2res), function(x) {

        ## Create module variable
        mods <- as.factor(c(rep(1, length(x$obs[[1]])), rep(2,
            length(x$obs[[2]]))))
        names(mods) <- c(x$obs[[1]], x$obs[[2]])

        ## Perform differential analysis
        dsseRes <- .signatureWrapper(K2gSet(K2res), K2meta(K2res)$cohorts,
            mods, K2meta(K2res)$vehicle, K2meta(K2res)$covariates,
            K2meta(K2res)$block)
        x$dsse <- dsseRes$modStats
        x$dsseFormula <- dsseRes$formula
        if (!is.null(x$dsse)) {
            x$dsse$category <- rownames(x$dsse)
            x$dsse <- x$dsse[, c(ncol(x$dsse), seq_len(ncol(x$dsse) -
                1))]
        }

        return(x)
    })

    ## Fix FDR values
    K2res <- .fixFDR(K2res, "dsse")

    return(K2res)
}
