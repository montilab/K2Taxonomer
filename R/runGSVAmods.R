#' Perform ssGSVA on gene sets and get set signature hits
#'
#' Adds hyperenrichment analysis results to the output of runDGEmods().
#' @param K2res An object of class K2. The output of runDGEmods().
#' @param genesets A named list of feature IDs
#' @param qthresh A numeric value between 0 and 1 of the FDR cuttoff to define
#' feature sets.
#' @param cthresh A positive value for the coefficient cuttoff to define
#' feature sets.
#' @param ... Additional arguments passed onto GSVA::gsva()
#' @return An object of class K2.
#' @keywords clustering
#' @export
#' @import GSVA
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
#'                stabThresh = 0.5)
#' 
#' ## Run differential analysis on each partition
#' K2res <- runDGEmods(K2res)
#' 
#' ## Create dummy set of gene sets
#' DGEtable <- getDGETable(K2res)
#' genes <- unique(DGEtable$gene)
#' genesetsMadeUp <- list(
#'     GS1 = genes[1:50],
#'     GS2 = genes[51:100],
#'     GS3 = genes[101:150]
#' )
#'
#' ## Run gene set hyperenrichment
#' K2res <- runGSEmods(K2res, 
#'                     genesets = genesetsMadeUp,
#'                     qthresh = 0.1)
#'                     
#' ## Run GSVA on genesets           
#' K2res <- runGSVAmods(K2res, 
#'                      ssGSEAalg = "gsva",
#'                      ssGSEAcores = 1,
#'                      verbose = FALSE)
#'

runGSVAmods <- function(K2res, ssGSEAalg = NULL, ssGSEAcores = NULL, ...) {
    
    ## Run checks
    .isK2(K2res)
    
    ## Change meta data if new value is specific
    K2meta(K2res)$ssGSEAalg <- .checkK2(K2res, "ssGSEAalg", ssGSEAalg)
    K2meta(K2res)$ssGSEAcores <- .checkK2(K2res, "ssGSEAcores", ssGSEAcores)
    
    ## Check K2 object
    k2Check <- .checkK2(K2res)

    ## Run GSVA
    K2gSet(K2res) <- gsva(K2eSet(K2res), method = K2meta(K2res)$ssGSEAalg, gset.idx.list = K2genesets(K2res),
        parallel.sz = K2meta(K2res)$ssGSEAcores, ...)

    return(K2res)
}
