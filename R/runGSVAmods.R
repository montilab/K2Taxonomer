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
#' runGSVAmods(K2res)

runGSVAmods <- function(K2res, ssGSEAalg = NULL, ssGSEAcores = NULL, ...) {

    ## Change meta data if new value is specific
    K2meta(K2res)$ssGSEAalg <- .checkMeta(K2res, "ssGSEAalg", ssGSEAalg)
    K2meta(K2res)$ssGSEAcores <- .checkMeta(K2res, "ssGSEAcores", ssGSEAcores)

    ## Run GSVA
    K2gSet(K2res) <- gsva(K2eSet(K2res), method = K2meta(K2res)$ssGSEAalg, gset.idx.list = K2genesets(K2res),
        parallel.sz = K2meta(K2res)$ssGSEAcores, ...)

    return(K2res)
}
