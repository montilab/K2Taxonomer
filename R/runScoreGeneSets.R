#' Perform ssGSVA on gene sets and get set signature hits
#'
#' Adds hyperenrichment analysis results to the output of runDGEmods().
#' @param K2res An object of class K2. The output of runDGEmods().
#' @param useCors Number of cores to use for running gsva() from the GSVA
#' package. Default is 1.
#' @param ... Additional arguments passed onto GSVA::gsva()
#' @return An object of class K2.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{gsva}{K2Taxonomer}
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
#'                 useCors=1,
#'                 verbose=FALSE)
#'

runScoreGeneSets <- function(K2res, ScoreGeneSetMethod = NULL, useCors=NULL) {

    ## Run checks
    .isK2(K2res)

    ## Change meta data if new value is specific
    K2meta(K2res)$useCors <- .checkK2(K2res, "useCors",
        useCors)
    K2meta(K2res)$ScoreGeneSetMethod <- .checkK2(K2res, "ScoreGeneSetMethod",
        ScoreGeneSetMethod)

    ## Check K2 object
    k2Check <- .checkK2(K2res)
    
    # Get function for expression data matrix
    if(nrow(K2eMatDS(K2res)) != 0) {
      EXPFUNC <- K2eMatDS
    } else {
      EXPFUNC <- K2eMat
    }
    
    if(K2meta(K2res)$useCors > 1) {
      bF <- get(class(bpparam())[[1]])
      bpp <- bF(workers = K2meta(K2res)$useCors)
    } else {
      bpp <- SerialParam()
    }
    
    if(K2meta(K2res)$ScoreGeneSetMethod == "GSVA") {
      cat("Scoring gene sets with GSVA.\n")
      gP <- gsvaParam(EXPFUNC(K2res), K2genesets(K2res))
      K2gMat(K2res) <- suppressWarnings(gsva(gP, BPPARAM = bpp))
    }
    
    if(K2meta(K2res)$ScoreGeneSetMethod == "AUCELL") {
      cat("Scoring gene sets with AUCell.\n")
      K2gMat(K2res) <- log2(AUCell_run(EXPFUNC(K2res), K2genesets(K2res), BPPARAM = bpp)@assays@data$AUC * 100 + 1)
    }

    return(K2res)
}
