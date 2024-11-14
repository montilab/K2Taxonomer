#' Enrichment scoring of gene sets on expression data
#'
#' Performs single-sample enrichment scoring using GSVA or AUCell. When AUCell
#' is specified, output values reflect the AUC levels scaled by 100 and log2
#' transformed.
#' @return An object of class K2.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{gsva}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2preproc
#' @export
#' @import GSVA
#' @import Biobase
#' @import AUCell

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
