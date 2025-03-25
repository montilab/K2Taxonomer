#' K2 object
#'
#' This class is the output of K2preproc().
#' @slot eMat A numeric matrix comprising expression data used for partitioning
#' @slot eMatDS A numeric matrix comprising expression data used for differential expression
#' @slot gMat A numeric matrix comprising the output of runGeneSetScores().
#' @slot results A named list object containing results from each partition
#' @slot colData A data frame containing observation level data
#' @slot meta A named list containing different parameters used in workflow
#' @slot genesets A named list containing gene sets
#' @slot gene2Pathway Vector of collapsed gene set for which each gene belongs
#' @slot dataMatrix A numeric matrix used for different partitioning tasks
#' @slot geneURL A named list containing URLs for each gene
#' @slot genesetURL A named list containing URLs for each gene set
#' @keywords clustering
#' @export

setClassUnion("matrixORdgCMatrix", c("matrix", "dgCMatrix"))

setClass("K2", slots = c(meta="list",
                              eMat="matrixORdgCMatrix",
                              eMatDS ="matrixORdgCMatrix",
                              gMat="matrixORdgCMatrix",
                              colData = "data.frame",
                              dataMatrix="matrixORdgCMatrix",
                              results="list",
                              genesets="list", 
                              gene2Pathway="character", 
                              geneURL="character", 
                              genesetURL="character"))


setMethod("show", "K2", function(object) {
    cat("K2Taxonomer results:", "\n",
        " Expression Data: ", ifelse(
            ncol(object@eMat) == 0, "FALSE", "TRUE"), "\n",
        " Pre-processing: ", ifelse(
                length(object@meta) == 0 &
                ncol(object@dataMatrix) == 0 &
                ncol(object@colData) == 0,
            "FALSE", "TRUE"), "\n",
        " Recursive partitioning: ", ifelse(
            length(object@results) == 0, "FALSE", "TRUE"), 
        "\n",
        " Differential Gene Expression: ", ifelse(
            is.null(object@results$A$dge), "FALSE", "TRUE"), 
        "\n",
        " Gene Set Hyperenrichment: ", ifelse(
            is.null(object@results$A$gse), "FALSE", "TRUE"), 
        "\n",
        " Single-sample Enrichment: ", ifelse(
            ncol(object@gMat) == 0, "FALSE", "TRUE"), 
        "\n",
        " Differential Enrichment Analysis: ", ifelse(
            is.null(object@results$A$dsse), "FALSE", "TRUE"), 
        "\n", sep = ""
    )
})
