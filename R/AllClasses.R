#' K2 object
#'
#' This class is the output of runK2Taxonomer() and K2preproc().
#' @slot eSet An Expression Set object
#' @slot meta A named list of K2Taxonomer information
#' @slot dataMatrix An P x N numeric matrix of data
#' @slot results A named list of K2Taxonomer partitioning reuslts
#' @slot info A data frame with rownames that match column names
#' in dataMatrix
#' @slot genesets A named list of features in row names of dataMatrix
#' @slot gene2Pathway A vector of collapsed genesets names, mapping
#' features to genesets
#' @slot gSet An Expression Set object
#' @slot geneURL A named vector of gene URLs
#' @slot genesetURL A named vector of geneset URLs
#'
#' @keywords clustering
#' @export

setClass("K2", representation(eSet="ExpressionSet", meta="list",
    dataMatrix="matrix", info="data.frame", results="list",
    genesets="list", gene2Pathway="character", gSet="ExpressionSet",
    geneURL="character", genesetURL="character"))

setMethod("show", "K2", function(object) {
    cat("K2Taxonomer results:", "\n",
        " Expression Data: ", ifelse(
            is.null(object@eSet), "FALSE", "TRUE"), "\n",
        " Pre-processing: ", ifelse(
            is.null(object@meta) &
                is.null(object@meta) &
                is.null(object@dataMatrix) &
                is.null(object@info),
            "FALSE", "TRUE"), "\n",
        " Recursive partitioning: ", ifelse(
            is.null(object@results), "FALSE", "TRUE"), 
        "\n",
        " Differential Gene Expression: ", ifelse(
            is.null(object@results$A$dge), "FALSE", "TRUE"), 
        "\n",
        " Gene Set Hyperenrichment: ", ifelse(
            is.null(object@results$A$gse), "FALSE", "TRUE"), 
        "\n",
        " Single-sample Enrichment: ", ifelse(
            is.null(object@gSet), "FALSE", "TRUE"), 
        "\n",
        " Differential Enrichment Analysis: ", ifelse(
            is.null(object@results$A$dsse), "FALSE", "TRUE"), 
        "\n",
        " Phenotypic Variable Testing: ", ifelse(
            is.null(object@results$A$modTests), "FALSE", "TRUE"), 
        "\n", sep = ""
    )
})
