#' K2 object
#'
#' This class is the output of runK2Taxonomer().
#' @slot dataMatrix An P x N numeric matrix of data
#' @slot info A data frame with rownames that match column names in dataMatrix
#' @slot genesets A named list of features in row names of dataMatrix
#' @slot gene2Pathway A vector of collapsed genesets names, mapping features to genesets
#' @slot eSet An Expression Set object
#' @slot gSet An Expression Set object
#' @slot meta A named list of K2 Taxonomer information
#' @slot geneURL A named vector of gene URLs
#' @slot genesetURL A named vector of geneset URLs
#' 
#' @keywords clustering
#' @export

setClass("K2", representation(eSet = "ExpressionSet",      # Expression Set
                              meta = "list",               # Named list of K2 Taxonomer parameters
                              dataMatrix = "matrix",       # Matrix of Z-scored data
                              info = "data.frame",         # Data frame of observation information
                              results = "list",            # K2 Taxonomer results 
                              genesets = "list",           # Named list of gene sets for hyperenrichment analysis
                              gene2Pathway = "character",  # Named list mapping genes to gene sets
                              gSet = "ExpressionSet",      # Expression Set (from GSVA)
                              geneURL = "character",       # Names vector of gene URLs
                              genesetURL = "character"    # Namec vector of geneset URLs
                              ))