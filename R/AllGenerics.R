#' @title Data matrix from K2 object
#' @description Retrieve or assign matrix of data used in K2Taxonomer run.
#' @param K2res K2 class object.
#' @return Matrix. Numeric data used for K2Taxonomer run.
#' @examples
#' K2data(K2res)
#' @export
setGeneric("K2data", function(K2res) standardGeneric("K2data"))

#' @export K2data<-
setGeneric("K2data<-", function(K2res, value) standardGeneric("K2data<-"))
setMethod("K2data", "K2", function(K2res) K2res@dataMatrix)
setMethod("K2data<-", "K2", function(K2res, value) {
    K2res@dataMatrix <- value
    K2res
})

#' @title Data frame of phenotypic information from K2 object
#' @description Retrieve or assign observation information in K2Taxonomer run.
#' @param K2res K2 class object.
#' @return Data frame. Observation information from K2Taxonomer run.
#' @examples
#' K2info(K2res)
#' @export
setGeneric("K2info", function(K2res) standardGeneric("K2info"))

#' @export K2info<-
setGeneric("K2info<-", function(K2res, value) standardGeneric("K2info<-"))
setMethod("K2info", "K2", function(K2res) K2res@info)
setMethod("K2info<-", "K2", function(K2res, value) {
    K2res@info <- value
    K2res
})

#' @title Annotated list of partition assignments from K2Taxonomer run
#' @description Retrieve or assign annotated list of subgroups from K2Taxonomer
#' run.
#' @param K2res K2 class object.
#' @return List. A list with an element for each partition.
#' @examples
#' K2results(K2res)
#' @export
setGeneric("K2results", function(K2res) standardGeneric("K2results"))

#' @export K2results<-
setGeneric("K2results<-", function(K2res, value) standardGeneric("K2results<-"))
setMethod("K2results", "K2", function(K2res) K2res@results)
setMethod("K2results<-", "K2", function(K2res, value) {
    K2res@results <- value
    K2res
})

#' @title Named list of genesets used in K2Taxonomer run
#' @description Retrieve or assign named list of genesets for K2Taxonomer run.
#' @param K2res K2 class object.
#' @return Named list. A list with an element for each gene set.
#' @examples
#' K2genesets(K2res)
#' @export
setGeneric("K2genesets", function(K2res) standardGeneric("K2genesets"))

#' @export K2genesets<-
setGeneric("K2genesets<-", function(K2res, value) {
    standardGeneric("K2genesets<-")
})
setMethod("K2genesets", "K2", function(K2res) K2res@genesets)
setMethod("K2genesets<-", "K2", function(K2res, value) {
    K2res@genesets <- value
    K2res
})

#' @title Vector of collapsed pathways for which each gene belongs
#' @description Retrieve or assign named vector of pathways for each gene.
#' @param K2res K2 class object.
#' @return Named vector. A vector with collapse pathway names for each gene.
#' @examples
#' K2gene2Pathway(K2res)
#' @export
setGeneric("K2gene2Pathway", function(K2res) standardGeneric("K2gene2Pathway"))

#' @export K2gene2Pathway<-
setGeneric("K2gene2Pathway<-", function(K2res, value) {
    standardGeneric("K2gene2Pathway<-")
})
setMethod("K2gene2Pathway", "K2", function(K2res) K2res@gene2Pathway)
setMethod("K2gene2Pathway<-", "K2", function(K2res, value) {
    K2res@gene2Pathway <- value
    K2res
})

#' @title ExpressionSet object used in K2Taxonomer run
#' @description Retrieve or assign ExpressionSet object.
#' @param K2res K2 class object.
#' @return ExpressionSet.
#' @examples
#' K2eSet(K2res)
#' @export
setGeneric("K2eSet", function(K2res) standardGeneric("K2eSet"))

#' @export K2eSet<-
setGeneric("K2eSet<-", function(K2res, value) standardGeneric("K2eSet<-"))
setMethod("K2eSet", "K2", function(K2res) K2res@eSet)
setMethod("K2eSet<-", "K2", function(K2res, value) {
    K2res@eSet <- value
    K2res
})

#' @title ExpressionSet object of of GSVA output
#' @description Retrieve or assign ExpressionSet object of GSVA output.
#' @param K2res K2 class object.
#' @return ExpressionSet.
#' @examples
#' K2gSet(K2res)
#' @export
setGeneric("K2gSet", function(K2res) standardGeneric("K2gSet"))

#' @export K2gSet<-
setGeneric("K2gSet<-", function(K2res, value) standardGeneric("K2gSet<-"))
setMethod("K2gSet", "K2", function(K2res) K2res@gSet)
setMethod("K2gSet<-", "K2", function(K2res, value) {
    K2res@gSet <- value
    K2res
})

#' @title Meta data defining parameters of K2Taxonomer run
#' @description Retrieve or assign meta data variables.
#' @param K2res K2 class object.
#' @return Named list. Parameters used in K2Taxonomer run.
#' @examples
#' K2meta(K2res)
#' @export
setGeneric("K2meta", function(K2res) standardGeneric("K2meta"))

#' @export K2meta<-
setGeneric("K2meta<-", function(K2res, value) standardGeneric("K2meta<-"))
setMethod("K2meta", "K2", function(K2res) K2res@meta)
setMethod("K2meta<-", "K2", function(K2res, value) {
    K2res@meta <- value
    K2res
})

#' @title Assign URLs to genes included in K2Taxonomer run
#' @description Retrieve or assign URLs to gene pages.
#' @param K2res K2 class object.
#' @return Named vector. A vector with a URL string for each gene.
#' @examples
#' K2geneURL(K2res)
#' @export
setGeneric("K2geneURL", function(K2res) standardGeneric("K2geneURL"))

#' @export K2geneURL<-
setGeneric("K2geneURL<-", function(K2res, value) standardGeneric("K2geneURL<-"))
setMethod("K2geneURL", "K2", function(K2res) K2res@geneURL)
setMethod("K2geneURL<-", "K2", function(K2res, value) {
    K2res@geneURL <- value
    K2res
})

#' @title Assign URLs to genesets included in K2Taxonomer run
#' @description Retrieve or assign URLs to geneset pages.
#' @param K2res K2 class object.
#' @return Named vector. A vector with a URL string for each geneset.
#' @examples
#' K2genesetURL(K2res)
#' @export
setGeneric("K2genesetURL", function(K2res) standardGeneric("K2genesetURL"))

#' @export K2genesetURL<-
setGeneric("K2genesetURL<-", function(K2res, value) {
    standardGeneric("K2genesetURL<-")
})
setMethod("K2genesetURL", "K2", function(K2res) K2res@genesetURL)
setMethod("K2genesetURL<-", "K2", function(K2res, value) {
    K2res@genesetURL <- value
    K2res
})
