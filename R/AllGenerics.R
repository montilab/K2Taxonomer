#' @title Data matrix from K2 object
#' @description Retrieve or assign matrix of data used in K2Taxonomer run.
#' @param K2res K2 class object.
#' @param value Matrix. Numeric data used for K2Taxonomer run.
#' @return Matrix. Numeric data used for K2Taxonomer run.
#' @examples
#' data(K2res)
#' K2data(K2res)
#'
#' @export
setGeneric("K2data", function(K2res) standardGeneric("K2data"))

#' @rdname K2data
setGeneric("K2data<-", function(K2res, value) standardGeneric("K2data<-"))

#' @rdname K2data
setMethod("K2data", "K2", function(K2res) K2res@dataMatrix)

#' @rdname K2data
setMethod("K2data<-", "K2", function(K2res, value) {
    K2res@dataMatrix <- value
    K2res
})

#' @title Data frame of phenotypic information from K2 object
#' @description Retrieve or assign observation information in K2Taxonomer run.
#' @param K2res K2 class object.
#' @param value Data frame. Observation information from K2Taxonomer run.
#' @return Data frame. Observation information from K2Taxonomer run.
#' @examples
#' data(K2res)
#' head(K2colData(K2res))
#'
#' @export
setGeneric("K2colData", function(K2res) standardGeneric("K2colData"))

#' @rdname K2colData
setGeneric("K2colData<-", function(K2res, value) standardGeneric("K2colData<-"))

#' @rdname K2colData
setMethod("K2colData", "K2", function(K2res) K2res@colData)

#' @rdname K2colData
setMethod("K2colData<-", "K2", function(K2res, value) {
    K2res@colData <- value
    K2res
})

#' @title Annotated list of partition assignments from K2Taxonomer run
#' @description Retrieve or assign annotated list of subgroups from K2Taxonomer
#' run.
#' @param K2res K2 class object.
#' @param value List. A list with an element for each partition.
#' @return List. A list with an element for each partition.
#' @examples
#' data(K2res)
#' K2resList <- K2results(K2res)
#'
#' @export
setGeneric("K2results", function(K2res) standardGeneric("K2results"))

#' @rdname K2results
setGeneric("K2results<-", function(K2res, value) standardGeneric("K2results<-"))

#' @rdname K2results
setMethod("K2results", "K2", function(K2res) K2res@results)

#' @rdname K2results
setMethod("K2results<-", "K2", function(K2res, value) {
    K2res@results <- value
    K2res
})

#' @title Named list of genesets used in K2Taxonomer run
#' @description Retrieve or assign named list of genesets for K2Taxonomer run.
#' @param K2res K2 class object.
#' @param value Named list. A list with an element for each ge
#' @return Named list. A list with an element for each gene set.
#' @examples
#' data(K2res)
#' head(K2genesets(K2res))
#'
#' @export
setGeneric("K2genesets", function(K2res) standardGeneric("K2genesets"))

#' @rdname K2genesets
setGeneric("K2genesets<-", function(K2res, value) {
    standardGeneric("K2genesets<-")
})

#' @rdname K2genesets
setMethod("K2genesets", "K2", function(K2res) K2res@genesets)

#' @rdname K2genesets
setMethod("K2genesets<-", "K2", function(K2res, value) {
    K2res@genesets <- value
    K2res
})

#' @title Vector of collapsed pathways for which each gene belongs
#' @description Retrieve or assign named vector of pathways for each gene.
#' @param K2res K2 class object.
#' @param value Named vector. A vector with collapse pathway names for each
#' gene.
#' @return Named vector. A vector with collapse pathway names for each gene.
#' @examples
#' data(K2res)
#' head(K2gene2Pathway(K2res))
#'
#' @export
setGeneric("K2gene2Pathway", function(K2res) standardGeneric("K2gene2Pathway"))

#' @rdname K2gene2Pathway
setGeneric("K2gene2Pathway<-", function(K2res, value) {
    standardGeneric("K2gene2Pathway<-")
})

#' @rdname K2gene2Pathway
setMethod("K2gene2Pathway", "K2", function(K2res) K2res@gene2Pathway)

#' @rdname K2gene2Pathway
setMethod("K2gene2Pathway<-", "K2", function(K2res, value) {
    K2res@gene2Pathway <- value
    K2res
})

#' @title Expression matrix object used in K2Taxonomer run
#' @description Retrieve or assign Matrix object.
#' @param K2res K2 class object.
#' @param value Matrix.
#' @return Matrix.
#' @examples
#' data(K2res)
#' K2eMat(K2res)
#'
#' @export
setGeneric("K2eMat", function(K2res) standardGeneric("K2eMat"))

#' @rdname K2eMat
setGeneric("K2eMat<-", function(K2res, value) standardGeneric("K2eMat<-"))

#' @rdname K2eMat
setMethod("K2eMat", "K2", function(K2res) K2res@eMat)

#' @rdname K2eMat
setMethod("K2eMat<-", "K2", function(K2res, value) {
    K2res@eMat <- value
    K2res
})

#' @title Expression matrix object used for differential expression
#' @description Retrieve or assign Matrix object.
#' @param K2res K2 class object.
#' @param value Matrix.
#' @return Matrix.
#' @examples
#' data(K2res)
#' K2eMatDS(K2res)
#'
#' @export
setGeneric("K2eMatDS", function(K2res) standardGeneric("K2eMatDS"))

#' @rdname K2eMatDS
setGeneric("K2eMatDS<-", function(K2res, value) standardGeneric("K2eMatDS<-"))

#' @rdname K2eMatDS
setMethod("K2eMatDS", "K2", function(K2res) K2res@eMatDS)

#' @rdname K2eMatDS
setMethod("K2eMatDS<-", "K2", function(K2res, value) {
  K2res@eMatDS <- value
  K2res
})

#' @title Matrix object of of GSVA output
#' @description Retrieve or assign Matrix object of GSVA output.
#' @param K2res K2 class object.
#' @param value Matrix.
#' @return Matrix.
#' @examples
#' data(K2res)
#' K2gMat(K2res)
#'
#' @export
setGeneric("K2gMat", function(K2res) standardGeneric("K2gMat"))

#' @rdname K2gMat
setGeneric("K2gMat<-", function(K2res, value) standardGeneric("K2gMat<-"))

#' @rdname K2gMat
setMethod("K2gMat", "K2", function(K2res) K2res@gMat)

#' @rdname K2gMat
setMethod("K2gMat<-", "K2", function(K2res, value) {
    K2res@gMat <- value
    K2res
})

#' @title Meta data defining parameters of K2Taxonomer run
#' @description Retrieve or assign meta data variables.
#' @param K2res K2 class object.
#' @param value Named list. Parameters used in K2Taxonomer run.
#' @return Named list. Parameters used in K2Taxonomer run.
#' @examples
#' data(K2res)
#' K2meta(K2res)
#'
#' @export
setGeneric("K2meta", function(K2res) standardGeneric("K2meta"))

#' @rdname K2meta
setGeneric("K2meta<-", function(K2res, value) standardGeneric("K2meta<-"))

#' @rdname K2meta
setMethod("K2meta", "K2", function(K2res) K2res@meta)

#' @rdname K2meta
setMethod("K2meta<-", "K2", function(K2res, value) {
    K2res@meta <- value
    K2res
})

#' @title Assign URLs to genes included in K2Taxonomer run
#' @description Retrieve or assign URLs to gene pages.
#' @param K2res K2 class object.
#' @param value Named vector. A vector with a URL string for each gene.
#' @return Named vector. A vector with a URL string for each gene.
#' @examples
#' data(K2res)
#' K2geneURL(K2res)
#'
#' @export
setGeneric("K2geneURL", function(K2res) standardGeneric("K2geneURL"))

#' @rdname K2geneURL
setGeneric("K2geneURL<-", function(K2res, value) standardGeneric("K2geneURL<-"))

#' @rdname K2geneURL
setMethod("K2geneURL", "K2", function(K2res) K2res@geneURL)

#' @rdname K2geneURL
setMethod("K2geneURL<-", "K2", function(K2res, value) {
    K2res@geneURL <- value
    K2res
})

#' @title Assign URLs to genesets included in K2Taxonomer run
#' @description Retrieve or assign URLs to geneset pages.
#' @param K2res K2 class object.
#' @param value Named vector. A vector with a URL string for each geneset.
#' @return Named vector. A vector with a URL string for each geneset.
#' @examples
#' data(K2res)
#' K2genesetURL(K2res)
#' @export
setGeneric("K2genesetURL", function(K2res) standardGeneric("K2genesetURL"))

#' @rdname K2genesetURL
setGeneric("K2genesetURL<-", function(K2res, value) {
    standardGeneric("K2genesetURL<-")
})

#' @rdname K2genesetURL
setMethod("K2genesetURL", "K2", function(K2res) K2res@genesetURL)

#' @rdname K2genesetURL
setMethod("K2genesetURL<-", "K2", function(K2res, value) {
    K2res@genesetURL <- value
    K2res
})
