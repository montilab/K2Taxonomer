#' @export K2data
#' @export K2data<-
setGeneric("K2data", function(x) standardGeneric("K2data"))
setGeneric("K2data<-", function(x, value) standardGeneric("K2data<-"))
setMethod("K2data", "K2", function(x) x@dataMatrix)
setMethod("K2data<-", "K2", function(x, value) {
  x@dataMatrix <- value
  x
})

#' @export K2info
#' @export K2info<-
setGeneric("K2info", function(x) standardGeneric("K2info"))
setGeneric("K2info<-", function(x, value) standardGeneric("K2info<-"))
setMethod("K2info", "K2", function(x) x@info)
setMethod("K2info<-", "K2", function(x, value) {
  x@info <- value
  x
})

#' @export K2results
#' @export K2results<-
setGeneric("K2results", function(x) standardGeneric("K2results"))
setGeneric("K2results<-", function(x, value) standardGeneric("K2results<-"))
setMethod("K2results", "K2", function(x) x@results)
setMethod("K2results<-", "K2", function(x, value) {
  x@results <- value
  x
})

#' @export K2genesets
#' @export K2genesets<-
setGeneric("K2genesets", function(x) standardGeneric("K2genesets"))
setGeneric("K2genesets<-", function(x, value) standardGeneric("K2genesets<-"))
setMethod("K2genesets", "K2", function(x) x@genesets)
setMethod("K2genesets<-", "K2", function(x, value) {
  x@genesets <- value
  x
})

#' @export K2gene2Pathway
#' @export K2gene2Pathway<-
setGeneric("K2gene2Pathway", function(x) standardGeneric("K2gene2Pathway"))
setGeneric("K2gene2Pathway<-", function(x, value) standardGeneric("K2gene2Pathway<-"))
setMethod("K2gene2Pathway", "K2", function(x) x@gene2Pathway)
setMethod("K2gene2Pathway<-", "K2", function(x, value) {
  x@gene2Pathway <- value
  x
})

#' @export K2eSet
#' @export K2eSet<-
setGeneric("K2eSet", function(x) standardGeneric("K2eSet"))
setGeneric("K2eSet<-", function(x, value) standardGeneric("K2eSet<-"))
setMethod("K2eSet", "K2", function(x) x@eSet)
setMethod("K2eSet<-", "K2", function(x, value) {
  x@eSet <- value
  x
})

#' @export K2gSet
#' @export K2gSet<-
setGeneric("K2gSet", function(x) standardGeneric("K2gSet"))
setGeneric("K2gSet<-", function(x, value) standardGeneric("K2gSet<-"))
setMethod("K2gSet", "K2", function(x) x@gSet)
setMethod("K2gSet<-", "K2", function(x, value) {
  x@gSet <- value
  x
})

#' @export K2meta
#' @export K2meta<-
setGeneric("K2meta", function(x) standardGeneric("K2meta"))
setGeneric("K2meta<-", function(x, value) standardGeneric("K2meta<-"))
setMethod("K2meta", "K2", function(x) x@meta)
setMethod("K2meta<-", "K2", function(x, value) {
  x@meta <- value
  x
})

#' @export K2geneURL
#' @export K2geneURL<-
setGeneric("K2geneURL", function(x) standardGeneric("K2geneURL"))
setGeneric("K2geneURL<-", function(x, value) standardGeneric("K2geneURL<-"))
setMethod("K2geneURL", "K2", function(x) x@geneURL)
setMethod("K2geneURL<-", "K2", function(x, value) {
  x@geneURL <- value
  x
})

#' @export K2genesetURL
#' @export K2genesetURL<-
setGeneric("K2genesetURL", function(x) standardGeneric("K2genesetURL"))
setGeneric("K2genesetURL<-", function(x, value) standardGeneric("K2genesetURL<-"))
setMethod("K2genesetURL", "K2", function(x) x@genesetURL)
setMethod("K2genesetURL<-", "K2", function(x, value) {
  x@genesetURL <- value
  x
})



