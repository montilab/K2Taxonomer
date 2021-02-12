#' Wrapper for constrained K-means on subsampled data
#'
#' This fuction is a wrapper for the constrained Kmeans algorithm using
#' `lcvqe` from the `conclust` package. This function will subset each
#' cohort down to that with the smallest number of observations.This 
#' function is not meant to be run individually, but as a 'clustFunc' 
#' argument for running `K2preproc()`, `runK2Taxonomer()`, and 
#' `K2tax()`.
#'
#' @param dataMatrix An P x C numeric matrix of data. Where C is the number of cohort labels. 
#' @param clustList List of objects to use for clustering procedure.
#' \itemize{
#'  \item{"eMat"}{P x N Expression matrix of data set. See `?Biobase::exprs()`.}
#'  \item{"labs"}{Vector of cohort labels of observations in data set, corresponding to columns of `clustList$eMat`.}
#' }
#' @return A character string of concatenated 1's and 2's pertaining to the
#' cluster assignment of each column in dataMatrix.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{cKm}{K2Taxonomer}
#' @export
#' @import conclust
#' @examples
#' 
#' dat <- scRNAseq::ReprocessedAllenData(assays = "rsem_tpm")[seq_len(50),]
#' 
#' eSet <- ExpressionSet(assayData = assay(dat))
#' pData(eSet) <- as.data.frame(colData(dat))
#' exprs(eSet) <- log2(exprs(eSet) + 1)
#' 
#' ## Subset for fewer cluster labels for this example
#' eSet <- eSet[, !is.na(eSet$Primary.Type) & 
#'                  eSet$Primary.Type %in% c("L4 Arf5", 
#'                      "L4 Ctxn3", "L4 Scnn1a", "L5 Ucma", "L5a Batf3")]
#' 
#' ## Create cell type variable with spaces
#' eSet$celltype <- gsub(" ", "_", eSet$Primary.Type)
#' 
#' ## Create clustList
#' cL <- list(
#'     eMat = exprs(eSet),
#'     labs = eSet$celltype,
#'     maxIter = 10
#' )
#' 
#' ## Run K2preproc to generate generate data matrix
#' ## with a column for each celltype.
#' K2res <- K2preproc(eSet,
#'                    cohorts = "celltype",
#'                    featMetric = "F",
#'                    logCounts = TRUE)
#' dm <- K2data(K2res)
#' 
#' ## Generate K=2 split with constrained K-means
#' cKmeansWrapperSubsample(dm, cL)
#' 

## Create wrapper to subsample
cKmeansWrapperSubsample <- function(dataMatrix, clustList) {
    
    clustList$labs <- as.character(clustList$labs)
    
    eMatSub <- clustList$eMat[rownames(dataMatrix), clustList$labs %in% colnames(dataMatrix)]
    labsSub <- clustList$labs[clustList$labs %in% colnames(dataMatrix)]
    
    ## Subsample the data
    minSize <- min(table(labsSub))
    sVec <- unlist(lapply(unique(labsSub), function(x, minSize) {
        sample(which(labsSub == x), minSize)
    }, minSize))
    eMatSub <- eMatSub[, sVec]
    labsSub <- labsSub[sVec]
    
    ## Get constraints
    mustLink <- outer(labsSub, labsSub, "==")
    mustLink[upper.tri(mustLink, diag = TRUE)] <- FALSE
    mustLink <- which(mustLink, arr.ind = TRUE)
    
    ## Cluster data
    dClust = factor(lcvqe(t(eMatSub), k = 2, mustLink = mustLink, cantLink = matrix(c(1, 
        1), nrow = 1), maxIter = clustList$maxIter), levels = c(1, 2))
    
    ## Get label-level clusters
    dMat <- as.matrix(table(dClust, labsSub))[, colnames(dataMatrix)]
    modVec <- apply(dMat, 2, which.max)
    mods <- paste(modVec, collapse = "")
    
    return(mods)
}
