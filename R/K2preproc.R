#' Function to create K2 object for pre-processing
#'
#' This function will generate an object of class, K2.  This will run
#' pre-processing functions for running K2 Taxonomer procedure.
#' @param eSet An expression set object.
#' @param cohorts The column in phenotype data of eSet that has cohort ID's.
#' Default NULL if no pre-processing of data.
#' @param vehicle The value in the cohort variable that contains the vehicle
#' ID. Default NULL if no vehicle to be used.
#' @param covariates Covariates in phenotype data of eSet to control for in
#' differential analysis.
#' @param block Block parameter in limma for modelling random-like effects.
#' @param logCounts Logical. Whether or not expression values are log-scale
#' counts or log normalized counts from RNA-seq. Default is FALSE.
#' @param use Character string. Options are 'Z' to generate test statistics or
#' 'MEAN' to use means from differential analysis for clustering.
#' @param nFeats 'sqrt' or a numeric value <= number of features to subset the
#' features for each partition.
#' @param featMetric Metric to use to assign variance/signal score. Options are
#' 'square' (default) use square values and 'mad' to use MAD scores.
#' @param recalcDataMatrix Logical. Recalculate dataMatrix for each partion?
#' Default is FALSE.
#' @param nBoots A numeric value of the number of bootstraps to run at each
#' split.
#' @param clustFunc Wrapper function to cluster a P x N (See details).
#' @param useCors Number of cores to use for clustering.
#' @param clustList List of objects to use for clustering procedure.
#' @param linkage Linkage criteria for splitting cosine matrix ('method' in
#' hclust). 'average' by default.
#' @param info A data frame with rownames that match column names in dataMatrix.

#' @param genesets A named list of features in row names of dataMatrix.
#' @param qthresh A numeric value between 0 and 1 of the FDR cuttoff to define
#' feature sets.
#' @param cthresh A positive value for the coefficient cuttoff to define
#' feature sets.
#' @param ntotal A positive value to use as the background feature count. 20000
#' by default.
#' @param oneoff Logical. Allow 1 member clusters?
#' @param stabThresh Threshold for ending clustering.
#' @param geneURL Optional. Named list of URLs to gene information.
#' @param genesetURL Optional. Named list of URLs to geneset information.
#' @return An object of class, `K2`.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{limma}{K2Taxonomer}
#' @keywords clustering
#' @export
#' @import limma
#' @import Biobase
#' @examples
#' ## Read in ExpressionSet object
#' library(Biobase)
#' data(sample.ExpressionSet)
#'
#' ## Pre-process and create K2 object
#' K2res <- K2preproc(sample.ExpressionSet)
#'

K2preproc <- function(object, cohorts=NULL, eMatDS = NULL, vehicle=NULL, covariates=NULL,
    seuAssay = "RNA", seuAssayDS = "RNA",
    sceAssay = "logcounts", sceAssayDS = NULL,
    block=NULL, logCounts=FALSE, use=c("Z", "MEAN"), nFeats="sqrt",
    featMetric=c("mad", "sd", "Sn", "Qn", "F", "square"),
    DGEmethod = c("limma", "mast"), DGEexpThreshold = 0.25, 
    recalcDataMatrix=TRUE, nBoots=500, clustFunc="cKmeansDownsampleSqrt",
    useCors=1, clustList=NULL, linkage=c("mcquitty", "ward.D",
    "ward.D2", "single", "complete", "afverage", "centroid"), info=NULL, 
    genesets=NULL, qthresh=0.05, cthresh=0, ntotal=20000, 
    ScoreGeneSetMethod = c("GSVA", "AUCELL"), oneoff=TRUE, 
    stabThresh=0, geneURL=NULL, genesetURL=NULL) {

    ## Match arguments
    use <- match.arg(use)
    featMetric <- match.arg(featMetric)
    linkage <- match.arg(linkage)
    ScoreGeneSetMethod <- match.arg(ScoreGeneSetMethod)
    DGEmethod <- match.arg(DGEmethod)

    ## Initialize K2 object with expression data
    if(!class(object) %in% c("Seurat", "ExpressionSet", "SingleCellExperiment")) {
      stop("'object' must be one of classes 'Seurat', 'ExpressionSet', 'SingleCellExperiment'")
    }
    
    ## Set cluster function
    if(class(clustFunc) == "character") {
      clustFunc <- get(clustFunc)
    }
    
    ## Initialize downstream expression matrix
    if(is.null(eMatDS)) {
      eMatDS <- new("matrix")
    }
    
    ## Initialize clustList
    if(is.null(clustList)) {
      clustList <- list()
    }
    
    if(class(object) == "ExpressionSet") {
      K2res <- new("K2", 
                   eMat=object@assayData$exprs,
                   eMatDS=eMatDS,
                   colData=object@phenoData@data)
    }
    
    if(class(object) == "Seurat") {
      if(nrow(eMatDS) == 0) {
        K2res <- new("K2", 
                     eMat=object@assays[[seuAssay]]@scale.data,
                     eMatDS=object@assays[[seuAssayDS]]@data,
                     colData=object@meta.data)
      } else {
        K2res <- new("K2", 
                     eMat=object@assays[[seuAssay]]@scale.data,
                     eMatDS=eMatDS,
                     colData=object@meta.data)
      }
    
    }
    
    if(class(object) == "SingleCellExperiment") {
      if(is.null(sceAssayDS)) {
        K2res <- new("K2", 
                     eMat=object@assays@data[[sceAssay]],
                     eMatDS=eMatDS,
                     colData=as.data.frame(object@colData))
      } else {
        K2res <- new("K2", 
                     eMat=object@assays@data[[sceAssay]],
                     eMatDS=object@assays@data[[sceAssayDS]],
                     colData=as.data.frame(object@colData))
      }
    }
    
    ## Set nFeats if argument == 'sqrt'
    if (nFeats == "sqrt") {
      nFeats <- ceiling(sqrt(nrow(K2eMat(K2res))))
    }

    ## Add genesets to K2res
    if (!is.null(genesets)) {

        K2genesets(K2res) <- genesets

        ## Add genesets and pathway maps to analysis
        K2gene2Pathway(K2res) <- getGenePathways(genesets)
    }

    ## Add URLs
    if (!is.null(geneURL)) {
        K2geneURL(K2res) <- geneURL
    }
    if (!is.null(genesetURL)) {
        K2genesetURL(K2res) <- genesetURL
    }

    ## Add meta information
    K2meta(K2res) <- list(cohorts=cohorts, vehicle=vehicle,
        covariates=covariates, block=block, logCounts=logCounts,
        use=use, nFeats=nFeats, featMetric=featMetric,
        DGEmethod=DGEmethod, DGEexpThreshold=DGEexpThreshold, 
        recalcDataMatrix=recalcDataMatrix, nBoots=nBoots,
        clustFunc=clustFunc, useCors=useCors, clustList=clustList,
        linkage=linkage, qthresh=qthresh, cthresh=cthresh,
        ntotal=ntotal, ScoreGeneSetMethod=ScoreGeneSetMethod,
        oneoff=oneoff, stabThresh=stabThresh, info=info)

    # Check inputs
    k2Check <- .checkK2(K2res, inputsOnly=TRUE)

    ## Perform differential analysis if cohort information is
    ## given
    if (is.null(K2meta(K2res)$cohorts)) {

        dataMatrix <- K2eMat(K2res)

    } else {

        dataMatrix <- suppressWarnings(.dgeWrapper(K2res))

    }
    
    ## Add dataMatrix and info to K2 object
    K2data(K2res) <- dataMatrix
    
    ## Format sample/cohort information
    K2meta(K2res)$info <- .formatInfo(K2res)

    ## Check K2 object
    k2Check <- .checkK2(K2res)

    return(K2res)
}
