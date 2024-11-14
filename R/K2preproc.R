#' Function to create K2 object for pre-processing
#'
#' This function will generate an object of class, K2.  This will run
#' pre-processing functions for running K2 Taxonomer procedure.
#' @param object An object of class Seurat, SingleCellExperiment, or ExpressionSet.
#' @param cohorts Character. The column in meta data of 'object' that has cohort IDs. Default NULL if no cohorts in data.
#' @param eMatDS Numeric matrix. A matrix with the same number of observations as 'object' containing normalized expression data to be used in analyses downstream of partitioning algorithm.
#' @param vehicle The value in the cohort variable that contains the ID of observation to use as control. Default NULL if no vehicle to be used.
#' @param variables Character. Columns in meta data of 'object' to control for in differential analyses.
#' @param seuAssay Character. Name of assay in Seurat object containing expression data for running partitioning algorithm. If cohorts based on clustering, this should be the assay used.
#' @param seuAssayDS Character. Name of assay in Seurat object containing expression data normalized expression data to be used in analyses downstream of partitioning algorithm.
#' @param sceAssay Character. Name of assay in SingleCellExperimen object containing expression data for running partitioning algorithm. If cohorts based on clustering, this should be the assay used.
#' @param sceAssayDS Character. Name of assay in SingleCellExperiment object containing expression data normalized expression data to be used in analyses downstream of partitioning algorithm.
#' @param block Character. Block parameter in limma for modeling higherarchical data structure, such as multiple observations per individual.
#' @param logCounts Logical. Whether or not expression values are log-scale counts or log normalized counts from RNA-seq. Default is TRUE.
#' @param use Character. Options are 'Z' to generate test statistics or 'MEAN' to use means from differential analysis for clustering.
#' @param nFeats 'sqrt' or a numeric value <= number of features to subset for each partition.
#' @param featMetric Character. Metric to use to assign gene-level variance/signal score.
#' \itemize{
#'  \item{F: }{F-statistic from evaluating differences in means across cohort}
#'  \item{mad: }{Median absolute deviation}
#'  \item{sd: }{Standard deviation}
#'  \item{Sn: }{Robust scale estimator}
#' }
#' @param DGEmethod Character. Method for running differential gene expression analyses. Use one of either 'limma' (default) or 'mast'.
#' @param DGEexpThreshold Numeric. A value between 0 and 1 indicating for filtering lowly expressed genes for partition-specific differential gene expression. Proportion of observations with counts > 0 in at least one subgroup at a specific partition.
#' @param recalcDataMatrix Logical. Recalculate dataMatrix for each partion? Default is FALSE.
#' @param nBoots nBoots A value of the number of bootstraps to run at each partition. Default is 500.
#' @param clustFunc Character. Wrapper function to be used in recursive partitioning.
#' \itemize{
#'  \item{cKmeansDownsampleSqrt: }{Perform constrained K-means clustering after
#'  subsampling each cohort by the square root of the number of observations}
#'  \item{cKmeansDownsampleSmallest: }{Perform constrained K-means clustering
#'  after subsampling each cohort by the size of the smallest cohort}
#'  \item{hclustWrapper: }{Perform hierarchical clustering}
#' }
#' @param useCors Numeric. Number of cores to use for parallelizable processes.
#' @param clustList Optional named list of parameters to use with clustFunc.
#' \itemize{
#'  \item{cKmeansDownsampleSqrt: }{
#'    \itemize{
#'      \item{maxIter: }{The maximum number of iterations to use with lcvqe()}
#'    }
#'  }
#'  \item{cKmeansDownsampleSmallest: }{
#'    \itemize{
#'      \item{maxIter: }{The maximum number of iterations to use with lcvqe()}
#'    }
#'  }
#'  \item{hclustWrapper: }{
#'    \itemize{
#'      \item{aggMethod: }{One of the hierarchichal methods specified by hclust() function}
#'      \item{distMetric: }{One of the distance metrics specified by dist() function}
#'    }
#'  }
#' }
#' @param linkage Character. Linkage criteria for splitting cosine matrix ('method' in hclust). 'average' by default.
#' @param genesets Named list. Feature sets to be includes in enrichment-based analyses.
#' @param qthresh Numeric. A value between 0 and 1 indicating the FDR cuttoff to define feature sets.
#' @param cthresh Numeric. A positive value for the coefficient cuttoff to define feature sets.
#' @param ntotal Numeric. A positive value to use as the background feature count. 20000 by default.
#' @param ScoreGeneSetMethod Character. Method for gene set scoring. Use one of either 'GSVA' (default) or 'AUCELL'.
#' @param oneoff Logical. Allow 1 observation partition groups? Default is TRUE.
#' @param stabThresh Numeric. A value between 0 and 1 indicatingThreshold for ending clustering.
#' @param info Character. A vector of column names in meta data of 'object' that contain information to be used in cohort annotation of dashboard visualization
#' @param geneURL Named list. URLs linking genes to external resources.
#' @param genesetURL Named list. URLs linking gene set to external resources.
#' @return An object of class, `K2`.
#' @references
#' \insertRef{reed_2020}{K2Taxonomer}
#' \insertRef{limma}{K2Taxonomer}
#' \insertRef{Rousseeuw1993}{K2Taxonomer}
#' @keywords clustering
#' @export
#' @import limma
#' @import Biobase

K2preproc <- function(object, cohorts=NULL, eMatDS = NULL, vehicle=NULL, variables=NULL,
    seuAssay = "RNA", seuAssayDS = "RNA",
    sceAssay = "logcounts", sceAssayDS = NULL,
    block=NULL, logCounts=TRUE, use=c("Z", "MEAN"), nFeats="sqrt",
    featMetric=c("F", "mad", "sd", "Sn"),
    DGEmethod = c("limma", "mast"), DGEexpThreshold = 0.25, 
    recalcDataMatrix=TRUE, nBoots=500, useCors=1,
    clustFunc="cKmeansDownsampleSqrt", clustList=NULL, linkage=c("mcquitty", "ward.D",
    "ward.D2", "single", "complete", "average", "centroid"), info=NULL, 
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
        variables=variables, block=block, logCounts=logCounts,
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
