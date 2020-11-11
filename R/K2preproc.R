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
#' @param use Character string. Options are "Z" to generate test statistics or
#' "MEAN" to use means from differential analysis for clustering.
#' @param nFeats "sqrt" or a numeric value <= number of features to subset the
#' features for each partition.
#' @param featMatric Metric to use to assign variance/signal score. Options are
#' 'square' (default) use square values and 'mad' to use MAD scores.
#' @param recalcDataMatrix Logical. Recalculate dataMatrix for each partion?
#' Default is FALSE.
#' @param nBoots A numeric value of the number of bootstraps to run at each
#' split.
#' @param clustFunc Wrapper function to cluster a P x N (See details).
#' @param clustCors Number of cores to use for clustering.
#' @param clustList List of objects to use for clustering procedure.
#' @param linkage Linkage criteria for splitting cosine matrix ('method' in
#' hclust). 'average' by default.
#' @param info A data frame with rownames that match column names in dataMatrix.
#' @param infoClass A named vector denoted types of tests to run on
#' metavariables.
#' @param genesets A named list of features in row names of dataMatrix.
#' @param qthresh A numeric value between 0 and 1 of the FDR cuttoff to define
#' feature sets.
#' @param cthresh A positive value for the coefficient cuttoff to define
#' feature sets.
#' @param ntotal A positive value to use as the background feature count. 20000
#' by default.
#' @param ssGSEAalg A character string, specifying which algorithm to use for
#' running the gsva() function from the GSVA package. Options are "gsva",
#' "ssgsea", "zscore", and "plage". "gsva" by default.
#' @param ssGSEAcores Number of cores to use for running gsva() from the GSVA
#' package. Default is 1.
#' @param oneoff Logical. Allow 1 member clusters?
#' @param stabThresh Threshold for ending clustering.
#' @param geneURL Optional. Named list of URLs to gene information.
#' @param genesetURL Optional. Named list of URLs to geneset information.
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

K2preproc <- function(eSet, cohorts = NULL, vehicle = NULL, covariates = NULL, block = NULL,
    logCounts = FALSE, use = c("Z", "MEAN"), nFeats = "sqrt", featMetric = c("mad",
        "sd", "Sn", "Qn", "F", "square"), recalcDataMatrix = FALSE, nBoots = 500,
    clustFunc = hclustWrapper, clustCors = 1, clustList = list(), linkage = c("mcquitty",
        "ward.D", "ward.D2", "single", "complete", "average", "centroid"), info = NULL,
    infoClass = NULL, genesets = NULL, qthresh = 0.05, cthresh = 0, ntotal = 20000,
    ssGSEAalg = c("gsva", "ssgsea", "zscore", "plage"), ssGSEAcores = 1,
    oneoff = TRUE, stabThresh = 0, geneURL = NULL, genesetURL = NULL) {

    ## Match arguments
    use <- match.arg(use)
    featMetric <- match.arg(featMetric)
    linkage <- match.arg(linkage)
    ssGSEAalg <- match.arg(ssGSEAalg)

    ## Set nFeats if argument == "sqrt"
    if (nFeats == "sqrt") {
        nFeats <- sqrt(nrow(eSet))
    }

    ## Create K2 object from eSet
    K2res <- new("K2", eSet = eSet)

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
    K2meta(K2res) <- list(cohorts = cohorts, vehicle = vehicle, covariates = covariates,
        block = block, logCounts = logCounts, infoClass = infoClass, use = use, nFeats = nFeats,
        featMetric = featMetric, recalcDataMatrix = recalcDataMatrix, nBoots = nBoots,
        clustFunc = clustFunc, clustCors = clustCors, clustList = clustList, linkage = linkage,
        qthresh = qthresh, cthresh = cthresh, ntotal = ntotal, ssGSEAalg = ssGSEAalg,
        ssGSEAcores = ssGSEAcores, oneoff = oneoff, stabThresh = stabThresh)
    
    # Check inputs
    k2Check <- .checkK2(K2res, inputsOnly = TRUE)

    ## Perform differential analysis if cohort information is given
    if (is.null(cohorts)) {

        dataMatrix <- exprs(eSet)

        ## Format info
        if (is.null(info)) {
            info <- pData(eSet)
        }
        info <- data.frame(sampleID = colnames(dataMatrix), info, stringsAsFactors = FALSE)

    } else {

        cat("Collapsing group-level values with LIMMA.\n")
        dataMatrix <- suppressWarnings(.dgeWrapper(eSet, cohorts, vehicle, covariates, use, logCounts = logCounts))

        ## Format info
        if (is.null(info)) {
            info <- pData(eSet)
        }
        info <- info[!duplicated(info[, cohorts]), , drop = FALSE]
        rownames(info) <- info[, cohorts]
        info <- droplevels(info)

    }

    ## Add dataMatrix and info to K2 object
    K2data(K2res) <- dataMatrix
    K2info(K2res) <- info

    ## Check K2 object
    k2Check <- .checkK2(K2res)
    
    return(K2res)
}
