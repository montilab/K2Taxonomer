#' Wrapper function to run K2Taxonomer algorithm and annotation
#'
#' This function will generate an object of class, K2.  This will run the K2
#' Taxonomer procedure, differential analysis, and finally hyperenrichment on a
#' named list of feature sets.
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
#' @param featMetric Metric to use to assign variance/signal score. Options are:
#' 'mad' (default), 'mad', 'Sn', 'Qn', 'F', and 'square'.
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
#' running the gsva() function from the GSVA package. Options are 'gsva',
#' 'ssgsea', 'zscore', and 'plage'. 'gsva' by default.
#' @param ssGSEAcores Number of cores to use for running gsva() from the GSVA
#' package. Default is 1.
#' @param oneoff Logical. Allow 1 member clusters?
#' @param stabThresh Threshold for ending clustering.
#' @param geneURL Optional. Named list of URLs to gene information.
#' @param genesetURL Optional. Named list of URLs to geneset information.
#' @return An object of class, `K2`.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#' @export
#' @import limma
#' @import Biobase
#' @import GSVA
#' @examples
#' ## Read in ExpressionSet object
#' library(Biobase)
#' data(sample.ExpressionSet)
#'
#' ## Create dummy set of gene sets
#' genes <- rownames(sample.ExpressionSet)
#' genesetsMadeUp <- list(
#'     GS1=genes[1:50],
#'     GS2=genes[51:100],
#'     GS3=genes[101:150])
#'
#' ## Run K2 Taxonomer wrapper
#' K2Res <- runK2Taxonomer(sample.ExpressionSet,
#'                     genesets=genesetsMadeUp,
#'                     qthresh=0.1,
#'                     ssGSEAalg='gsva',
#'                     ssGSEAcores=1,
#'                     stabThresh=0.5)
#'

runK2Taxonomer <- function(eSet, cohorts=NULL, vehicle=NULL,
    covariates=NULL, block=NULL, logCounts=FALSE, use=c("Z",
        "MEAN"), nFeats="sqrt", featMetric=c("mad", "sd",
        "Sn", "Qn", "F", "square"), recalcDataMatrix=FALSE,
    nBoots=500, clustFunc=hclustWrapper, clustCors=1, clustList=list(),
    linkage=c("mcquitty", "ward.D", "ward.D2", "single", "complete",
        "average", "centroid"), info=NULL, infoClass=NULL,
    genesets=NULL, qthresh=0.05, cthresh=0, ntotal=20000,
    ssGSEAalg=c("gsva", "ssgsea", "zscore", "plage"), ssGSEAcores=1,
    oneoff=TRUE, stabThresh=0, geneURL=NULL, genesetURL=NULL) {

    ## Match arguments
    K2res <- K2preproc(eSet, cohorts=cohorts, vehicle=vehicle,
        covariates=covariates, block=block, logCounts=logCounts,
        use=use, nFeats=nFeats, featMetric=featMetric,
        recalcDataMatrix=recalcDataMatrix, nBoots=nBoots,
        clustFunc=clustFunc, clustCors=clustCors, clustList=clustList,
        linkage=linkage, info=info, infoClass=infoClass,
        genesets=genesets, qthresh=qthresh, cthresh=cthresh,
        ntotal=ntotal, ssGSEAalg=ssGSEAalg, ssGSEAcores=ssGSEAcores,
        oneoff=oneoff, stabThresh=stabThresh, geneURL=geneURL,
        genesetURL=genesetURL)

    ## Cluster data and create S4 object
    cat("Running K2Taxonomer bootstraps.\n")
    K2res <- K2tax(K2res)

    ## Check clusters for significance of specified variables
    if (!is.null(infoClass)) {
        cat("Checking clusters for significant associations with
        meta-variables.\n")
        K2res <- runTestsMods(K2res)
    }

    ## Add differential analysis results
    cat("Running differential analysis on genes.\n")
    K2res <- runDGEmods(K2res)

    if (!is.null(genesets)) {

        ## Perform hyperenrichment analysis
        cat("Running hyperenrichment analysis on genes with FDR<",
            qthresh, " & |coefficient|>", cthresh, ".\n", sep="")
        K2res <- runGSEmods(K2res)

        ## Perform ssGSEA or ssGVA
        cat("Running ssGSEA.\n")
        K2res <- runGSVAmods(K2res, min.sz=2)

        ## Perform differential analysis on gene sets
        cat("Running differential analysis on gene sets.\n")
        K2res <- runDSSEmods(K2res)

    }

    return(K2res)

}
