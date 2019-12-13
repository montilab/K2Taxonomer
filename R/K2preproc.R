#' Function to create K2 object for pre-processing
#'
#' This function will generate an object of class, K2.  This will run pre-processing functions for running K2 Taxonomer procedure.
#' @param eSet An expression set object.
#' @param cohorts The column in phenotype data of eSet that has cohort ID's. Default NULL if no pre-processing of data.
#' @param vehicle The value in the cohort variable that contains the vehicle ID. Default NULL if no vehicle to be used.
#' @param covariates Covariates in phenotype data of eSet to control for in differential analysis.
#' @param use Generate test statistics or use means from differential analysis for clustering.
#' @param nFeats A numeric value <= P of subsets of the data to use.
#' @param featMatric Metric to use to assign variance/signal score. Options are
#' "square" (default) use square values and "mad" to use MAD scores.
#' @param nBoots A numeric value of the number of bootstraps to run at each split.
#' @param clustFunc Wrapper function to cluster a P x N (See details).
#' @param linkage Linkage criteria for splitting cosine matrix ("method" in hclust). "average" by default.
#' @param info A data frame with rownames that match column names in dataMatrix.
#' @param infoClass = A named vector denoted types of tests to run on metavariables.
#' @param genesets A named list of features in row names of dataMatrix.
#' @param qthresh A numeric value between 0 and 1 of the FDR cuttoff to define feature sets.
#' @param cthresh A positive value for the coefficient cuttoff to define feature sets.
#' @param oneoff Logical. Allow 1 member clusters?
#' @param stabThresh Threshold for ending clustering.
#' @param geneURL Optional. Named list of URLs to gene information.
#' @param genesetURL Optional. Named list of URLs to geneset information.
#' @keywords clustering
#' @export
#' @import limma
#' @import Biobase
#' @examples
#' K2preproc(eSet)

K2preproc <- function(eSet,
                         cohorts = NULL,
                         vehicle = NULL,
                         covariates = NULL,
                         use = c("Z", "MEAN"),
                         nFeats = nrow(eSet)*0.02,
                         featMetric = c("square", "mad"),
                         nBoots = 500,
                         clustFunc = hclust_wrapper,
                         linkage = "mcquitty",
                         info = NULL,
                         infoClass = NULL,
                         genesets = NULL,
                         qthresh = 0.05,
                         cthresh = 0,
                         ntotal = 20000,
                         ssGSEAalg = c("ssgsea", "gsva", "zscore", "plage"),
                         ssGSEAminQvalue = 0.01,
                         ssGSEAcores = 0, 
                         oneoff = TRUE,
                         stabThresh = 0,
                         geneURL = NULL,
                         genesetURL = NULL
){
  
  # Match arguments
  use <- match.arg(use)
  featMetric <- match.arg(featMetric)
  ssGSEAalg <- match.arg(ssGSEAalg)

  # Create K2 object from eSet
  K2res <- new("K2",
               eSet = eSet)
  
  # Add meta information
  K2meta(K2res) <- list(cohorts = cohorts,
                        vehicle = vehicle,
                        covariates = covariates,
                        infoClass = infoClass,
                        use = use,
                        nFeats = nFeats,
                        featMetric = featMetric,
                        nBoots = nBoots,
                        clustFunc = clustFunc,
                        linkage = linkage,
                        qthresh = qthresh,
                        cthresh = cthresh,
                        ntotal = ntotal,
                        ssGSEAalg = ssGSEAalg,
                        ssGSEAminQvalue = ssGSEAminQvalue,
                        ssGSEAcores = ssGSEAcores,
                        oneoff = oneoff,
                        stabThresh = stabThresh)
  
  # Perform differential analysis if cohort information is given
  if  (is.null(cohorts)){
    
    dataMatrix <- exprs(eSet)
    
    # Format info
    if (is.null(info)) {
      info <- data.frame(row.names = colnames(dataMatrix))
    }
    info <- data.frame(sampleID = colnames(dataMatrix),
                       info,
                       stringsAsFactors = FALSE)
    
  } else {
    
    cat("Collapsing group-level values with LIMMA.\n")
    dataMatrix <- .dge_wrapper(eSet, cohorts, vehicle, covariates, use)
    
    # Format info
    if (is.null(info)) {
      info <- data.frame(pData(eSet)[,cohorts], row.names = colnames(eSet))
      colnames(info) <- cohorts
    }
    info <- info[!duplicated(info[,cohorts]), ,drop = FALSE]; rownames(info) <- info[,cohorts]
    info <- droplevels(info)
    
  }
  
  # Add dataMatrix and info to K2 object
  K2data(K2res) <- dataMatrix
  K2info(K2res) <- info
  
  return(K2res)
}