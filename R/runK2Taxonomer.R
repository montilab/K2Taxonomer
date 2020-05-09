#' Wrapper function to perform K=2 clustering and annotation
#'
#' This function will generate an object of class, K2.  This will run the K2 Taxonomer procedure, differential analysis, and finally hyperenrichment on a named list of feature sets.
#' @param eSet An expression set object.
#' @param cohorts The column in phenotype data of eSet that has cohort ID's. Default NULL if no pre-processing of data.
#' @param vehicle The value in the cohort variable that contains the vehicle ID. Default NULL if no vehicle to be used.
#' @param covariates Covariates in phenotype data of eSet to control for in differential analysis.
#' @param block Block parameter in limma for modelling random-like effects.
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
#' @param oneoff Logical. Allow 1 member clusters?
#' @param stabThresh Threshold for ending clustering.
#' @param geneURL Optional. Named list of URLs to gene information.
#' @param genesetURL Optional. Named list of URLs to geneset information.
#' @keywords clustering
#' @export
#' @import limma
#' @import Biobase
#' @import GSVA
#' @examples
#' runK2Taxonomer(eSet, cohorts = NULL, vehicle = NULL, nFeats = nrow(eSet)*0.02, nBoots = 200, clustFunction = hclust_wrapper_fast, info = NULL, genesets = NULL)

runK2Taxonomer <- function(eSet,
                         cohorts = NULL,
                         vehicle = NULL,
                         covariates = NULL,
                         block = NULL,
                         use = c("Z", "MEAN"),
                         nFeats = nrow(eSet)*0.02,
                         featMetric = c("square", "mad"),
                         nBoots = 200,
                         clustFunc = hclust_wrapper,
                         linkage = "mcquitty",
                         info = NULL,
                         infoClass = NULL,
                         genesets = NULL,
                         qthresh = 0.05,
                         cthresh = 0,
                         ntotal = 20000,
                         ssGSEAalg = c("ssgsea", "gsva", "zscore", "plage"),
                         ssGSEAminQvalue = 0.05,
                         ssGSEAcores = 0, 
                         oneoff = FALSE,
                         stabThresh = -1,
                         geneURL = NULL,
                         genesetURL = NULL
                         ){

  # Match arguments
  K2res <- K2preproc(eSet,
                     cohorts = cohorts,
                     vehicle = vehicle,
                     covariates = covariates,
                     block = block,
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

  # Cluster data and create S4 object
  cat("Running K2Taxonomer bootstraps.\n")
  K2res <- K2tax(K2res)

  # Check clusters for significance of specified variables
  cat("Checking clusters for significant associations with meta-variables.\n")
  if (!is.null(infoClass)) {
    K2res <- runTests_mods(K2res)
  }

  # Add differential analysis results
  cat("Running differential analysis.\n")
  K2res <- runDGE_mods(K2res)

  if (!is.null(genesets)) {
    
    # Add genesets and pathway maps to analysis
    K2genesets(K2res) <- genesets # Add gene sets
    K2gene2Pathway(K2res) <- getGenePathways(genesets) # Add gene set mapping
    
    # Perform hyperenrichment analysis
    cat("Running hyperenrichment analysis on genes with FDR<", qthresh, " & |coefficient|>", cthresh, ".\n", sep = "")
    K2res <- runGSE_mods(K2res)
    
    # Perform ssGSEA or ssGVA
    cat("Running ssGSEA.\n")
    K2res <- runGSVA_mods(K2res)
  }
  
  # Add URLs
  if (!is.null(geneURL)) {
    K2geneURL(K2res) <- geneURL
  }
  if (!is.null(genesetURL)) {
    K2genesetURL(K2res) <- genesetURL
  }

  return(K2res)

}
