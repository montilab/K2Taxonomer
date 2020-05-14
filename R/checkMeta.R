## Function to change meta data if new value is entered.
.checkMeta <- function(K2res, arg = NULL, argValue = NULL) {
    
    if (!is.null(arg)) {
        
        if (is.null(argValue)) {
            argValue <- K2meta(K2res)[[arg]]
        }
        if (is.null(argValue)) {
            stop(cat("No value of ", arg, " specified.\n"))
        }
        
    } else {
        
        ## Run Stopping criteria
        K2m <- K2meta(K2res)
        
        ## ssGSEAalg
        if (!K2m$ssGSEAalg %in% c("ssgsea", "gsva", "zscore", "plage")) {
            stop("Argument, ssGSEAalg, must be one of: 'ssgsea', 'gsva',
                'zscore', 'plage'.\n")
        }
        
        ## linkage
        if (!K2m$linkage %in% c("mcquitty", "ward.D", "ward.D2", "single", "complete", 
            "average", "centroid")) {
            stop("Argument, linkage, must be one of 'mcquitty', 'ward.D',
                'ward.D2', 'single', 'complete', 'average', 'centroid'.\n")
        }
        
        ## use
        if (!K2m$use %in% c("Z", "MEAN")) {
            stop("Argument, linkage, must be one of 'mcquitty', 'ward.D',
                'ward.D2', 'single', 'complete', 'average', 'centroid'.\n")
        }
        
        ## nFeats
        if (!K2m$nFeats <= 1) {
            stop("Argument, nFeats, must be one a positive value > 1.\n")
        }
        
        ## cohorts
        if (!is.null(K2m$cohorts) && !K2m$cohorts %in% colnames(pData(K2eSet(K2res)))) {
            stop("Argument, cohorts, must match column name in phenodata of
                expression set.\n")
        }
        
        ## vehicle
        if (!is.null(K2m$vehicle) && !K2m$vehicle %in% colnames(pData(K2eSet(K2res)))) {
            stop("Argument, vehicle, must match column name in phenoData of
                expressionSet.\n")
        }
        
        ## block
        if (!is.null(K2m$block) && !K2m$block %in% colnames(pData(K2eSet(K2res)))) {
            stop("Argument, block, must match column name in 
                phenoData of expressionSet.\n")
        }
        
        ## recalcDataMatrix = FALSE,
        if (K2m$recalcDataMatrix) {
            stop("Argument, recalcDataMatrix, must be logical.\n")
        }
        
        ## nBoots
        if (K2m$nboots <= 1) {
            stop("Argument, nboots, must be a value > 1.\n")
        }
        
        ## clustFunc
        if (class(K2m$clustFunc) != "function") {
            stop("Argument, clustFunc, must be of class, function.\n")
        }
        
        ## clustCors
        if (K2m$clustCors <= 1) {
            stop("Argument, clustCors, must be a value > 1.\n")
        }
        
        ## clustList
        if (class(K2m$clustList) != "list") {
            stop("Argument, clustList, must of class, list.\n")
        }
        
        ## info
        if (!is.null(K2m$info) && class(K2m$info) != "data.frame") {
            stop("Argument, , .\n")
        }
        
    }
    
    return(argValue)
    
}
