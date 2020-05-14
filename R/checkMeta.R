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

        argValue <- NULL

        ## Run Stopping criteria
        K2m <- K2meta(K2res)

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

        ## covariates
        if (!is.null(K2m$covariates) && !K2m$covariates %in% colnames(pData(K2eSet(K2res)))) {
            stop("Argument, covariates, must match column name in phenoData of
                expressionSet.\n")
        }

        ## block
        if (!is.null(K2m$block) && !K2m$block %in% colnames(pData(K2eSet(K2res)))) {
            stop("Argument, block, must match column name in
                phenoData of expressionSet.\n")
        }

        ## logCounts
        if (!is.logical(K2m$logCounts)) {
            stop("Argument, logCounts, must be logical.\n")
        }

        ## use
        if (!K2m$use %in% c("Z", "MEAN")) {
            stop("Argument, linkage, must be one of 'mcquitty', 'ward.D',
                'ward.D2', 'single', 'complete', 'average', 'centroid'.\n")
        }

        ## nFeats
        if (K2m$nFeats <= 1) {
            stop("Argument, nFeats, must be either 'sqrt' or a positive value > 1.\n")
        }

        ## featMetric
        if (!K2m$featMetric %in% c("sd", "mad", "Sn", "Qn", "F", "square")) {
            stop("Argument, featMetric, must be one of 'sd', 'mad', 'Sn', 'Qn',
            'F', 'square'.\n")
        }

        ## recalcDataMatrix
        if (!is.logical(K2m$recalcDataMatrix)) {
            stop("Argument, recalcDataMatrix, must be logical.\n")
        }

        ## nBoots
        if (K2m$nBoots <= 1) {
            stop("Argument, nBoots, must be a value > 1.\n")
        }

        ## clustFunc
        if (!is(K2m$clustFunc, "function")) {
            stop("Argument, clustFunc, must be of class, function.\n")
        }

        ## clustCors
        if (K2m$clustCors < 1) {
            stop("Argument, clustCors, must be a value >= 1.\n")
        }

        ## clustList
        if (!is(K2m$clustList, "list")) {
            stop("Argument, clustList, must of class, list.\n")
        }

        ## linkage
        if (!K2m$linkage %in% c("mcquitty", "ward.D", "ward.D2", "single", "complete",
            "average", "centroid")) {
            stop("Argument, linkage, must be one of 'mcquitty', 'ward.D',
                'ward.D2', 'single', 'complete', 'average', 'centroid'.\n")
        }

        ## info
        if (!is.null(K2m$info) && !is(K2m$info, "data.frame")) {
            stop("Argument, info, must be data frame.\n")
        }

        ## infoClass
        if (!is.null(K2m$infoClass) && mean(names(K2m$infoClass) %in% colnames(K2info(K2res))) != 1) {
            stop("Names of argument, infoClass, don't match column names of
            argument, info, or phenoData of expressionSet.\n")
        }

        ## qthresh
        if (K2m$qthresh > 1 | K2m$qthresh < 0) {
            stop("Argument, qthresh, must be a value between 0 and 1.\n")
        }

        ## cthresh
        if (K2m$cthresh < 0) {
            stop("Argument, cthresh, must be a value greater than 0.\n")
        }

        ## ntotal
        if (K2m$ntotal < 1) {
            stop("Argument, ntotal, must be a value greater than 1.\n")
        }

        ## ssGSEAalg
        if (!K2m$ssGSEAalg %in% c("ssgsea", "gsva", "zscore", "plage")) {
            stop("Argument, ssGSEAalg, must be one of: 'gsva', 'ssgsea',
                'zscore', 'plage'.\n")
        }

        ## ssGSEAcores
        if (K2m$ssGSEAcores < 1) {
            stop("Argument, ssGSEAcores, must be an integer greater than 0.\n")
        }

        ## oneoff
        if (!is.logical(K2m$oneoff)) {
            stop("Argument, oneoff, must be logical.\n")
        }

        ## stabThresh
        if (K2m$stabThresh > 1 | K2m$stabThresh < 0) {
            stop("Argument, stabThresh, must be a value between 0 and 1.\n")
        }

    }

    return(argValue)

}
