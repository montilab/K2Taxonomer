## Function to change meta data if new value is entered.
.checkK2 <- function(K2res, arg=NULL, argValue=NULL, inputsOnly=FALSE) {

    if (!is.null(arg)) {

        if (is.null(argValue)) {
            argValue <- K2meta(K2res)[[arg]]
        }
        if (is.null(argValue)) {
            stop(cat("No value of ", arg, " specified.\n"))
        }

    } else {

        argValue <- NULL

        ## Run stopping criteria for meta slot
        K2m <- K2meta(K2res)

        ## cohorts

        ## Check if column name is found
        if (!is.null(K2m$cohorts) && !K2m$cohorts %in%
            colnames(pData(K2eSet(K2res)))) {
            stop("Argument, cohorts, must match column name in phenodata of
                Expression set.\n")
        }

        ## Check for not-allowed characters in cohort factor levels
        if (!is.null(K2m$cohorts) && sum(grepl(" |-",
            as.character(pData(K2eSet(K2res))[,
            K2m$cohorts]))) > 0) {
            stop("Values in cohorts variable cannot contain spaces or '-'.\n")
        }

        ## Check for factor or character
        if (!is.null(K2m$cohorts) && !(is(pData(K2eSet(K2res))[,
            K2m$cohorts], "factor") | is(pData(K2eSet(K2res))[,
            K2m$cohorts], "character"))) {
            stop("Values in cohorts variable must be character or factor.\n")
        }

        ## Vehicle
        if (!is.null(K2m$vehicle) && !K2m$vehicle %in% pData(K2eSet(K2res))[,
            K2m$cohorts]) {
            stop("Argument, vehicle, not found in cohort variable.\n")
        }

        ## Vehicle and recalculating data matrix
        if (!is.null(K2m$vehicle) && (K2m$featMetric == "F" |
            K2m$recalcDataMatrix)) {
            stop("Specifying argument, vehicle, currently not supported if
                argument, featMetric='F' or if argument,
                recalcDataMatrix=TRUE.\n")
        }

        ## covariates
        if (!is.null(K2m$covariates) && !K2m$covariates %in%
            colnames(pData(K2eSet(K2res)))) {
            stop("Argument, covariates, must match column name in phenoData of
                ExpressionSet.\n")
        }

        # ## Covariated and recalculating data matrix
        # if (!is.null(K2m$covariates) && (K2m$featMetric == "F" |
        #     K2m$recalcDataMatrix)) {
        #     stop("Specifying argument, covariates, currently not supported if
        #         argument, featMetric='F' or if argument,
        #         recalcDataMatrix=TRUE.\n")
        # }

        ## block
        if (!is.null(K2m$block) && !K2m$block %in%
            colnames(pData(K2eSet(K2res)))) {
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
            stop("Argument, nFeats, must be either 'sqrt' or a positive value >
                1.\n")
        }

        ## featMetric
        if (!K2m$featMetric %in% c("sd", "mad", "Sn", "Qn", "F",
            "square")) {
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
        if (!K2m$linkage %in% c("mcquitty", "ward.D", "ward.D2",
            "single", "complete", "average", "centroid")) {
            stop("Argument, linkage, must be one of 'mcquitty', 'ward.D',
                'ward.D2', 'single', 'complete', 'average', 'centroid'.\n")
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
        if (!K2m$ssGSEAalg %in% c("ssgsea", "gsva", "zscore",
            "plage")) {
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

        ## Run stopping criteria for eSet slot

        ## Missing values
        if (sum(is.na(exprs(K2eSet(K2res)))) > 0) {
            stop("assayData slot in ExpressionSet object contains missing
                values.\n")
        }

        ## Run stopping criteria for geneURL split

        ## Missing names
        if (length(K2geneURL(K2res)) > 0 &&
            is.null(names(K2geneURL(K2res)))) {
            stop("Argument, geneURL, must be a named vector.\n")
        }

        ## Run stopping criteria for genesetURL split

        ## Missing names
        if (length(K2genesetURL(K2res)) > 0 &&
            is.null(names(K2genesetURL(K2res)))) {
            stop("Argument, genesetURL, must be a named vector.\n")
        }

        ## Run stopping criteria for genesets slot

        ## Gene names not found in expression set
        if (length(K2genesets(K2res)) > 0 &&
            sum(unique(unlist(K2genesets(K2res))) %in%
            rownames(K2eSet(K2res))) == 0) {
            stop("No features in argument, genesets, found in ExpressionSet.\n")
        }

        ## Geneset names not found
        if (length(K2genesets(K2res)) > 0 &&
            is.null(names(K2genesets(K2res)))) {
            stop("Argument, genesets, needs to be a named list.\n")
        }

        ## Geneset names contain ';'
        if (length(K2genesets(K2res)) > 0 &&
            sum(grepl(";", names(K2genesets(K2res)))) >
            0) {
            stop("Names in argument, genesets, cannot contain ';'.\n")
        }

        if (!inputsOnly) {

            ## Run stopping criteria for info slot

            ## Mismatch info and dataMatrix
            if (nrow(K2info(K2res)) > 0 && nrow(K2info(K2res)) !=
                ncol(K2data(K2res))) {
                stop("No. columns of info slot doesn't equal No. columns in
                    dataMatrix slot.\n")
            }

            ## infoClass names
            if (!is.null(K2m$infoClass) && mean(names(K2m$infoClass) %in%
                colnames(K2info(K2res))) != 1) {
                stop("Names of argument, infoClass, don't match column names of
            argument, info or phenoData of expressionSet.\n")
            }

            ## infoClass and cohorts
            if (!is.null(K2m$infoClass) && !is.null(K2m$cohorts)) {
                stop("Specifying argument, infoClass, currently not supported
                    if argument, cohorts, is specified.\n")
            }

            ## Run stopping criteria for dataMatrix slot

            ## Missing values
            if (sum(is.na(K2data(K2res))) > 0) {
                stop("dataMatrix slot contains missing values.\n")
            }

        }

    }

    return(argValue)

}
