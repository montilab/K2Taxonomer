#' Perform tests to test for associations between modules and meta-variables
#'
#' Adds statistical tests results to the output of K2tax()
#' based on numeric and factor variables in info.
#' @param K2res An object of class K2. The output of K2tax().
#' @return An object of class K2.
#' @keywords clustering
#' @export
#' @examples
#' runTestsMods(K2res)

runTestsMods <- function(K2res) {
    K2results(K2res) <- lapply(K2results(K2res), function(x) {
        
        ## Get module results
        obs <- x$obs
        
        ## Run test on both
        obsRes <- lapply(obs, function(mod) {
            
            ## Run tests for different data types
            do.call(rbind, lapply(names(K2meta(K2res)$infoClass), function(colName) {
                
                ## Get testID
                testID <- K2meta(K2res)$infoClass[colName]
                
                ## Create data.frame of meta variables
                cohort <- K2meta(K2res)$cohort
                if (is.null(cohort)) 
                  cohort <- "sampleID"
                
                if (testID != "survival") {
                  csDF <- K2info(K2res)[, c(cohort, colName)]
                  colnames(csDF) <- c("cohort", "value")
                } else {
                  csDF <- K2info(K2res)[, c(cohort, paste0(colName, "_time"), paste0(colName, 
                    "_status"))]
                  colnames(csDF) <- c("cohort", "time", "value")
                }
                csDF$cat <- csDF$cohort %in% mod
                
                ## Run test
                if (testID == "factor") 
                  out <- .runFisher2sided(csDF)
                if (testID == "factor1") 
                  out <- .runFisher1sided(csDF)
                if (testID == "numeric") 
                  out <- .runWilcox(csDF)
                if (testID == "numeric1") 
                  out <- .runWilcox(csDF, alternative = "greater")
                if (testID == "normal") 
                  out <- .runTtest(csDF)
                if (testID == "normal1") 
                  out <- .runTtest(csDF, alternative = "greater")
                if (testID == "survival") 
                  out <- .runSurvival(csDF)
                
                ## Add value ID
                out <- data.frame(value = colName, out, stringsAsFactors = FALSE)
                return(out)
                
            }))
        })
        
        x$modTests <- obsRes
        return(x)
    })
    
    ## Calculate and merge FDR values
    pValueDF <- data.frame(pval = unlist(lapply(K2results(K2res), function(x) lapply(x$modTests, 
        function(y) y$pval))))
    pValueDF$fdr <- p.adjust(pValueDF$pval, method = "BH")
    pValueDF <- unique(pValueDF)
    K2results(K2res) <- lapply(K2results(K2res), function(x, pValueDF) {
        x$modTests <- lapply(x$modTests, function(y, pValueDF) {
            y <- merge(y, pValueDF)
            if (nrow(y) > 0) {
                y <- y[, c("value", "pval", "fdr", "stat", "df", "obsMean", "altMean", 
                  "diffMean", "nhits", "ncase", "nalt", "ndrawn", "hits", "test")]
            }
            return(y)
        }, pValueDF)
        return(x)
    }, pValueDF)
    
    return(K2res)
}
