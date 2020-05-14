## Create data.frame of factor variable
.runSurvival <- function(csDF) {
    
    ## Run test
    SurvOut <- survdiff(Surv(csDF$time, csDF$value) ~ csDF$cat)
    pval <- 1 - pchisq(SurvOut$chisq, length(SurvOut$n) - 1)
    stat <- SurvOut$chisq
    
    ## Generate results
    data.frame(pval = pval, stat = stat, df = NA, obsMean = NA, altMean = NA, diffMean = NA, 
        nhits = NA, ncase = NA, nalt = NA, ndrawn = NA, hits = NA, test = "Survival", 
        stringsAsFactors = FALSE)
}
