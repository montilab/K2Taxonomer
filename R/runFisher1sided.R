## Create data.frame of factor variable
.runFisher1sided <- function(csDF) {

    if (sum(!is.na(csDF$value[csDF$cat])) > 0 &
        sum(!is.na(csDF$value[!csDF$cat])) > 0) {
        fisherOut <- fisher.test(csDF$cat, csDF$value, alternative="greater")
        pval <- fisherOut$p.value
        stat <- fisherOut$estimate
    } else {
        pval <- NA
        stat <- NA
    }

    if (is.null(stat)) {
        stat <- NA
    }

    ## Get hits
    hits <- paste(csDF$cohort[csDF$cat & csDF$value == levels(csDF$value)[2] &
        !is.na(csDF$value)], collapse=", ")

    ## Get contingency table
    cTab <- table(csDF$cat, csDF$value)

    ## Get some stats
    nhits <- cTab[2, 2]
    ncase <- sum(cTab[, 2])
    nalt <- sum(cTab[, 1])
    ndrawn <- sum(cTab[2, ])

    ## Generate results
    data.frame(pval=pval, stat=stat, df=NA, obsMean=NA,
        altMean=NA, diffMean=NA, nhits=nhits, ncase=ncase,
        nalt=nalt, ndrawn=ndrawn, hits=hits,
        test="one-sided Fishers Exact Test", stringsAsFactors=FALSE)

}
