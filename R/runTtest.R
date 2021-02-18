## Create data.frame of factor variable
.runTtest <- function(csDF, alternative="two.sided") {

    ## Run test
    if (sum(!is.na(csDF$value[csDF$cat])) > 2 &
        sum(!is.na(csDF$value[!csDF$cat])) > 2) {
        TtestOut <- t.test(csDF$value[csDF$cat], csDF$value[!csDF$cat],
            alternative=alternative)
        pval <- TtestOut$p.value
        stat <- TtestOut$statistic
        df <- TtestOut$parameter
    } else {
        pval <- NA
        stat <- NA
        df <- NA
    }

    ## Get means of each group and diff
    obsMean <- mean(csDF$value[csDF$cat], na.rm=TRUE)
    altMean <- mean(csDF$value[!csDF$cat], na.rm=TRUE)
    diffMean <- obsMean - altMean

    ## Get type of test
    test <- c("one-sided T-test", "two-sided T-test")[as.numeric(alternative ==
        "two.sided") + 1]

    ## Generate results
    data.frame(pval=pval, stat=stat, df=df, obsMean=obsMean,
        altMean=altMean, diffMean=diffMean, nhits=NA, ncase=NA,
        nalt=NA, ndrawn=NA, hits=NA, test=test, stringsAsFactors=FALSE)
}
