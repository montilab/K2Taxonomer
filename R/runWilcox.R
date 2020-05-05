# Create data.frame of factor variable
.runWilcox <- function(csDF, alternative = "two.sided") {
  
  # Run test
  if(sum(!is.na(csDF$value[csDF$cat])) > 2 & sum(!is.na(csDF$value[!csDF$cat])) > 2) {
    wilcoxOut <- wilcox.test(csDF$value[csDF$cat], 
                             csDF$value[!csDF$cat], 
                             alternative = alternative)
    pval <- wilcoxOut$p.value
    stat <- wilcoxOut$statistic
  } else {
    pval <- NA
    stat <- NA
  }
  
  # Get means of each group and diff
  obsMean <- mean(csDF$value[csDF$cat], na.rm = TRUE)
  altMean <- mean(csDF$value[!csDF$cat], na.rm = TRUE)
  diffMean <- obsMean - altMean
  
  # Get type of test
  test <- c("one-sided Wilcox RS Test", "two-sided Wilcox RS Test")[as.numeric(alternative == "two.sided") + 1]
  
  # Generate results
  data.frame(pval = pval,
             stat = stat,
             df = NA,
             obsMean = obsMean,
             altMean = altMean,
             diffMean = diffMean,
             nhits = NA,
             ncase = NA,
             nalt = NA,
             ndrawn = NA,
             hits = NA,
             test = "1-sided Wilcox Test",
             stringsAsFactors = FALSE)
}