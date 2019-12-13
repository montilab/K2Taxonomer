# Create data.frame of factor variable
.runTtest <- function(csDF, alternative = "two.sided") {
  
  # Run test
  TtestOut <- t.test(csDF$value[csDF$cat], 
                     csDF$value[!csDF$cat], 
                     alternative = alternative)
  pval <- TtestOut$p.value
  stat <- TtestOut$statistic
  df <- TtestOut$parameter
  
  # Get means of each group and diff
  obsMean <- mean(csDF$value[csDF$cat], na.rm = TRUE)
  altMean <- mean(csDF$value[!csDF$cat], na.rm = TRUE)
  diffMean <- obsMean - altMean
  
  # Generate results
  data.frame(pval = pval,
             stat = stat,
             df = df,
             obsMean = obsMean,
             altMean = altMean,
             diffMean = diffMean,
             nhits = NA,
             ncase = NA,
             nalt = NA,
             ndrawn = NA,
             hits = NA,
             stringsAsFactors = FALSE)
}