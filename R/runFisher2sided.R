# Create data.frame of factor variable
.runFisher2sided <- function(csDF) {
  
  # Run test
  if(sum(!is.na(csDF$value[csDF$cat])) > 0 & sum(!is.na(csDF$value[!csDF$cat])) > 0) {
    fisherOut <- fisher.test(csDF$cat, csDF$value, simulate.p.value = TRUE)
    pval <- fisherOut$p.value
    stat <- fisherOut$estimate
  } else {
    pval <- NA
    stat <- NA
  }
  
  if(is.null(stat)) {
    stat <- NA
  }
  
  # Generate results
  data.frame(pval = pval,
             stat = stat,
             df = NA,
             obsMean = NA,
             altMean = NA,
             diffMean = NA,
             nhits = NA,
             ncase = NA,
             nalt = NA,
             ndrawn = NA,
             hits = NA,
             test = "two-sided Fishers Exact Test",
             stringsAsFactors = FALSE)
}