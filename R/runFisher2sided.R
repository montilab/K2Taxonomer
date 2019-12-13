# Create data.frame of factor variable
.runFisher2sided <- function(csDF) {
  
  # Run test
  fisherOut <- fisher.test(csDF$cat, csDF$value, simulate.p.value = TRUE)
  pval <- fisherOut$p.value
  
  # Get estimate if possible
  if (!is.null(fisherOut$estimate)) {
    stat <- fisherOut$estimate
  } else {
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
             test = "2-sided Fisher Exact Test",
             stringsAsFactors = FALSE)
}