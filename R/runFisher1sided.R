# Create data.frame of factor variable
.runFisher1sided <- function(csDF) {
  
  fisherOut <- fisher.test(csDF$cat, csDF$value, alternative = "greater")
  pval <- fisherOut$p.value
  
  # Get estimate if possible
  stat <- fisherOut$estimate
  
  # Get hits
  hits <- paste(
    csDF$cohort[csDF$cat & csDF$value == levels(csDF$value)[2] & !is.na(csDF$value)],
    collapse = ", "
  )
  
  # Get contingency table
  cTab <- table(csDF$cat, csDF$value)
  
  # Get some stats
  nhits <- cTab[2,2]
  ncase <- sum(cTab[,2])
  nalt <- sum(cTab[,1])
  ndrawn <- sum(cTab[2,])
  nleft <- sum(cTab[1,])
  
  # Generate results
  data.frame(pval = pval,
             stat = stat,
             df = NA,
             obsMean = NA,
             altMean = NA,
             diffMean = NA,
             nhits = nhits,
             ncase = ncase,
             nalt = nalt,
             ndrawn = ndrawn,
             hits = hits,
             test = "1-sided Fisher Exact Test",
             stringsAsFactors = FALSE)

}