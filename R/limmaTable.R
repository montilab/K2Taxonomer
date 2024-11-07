.limmaTable <- function(x, f, l, cv = NULL, d = NULL, v = NULL) {
  
  if(is.null(cv) | !is.null(v)) {
    eF <- eBayes(f, trend = l, robust=l)
    contFit <- suppressWarnings(topTable(eF,
                                         coef = x, number=Inf, sort.by="none"))
    
    if(!is.null(v)) {
      contFitV <- contFit
    }
    
  }
    
  if(!is.null(cv) | !is.null(v)) {
      
    fC <- suppressWarnings(contrasts.fit(f, makeContrasts(contrasts = cv[x], 
                                                          levels = d)))
    eF <- eBayes(fC, trend = l, robust=l)
    contFit <- suppressWarnings(topTable(eF,
                                         number=Inf, sort.by="none"))
  }
  contFit$df <- eF$df.total
  contFit <- contFit[, c("logFC", "AveExpr",
                         "t", "df", "P.Value")]
  colnames(contFit) <- c("coef", "mean", "t", "df", "pval")
  contFit$feat <- rownames(contFit)
  contFit$fdr <- NA
  contFit$edge <- sub("mods", "", x)
  rownames(contFit) <- NULL
  contFit <- contFit[, c("feat", "coef", "mean", "t", "df", "pval", "fdr", "edge")]
  
  if(!is.null(v)) {
    contFit$pvalV <- contFitV$P.Value
    contFit$coefV <- contFitV$logFC
  }
  
  return(contFit)
}
