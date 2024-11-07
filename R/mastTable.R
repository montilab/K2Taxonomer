.mastTable <- function(x, f, em, cv = NULL, ms = NULL, d = NULL, v = NULL) {
  
  if(is.null(cv) | !is.null(v)) {
    
    # Get hurdle test result
    lrRes <- suppressMessages(lrTest(f, CoefficientHypothesis(x), parallel = TRUE))
    
    # Get fold change
    cs <- colnames(f@coefD)
    c1 <- rep(0, length(cs))
    names(c1) <- cs
    c1[1] <- c1[x] <- 1
    fcRes <- logFC(f, contrast1 = c1)[[1]]
    
    if(!is.null(v)) {
      lrResV <- lrRes
      fcResV <- fcRes
    }
    
  }
  
  if(!is.null(cv) | !is.null(v)) {
    
    # Fit contrasts
    cf <- suppressWarnings(makeContrasts(contrasts = cv[x], levels = d))
    lrRes <- suppressMessages(lrTest(f, cf, parallel = TRUE))
    
    # Get log fold changes
    
    # Create contrats
    cs <- colnames(f@coefD)
    c1 <- c0 <- rep(0, length(cs))
    names(c1) <- names(c0) <- cs
    if(is.null(v)) {
      msOut <- ms[ms != x]
    } else {
      msOut <- ms[!ms %in% c(x, v)]
    }
    c1[x] <- 1
    c0[msOut] <- 1/length(msOut)
    
    # Calculate fold changes
    fcRes <- logFC(f, c0, c1)[[1]]
    
  }
  
  # Compile results
  contFit <- data.frame(feat = rownames(fcRes), 
                        coef = fcRes[, 1],
                        mean = em[rownames(fcRes)],
                        x2 = lrRes[,,1][,3],
                        df = lrRes[,,2][,3],
                        pval = lrRes[,,3][,3])
  contFit$fdr <- NA
  contFit$edge <- sub("^mods", "", x)
  
  if(!is.null(v)) {
    contFit$pvalV <- lrResV[,,3][,3]
    contFit$coefV <- fcResV[, 1]
  }
  
  contFit <- contFit[!is.na(contFit$coef),]
  
  return(contFit)
  
}
