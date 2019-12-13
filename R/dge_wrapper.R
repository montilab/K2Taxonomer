.dge_wrapper <- function(eSet, cohorts, vehicle = NULL, covariates = NULL, use = c("Z", "MEAN")){
  
  # Get unique cohorts and remove vehicle
  gUnique <- as.character(unique(pData(eSet)[,cohorts]))
  
  if(!is.null(vehicle)){
    gUnique <- gUnique[gUnique != vehicle]
  }
  
  # Create new variable
  pData(eSet)$GROUP <- factor(pData(eSet)[,cohorts], levels = c(vehicle, gUnique))
  
  # Drop levels
  pData(eSet) <- droplevels(pData(eSet))
  
  # Create formula
  if(!is.null(vehicle)){
    design <- model.matrix(as.formula(paste0("~", "GROUP", .formatCov(covariates))), data = pData(eSet))
  } else {
    design <- model.matrix(as.formula(paste0("~ 0+", "GROUP", .formatCov(covariates))), data = pData(eSet))
  }
  colnames(design) <- sub("GROUP", "X", colnames(design))
  
  # Fit model
  fit <- eBayes(lmFit(eSet, design))
  
  # Get model log fold-changes (means if no intercept)
  modStats <- topTable(fit, number = Inf, sort.by = "none")
  modStats <- modStats[,colnames(modStats) %in% paste0("X", gUnique)]
  
  # If use == Z get the test statistics instead
  if(use == "Z"){
    modZ <- do.call(cbind, lapply(colnames(modStats), function(g, fit){
      topTable(fit, number = Inf, sort.by = "none", coef = g)$t
    }, fit))
    colnames(modZ) <- colnames(modStats); rownames(modZ) <- rownames(modStats)
    modStats <- modZ
  }
  colnames(modStats) <- sub("X", "", colnames(modStats))
  
  # Add vehicle 0's
  if(!is.null(vehicle)){
    modStats <- cbind(modStats, 0)
    colnames(modStats)[ncol(modStats)] <- "Vehicle"
  } else {
    modStats <- modStats - rowMeans(modStats)
  }
  
  # Return
  return(modStats)
}