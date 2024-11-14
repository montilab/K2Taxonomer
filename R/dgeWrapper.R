.dgeWrapper <- function(K2res, keep, Fstat=FALSE) {
  
    pD <- K2colData(K2res)[keep,]
    cohorts <- K2meta(K2res)$cohorts
    vehicle <- K2meta(K2res)$vehicle
    variables <- K2meta(K2res)$variables
    logCounts <- K2meta(K2res)$logCounts
    use <- K2meta(K2res)$use
    
    ## Get unique cohorts and remove vehicle
    gUnique <- as.character(unique(pD[, cohorts]))
    
    condVehicle <- FALSE
    if (!is.null(vehicle)) {
      if(vehicle %in% gUnique) {
        gUnique <- gUnique[gUnique != vehicle]
        condVehicle <- TRUE
      }
    }
    
    ## Create new variable
    pD$GROUP <- factor(pD[, cohorts], levels=c(vehicle, gUnique))
    
    ## Drop levels
    pD <- droplevels(pD)
    
    ## Create formula
    design <- model.matrix(as.formula(paste0("~", "GROUP",
                                             .formatCov(variables))), data=pD)
    
    ## Check that model is full rank with variables, if not model w/o covariate
    if(!is.fullrank(design)) {
        design <- model.matrix(as.formula(paste0("~", "GROUP",
                                                 .formatCov(NULL))), data=pD) 
    }
    colnames(design) <- sub("GROUP", "X", colnames(design))
    
    ## Fit model
    fit <- eBayes(lmFit(K2eMat(K2res)[, rownames(pD)], design), 
                  trend=logCounts, robust=logCounts)
    
    ## Get gene-level results
    gUniqueX <- paste0("X", gUnique)
    modS <- topTable(fit, number=Inf, sort.by="none",
                         coef=gUniqueX[gUniqueX %in% colnames(design)])
    
    ## Get Fstat
    if (Fstat) {
        Fvec <- modS$F
        names(Fvec) <- rownames(modS)
    }
    
    ## Get model log fold-changes (means if no intercept)
    modS <- as.matrix(modS[, colnames(modS) %in% gUniqueX])
    
    ## If use == Z get the test statistics instead
    if (use == "Z") {
        modZ <- do.call(cbind, lapply(colnames(modS), function(g, fit) {
            topTable(fit, number=Inf, sort.by="none", coef=g)$t
        }, fit))
        colnames(modZ) <- colnames(modS)
        rownames(modZ) <- rownames(modS)
        modS <- modZ
    }
    colnames(modS) <- sub("X", "", colnames(modS))
    
    # Remove fit object
    rm(fit)
    
    ## Add 0's
    modS <- cbind(modS, 0)
    if (condVehicle) {
        colnames(modS)[ncol(modS)] <- vehicle
    } else {
        colnames(modS)[ncol(modS)] <- gUnique[!gUnique %in% colnames(modS)]
    }
    
    ## Return
    if (Fstat) {
        return(list(dataMatrix=modS, Fvec=Fvec))
    } else {
        return(modS)
    }
}
