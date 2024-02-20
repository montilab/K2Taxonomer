.dgeWrapper <- function(eSet, cohorts, vehicle=NULL, covariates=NULL, 
                        use=c("Z", "MEAN"), Fstat=FALSE, logCounts=FALSE) {
    
    ## Get unique cohorts and remove vehicle
    gUnique <- as.character(unique(pData(eSet)[, cohorts]))
    
    if (!is.null(vehicle)) {
        gUnique <- gUnique[gUnique != vehicle]
    }
    
    ## Create new variable
    pData(eSet)$GROUP <- factor(pData(eSet)[, cohorts], levels=c(vehicle,
                                                                 gUnique))
    
    ## Drop levels
    pData(eSet) <- droplevels(pData(eSet))
    
    ## Create formula
    design <- model.matrix(as.formula(paste0("~", "GROUP",
                                             .formatCov(covariates))), data=pData(eSet))
    
    ## Check that model is full rank with covariates, if not model w/o covariates
    if(!is.fullrank(design)) {
        design <- model.matrix(as.formula(paste0("~", "GROUP",
                                                 .formatCov(NULL))), data=pData(eSet)) 
    }
    colnames(design) <- sub("GROUP", "X", colnames(design))
    
    ## Fit model
    fit <- eBayes(lmFit(eSet, design), trend=logCounts, robust=logCounts)
    
    ## Get gene-level results
    gUniqueX <- paste0("X", gUnique)
    modStats <- topTable(fit, number=Inf, sort.by="none",
                         coef=gUniqueX[gUniqueX %in% colnames(design)])
    
    ## Get Fstat
    if (Fstat) {
        Fvec <- modStats$F
        names(Fvec) <- rownames(modStats)
    }
    
    ## Get model log fold-changes (means if no intercept)
    modStats <- as.matrix(modStats[, colnames(modStats) %in%
                                       gUniqueX])
    
    ## If use == Z get the test statistics instead
    if (use == "Z") {
        modZ <- do.call(cbind, lapply(colnames(modStats), function(g,
                                                                   fit) {
            topTable(fit, number=Inf, sort.by="none", coef=g)$t
        }, fit))
        colnames(modZ) <- colnames(modStats)
        rownames(modZ) <- rownames(modStats)
        modStats <- modZ
    }
    colnames(modStats) <- sub("X", "", colnames(modStats))
    
    # Remove fit object
    rm(fit)
    
    ## Add 0's
    modStats <- cbind(modStats, 0)
    if (!is.null(vehicle)) {
        colnames(modStats)[ncol(modStats)] <- vehicle
    } else {
        colnames(modStats)[ncol(modStats)] <- gUnique[!gUnique %in% colnames(modStats)]
    }
    
    ## Return
    if (Fstat) {
        return(list(dataMatrix=modStats, Fvec=Fvec))
    } else {
        return(modStats)
    }
}
