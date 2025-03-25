.signatureWrapper <- function(K2r, mods, GENE = TRUE) {
  
    # Get objects from K2r object
    pD <- K2colData(K2r)
    cohorts <- K2meta(K2r)$cohorts
    vehicle <- K2meta(K2r)$vehicle
    variables <- K2meta(K2r)$variables
    block <- K2meta(K2r)$block
    logCounts <- K2meta(K2r)$logCounts
    DGEexpThreshold <- K2meta(K2r)$DGEexpThreshold

    # Get function for expression data matrix
    if(GENE) {
      if(nrow(K2eMatDS(K2r)) != 0) {
        EXPFUNC <- K2eMatDS
      } else {
        EXPFUNC <- K2eMat
      }
    } else {
      EXPFUNC <- K2gMat
    }
    
    # IF !GENE don't assume log counts
    if(!GENE) {
      logCounts <- FALSE
    }
    
    ## Remove vehicle from mods and make a data frame
    if (!is.null(vehicle)) {
        mods <- mods[names(mods) != vehicle]
    }
    modDF <- data.frame(mods=as.character(mods), GROUP=names(mods),
        stringsAsFactors=FALSE)
    
    ## Get min module size
    minModSize <- min(table(mods))

    modS <- NULL
    desForm <- NULL

    ## Run if there are subgroups to compare and at least three
    ## observations
    if (length(unique(modDF$mods)) > 1 && (minModSize > 2 | !is.null(cohorts))) {

        ## If replicates in data get unique cohorts
        if (is.null(cohorts)) {
            cohorts <- "GROUP"
            pD[, cohorts] <- colnames(EXPFUNC(K2r))
        }

        ## Get unique groups
        gUnique <- unique(pD[, cohorts])

        if (!is.null(vehicle)) {
            gUnique <- gUnique[gUnique != vehicle]
        }

        ## Subset data for mods
        eSub <- EXPFUNC(K2r)[, pD[, cohorts] %in% c(vehicle,
            modDF$GROUP)]

        ## Drop levels
        pD <- droplevels(pD)

        ## Create new variable in pData by merging with mods
        ## data.frame
        pD$GROUP <- factor(pD[, cohorts], levels=c(vehicle,
            modDF$GROUP))
        pD$Rownames <- rownames(pD)
        pD <- merge(pD, modDF, all.x=TRUE)
        rownames(pD) <- pD$Rownames
        pD <- pD[colnames(eSub), , drop=FALSE]
        
        # Get minimum number of observations per module
        minObs <- min(table(pD$mods))
        if(minObs > 2) {
          
          ## Add vehicle mod
          pD$mods[is.na(pD$mods)] <- "0"
          pD$mods <- as.factor(pD$mods)
          pDmods <- paste0("mods", as.character(sort(unique(pD$mods))))
          
          ## If block add to pD
          if (!is.null(block)) {
            pD$block <- as.character(pD[, block])
          }
          
          ## Remove lowly expressed genes
          if(DGEexpThreshold > 0 & GENE) {
            eSub <- .K2filterGenes(eSub, pD$mods, DGEexpThreshold)
          }
  
          ## Create design matrix
          if(!(length(pDmods) > 2 & is.null(vehicle))) {
          
            desForm <- paste0("~mods", .formatCov(variables))
            design <- model.matrix(as.formula(desForm), data=pD)
            
            ## Check that model is full rank with variables, if not model w/o variables
            if(!is.fullrank(design)) {
                desForm <- paste0("~mods", .formatCov(NULL))
                design <- model.matrix(as.formula(desForm), data=pD)
            }
            
            conVec <- sapply(pDmods[-1], function(x) {
              paste0(x, "-(", paste(pDmods[pDmods !=x][-1], collapse="+"), ")/",
                     sum(pDmods !=x) - 1)
            })
            
          } else {
            
            desForm <- paste0("~0+mods", .formatCov(variables))
            design <- model.matrix(as.formula(desForm), data=pD)
            
            ## Check that model is full rank with variables, if not model w/o variables
            if(!is.fullrank(design)) {
              desForm <- paste0("~0+mods", .formatCov(NULL))
              design <- model.matrix(as.formula(desForm), data=pD)
            }
            
            conVec <- sapply(pDmods, function(x) {
              paste0(x, "-(", paste(pDmods[pDmods !=x], collapse="+"), ")/",
                     sum(pDmods !=x))
            })
            
          }
          
          ## Run duplicateCorrelation
          if (!is.null(block)) {
            
              ## Create contrasts between chemicals
              corfit <- duplicateCorrelation(eSub, design, 
                                             block=pD$block)
              
              if (corfit$consensus.correlation %in% c(-1, 1) |
                  is.na(corfit$consensus.correlation)) {
                  fit <- lmFit(eSub, design)
              } else {
                  fit <- lmFit(eSub, design,
                      correlation=corfit$consensus.correlation,
                      block=pD[, block])
              }
              
          } else {
              ## Fit model
              fit <- lmFit(eSub, design)
          }
          
          if(length(pDmods) == 2 & is.null(vehicle)) {
            modS <- .limmaTable(pDmods[2], fit, logCounts)
            modVec <- sub("mods", "", pDmods)
            modS$edge <- ifelse(modS$t > 0, modVec[2], modVec[1])
            modS$coef <- abs(modS$coef)
            modS$t <- abs(modS$t)
            modS <- modS[order(modS$pval),]
          }
          
          if(!is.null(vehicle)) {
            modS <- do.call(rbind, lapply(pDmods[-1], function(x) {
              modM <- .limmaTable(x, fit, logCounts, conVec, design, "mods0")
              return(modM)
            }))
            modS <- modS[order(modS$pvalV, partial = -abs(modS$coefV)),]
            modS <- modS[sign(modS$coef) == sign(modS$coefV),]
            modS <- modS[!duplicated(modS$feat),]
            modS <- modS[order(modS$pval),-seq(ncol(modS)-1, ncol(modS))]
          }
          
          if(length(pDmods) > 2 & is.null(vehicle)) {
            modS <- do.call(rbind,
                            lapply(pDmods, 
                                   .limmaTable, fit, logCounts, conVec, design))
            modS <- modS[order(modS$pval, partial = -abs(modS$coef)),]
            modS <- modS[!duplicated(modS$feat),]
            modS <- modS[order(modS$pval),]
          }
  
          ## Remove large objects
          rm(fit, design, modDF, gUnique, eSub, pD)
  
          ## Save mods as character
          modS$edge <- as.character(modS$edge)
      }
    }

    ## Return
    return(list(
        modStats = modS,
        formula = desForm
    ))
}
