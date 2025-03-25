.signatureWrapper_MAST <- function(K2r, mods, GENE = TRUE) {
  
    # Get objects from K2r object
    pD <- K2colData(K2r)
    cohorts <- K2meta(K2r)$cohorts
    vehicle <- K2meta(K2r)$vehicle
    variables <- K2meta(K2r)$variables
    block <- K2meta(K2r)$block
    useCors <- K2meta(K2r)$useCors
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
    
    ## Remove vehicle from mods and make a data frame
    if (!is.null(vehicle)) {
        mods <- mods[names(mods) != vehicle]
    }
    modDF <- data.frame(mods=as.character(mods), GROUP=names(mods),
        stringsAsFactors=FALSE)

    modS <- NULL
    desForm <- NULL

    ## Run if there are subgroups to compare and at least three
    ## observations
    if (length(unique(modDF$mods)) > 1 && (nrow(modDF) > 2 | !is.null(cohorts))) {

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

        ## Add vehicle mod
        pD$mods[is.na(pD$mods)] <- "0"
        pD$mods <- as.factor(pD$mods)
        pDmods <- paste0("mods", as.character(sort(unique(pD$mods))))
        
        ## If block add to pD
        if (!is.null(block)) {
          pD$block <- as.character(pD[, block])
        }
        
        ## Remove lowly expressed genes
        if(DGEexpThreshold > 0) {
          eSub <- .K2filterGenes(eSub, pD$mods, DGEexpThreshold)
        }
        
        ## Get mean expression values
        eMeans <- rowMeans(eSub)
        
        ## Create single cell object
        SCA <- suppressMessages(FromMatrix(as.matrix(eSub), pD))

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
        
        options(mc.cores = useCors)
        
        ## If block run mixed model
        if (!is.null(block)) {
          
          # Add random effect
          desFormMM <- paste0(desForm, "+(1|block)")
          fit <- suppressMessages(zlm(as.formula(desFormMM), SCA, 
                                      parallel = TRUE,
                                      method='glmer', ebayes=FALSE))
        
        } else {
        ## Fit model
          fit <- suppressMessages(zlm(as.formula(desForm), SCA, parallel = TRUE))
        }
        
        if(length(pDmods) == 2 & is.null(vehicle)) {
          modS <- .mastTable(pDmods[2], fit, eMeans)
          modVec <- sub("mods", "", pDmods)
          modS$edge <- ifelse(modS$coef > 0, modVec[2], modVec[1])
          modS$coef <- abs(modS$coef)
          modS <- modS[order(modS$pval),]
        }
        
        if(!is.null(vehicle)) {
          modS <- do.call(rbind, lapply(pDmods[-1], function(x) {
            modM <- .mastTable(x, fit, eMeans, conVec, pDmods, design, "mods0")
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
                                 .mastTable, 
                                 fit, eMeans, conVec, pDmods, design))
          modS <- modS[order(modS$pval),]
          modS <- modS[!duplicated(modS$feat),]
        }
        rownames(modS) <- NULL

        ## Remove large objects
        rm(fit, design, modDF, gUnique, eSub, pD)
        
        options(mc.cores = 1)
    }

    ## Return
    return(list(
        modStats = modS,
        formula = desForm
    ))
}
