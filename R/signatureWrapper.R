.signatureWrapper <- function(eSet, cohorts, mods, vehicle=NULL,
    covariates=NULL, block=NULL, logCounts=FALSE) {

    ## Remove vehicle from mods and make a data frame
    if (!is.null(vehicle)) {
        mods <- mods[names(mods) != vehicle]
    }
    mods <- data.frame(mods=as.character(mods), GROUP=names(mods),
        stringsAsFactors=FALSE)

    modStats <- NULL

    ## Run if there are subgroups to compare and at least three
    ## observations
    if (length(unique(mods$mods)) > 1 && (nrow(mods) > 2 | !is.null(cohorts))) {

        ## If replicates in data get unique cohorts
        if (is.null(cohorts)) {
            cohorts <- "GROUP"
            pData(eSet)[, cohorts] <- colnames(eSet)
        }

        ## Get unique groups
        gUnique <- unique(pData(eSet)[, cohorts])

        if (!is.null(vehicle)) {
            gUnique <- gUnique[gUnique != vehicle]
        }

        ## Subset data for mods
        eSub <- eSet[, pData(eSet)[, cohorts] %in% c(vehicle,
            mods$GROUP)]

        ## Drop levels
        pData(eSub) <- droplevels(pData(eSub))

        ## Create new variable in pData by merging with mods
        ## data.frame
        pData(eSub)$GROUP <- factor(pData(eSub)[, cohorts], levels=c(vehicle,
            mods$GROUP))
        pData(eSub)$Rownames <- rownames(pData(eSub))
        pD <- merge(pData(eSub), mods, all.x=TRUE)
        rownames(pD) <- pD$Rownames
        pD <- pD[colnames(eSub), , drop=FALSE]

        ## Add vehicle mod
        pD$mods[is.na(pD$mods)] <- "0"
        pD$mods <- as.factor(pD$mods)
        pDmods <- paste0("X", as.character(unique(pD$mods)))

        ## Add back to eSet
        pData(eSub) <- pD

        ## Need to run different analysis if there are cohorts or not
        if (!is.null(cohorts)) {

            ## Create design matrix
            design <- model.matrix(as.formula(paste0("~ 0 +",
                                                     "mods", .formatCov(covariates))), data=pData(eSub))
            
            ## Check that model is full rank with covariates, if not model w/o covariates
            if(!is.fullrank(design)) {
                design <- model.matrix(as.formula(paste0("~ 0 +",
                                                         "mods", .formatCov(NULL))), data=pData(eSub))
            }
            colnames(design) <- sub("mods", "X", colnames(design))

            ## Run duplicateCorrelation
            if (!is.null(block)) {
                ## Create contrasts between chemicals
                corfit <- duplicateCorrelation(eSub, design,
                    block=pData(eSub)[, block])

                if (corfit$consensus.correlation %in% c(-1, 1) |
                    is.na(corfit$consensus.correlation)) {
                    fit <- lmFit(eSub, design)
                } else {
                    fit <- lmFit(eSub, design,
                        correlation=corfit$consensus.correlation,
                        block=pData(eSub)[, block])
                }
            } else {
                ## Fit model
                fit <- lmFit(eSub, design)
            }

            ## Fit contrasts
            modFit <- lapply(paste0("X", unique(mods$mods)),
                function(x, design, fit) {
                    conString <- paste0(x, " - (", paste(pDmods[pDmods !=
                    x], collapse="+"), ")/", sum(pDmods !=
                    x))
                    contrasts <- makeContrasts(contrasts=conString,
                    levels=design)
                    contFit <- suppressWarnings(topTable(eBayes(
                        contrasts.fit(fit, contrasts),
                        trend=logCounts, robust=logCounts),
                    number=Inf, sort.by="none"))
                    contFit <- contFit[, c("logFC", "AveExpr",
                    "t", "P.Value", "adj.P.Val", "B")]
                    contFit$edge <- sub("X", "", x)
                    return(contFit)
                }, design, fit)

        } else {

            design <- model.matrix(as.formula(paste0("~ 0 + ",
                                                     "GROUP", .formatCov(covariates))), data=pData(eSub))
            
            ## Check that model is full rank with covariates, if not model w/o covariates
            if(!is.fullrank(design)) {
                design <- model.matrix(as.formula(paste0("~ 0 +",
                                                         "GROUP", .formatCov(NULL))), data=pData(eSub))
            }
            colnames(design) <- sub("GROUP", "X", colnames(design))

            ## Run duplicateCorrelation
            if (!is.null(block)) {
                ## Create contrasts between chemicals
                corfit <- duplicateCorrelation(eSub, design,
                    block=pData(eSub)[, block])

                if (corfit$consensus.correlation %in% c(-1, 1) |
                    is.na(corfit$consensus.correlation)) {
                    fit <- lmFit(eSub, design)
                } else {
                    fit <- lmFit(eSub, design,
                        correlation=corfit$consensus.correlation,
                        block=pData(eSub)[, block])
                }
            } else {
                ## Fit model
                fit <- lmFit(eSub, design)
            }

            ## Create contrasts strings
            modsFull <- unique(pD[, c("GROUP", "mods")])
            modsTable <- table(modsFull$mods)
            cVec <- vapply(names(modsTable), function(x) {
                paste0("(", paste(paste0("X", modsFull$GROUP[modsFull$mods ==
                    x], collapse="+")), ")/", sum(modsFull$mods ==
                    x))
            }, FUN.VALUE=integer(1))

            ## Run each contrast
            modFit <- lapply(as.character(unique(mods$mods)),
                function(x, cVec, design, fit) {
                    conString <- paste0(cVec[x], " - (",
                    paste(cVec[names(cVec) !=x], collapse="+"), ")/",
                    sum(names(cVec) != x))
                    contrasts <- makeContrasts(contrasts=conString,
                    levels=design)
                    contFit <- topTable(eBayes(contrasts.fit(fit,
                    contrasts), trend=logCounts, robust=logCounts),
                    number=Inf, sort.by="none")
                    contFit <- contFit[, c("logFC", "AveExpr",
                    "t", "P.Value", "adj.P.Val", "B")]
                    contFit$edge <- x
                    return(contFit)
                }, cVec, design, fit)
        }

        ## Create vector of where to assign result
        if (is.null(vehicle)) {
            if (length(modFit) == 2) {
                one2 <- c(1, 2)[as.numeric(modFit[[1]]$t < 0) +
                    1]
            } else {
                one2 <- vapply(seq(nrow(modFit[[1]])), function(row,
                    modFit) {
                    which.max(vapply(seq(length(modFit)), function(g,
                    modFit, row) {
                    modFit[[g]][row, "t"]
                    }, modFit, row, FUN.VALUE=numeric(1)))
                }, modFit, FUN.VALUE=integer(1))
            }
        } else {
            one2 <- c(1, 2)[as.numeric(vapply(seq(nrow(modFit[[1]])),
                function(row, modFit) {
                    modFit[[1]]$P.Value[row] > modFit[[2]]$P.Value[row]
                }, modFit, FUN.VALUE=logical(1))) + 1]
        }

        ## Create just one data.frame
        modStats <- as.data.frame(t(vapply(seq_len(nrow(modFit[[1]])),
            function(row, one2, modFit) {
                unlist(as.numeric(modFit[[one2[row]]][row, ]))
            }, one2, modFit, FUN.VALUE=double(7))))
        colnames(modStats) <- colnames(modFit[[1]])
        rownames(modStats) <- rownames(modFit[[1]])

        ## Order by p-value
        modStats <- modStats[order(modStats$P.Value), ]
        modStats$adj.P.Val <- p.adjust(modStats$P.Value, method="BH")

        ## Remove large objects
        rm(fit, modFit, design, mods, gUnique, one2, eSet, eSub,
            pD)

        ## Save mods as character
        modStats$edge <- as.character(modStats$edge)

        ## Change column names
        colnames(modStats) <- c("coef", "mean", "t", "pval",
            "fdr", "B", "edge")
    }

    ## Return
    return(modStats)
}
