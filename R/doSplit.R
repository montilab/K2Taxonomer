.doSplit <- function(dataMatrix, nBoots, nFeats, featMetric,
    recalcDataMatrix, clustFunc, clustCors, clustList, linkage,
    K2res) {

    ## Order the chemicals in alphbetical order to generate
    ## predictable splits
    dataMatrix <- dataMatrix[, order(colnames(dataMatrix))]

    ## If there are only two observation set known values for
    ## split
    if (ncol(dataMatrix) == 2) {

        node_stability=0.5
        mod_stability=c(0, 0)
        Dn=0.5 * sqrt(2)
        mat <- matrix(c(1, -1, -1, 1), 2, 2)
        colnames(mat) <- rownames(mat) <- colnames(dataMatrix)
        mods <- c("1", "2")
        names(mods) <- colnames(dataMatrix)
        propBoots <- 1

    } else {

        if (recalcDataMatrix | featMetric == "F") {
            eSetSub <- K2eSet(K2res)[, pData(K2eSet(K2res))[,
                K2meta(K2res)$cohorts] %in% colnames(dataMatrix)]
            dgeSeg <- suppressWarnings(.dgeWrapper(eSetSub,
                    cohorts=K2meta(K2res)$cohorts,
                    vehicle=K2meta(K2res)$vehicle,
                    covariates=K2meta(K2res)$covariates,
                    use=K2meta(K2res)$use,
                    Fstat=TRUE))
            dataMatrix <- dgeSeg$dataMatrix
        }

        ## Get consensus module
        if (featMetric == "mad") {
            SCORE <- apply(dataMatrix, 1, mad)
        }
        if (featMetric == "sd") {
            SCORE <- apply(dataMatrix, 1, sd)
        }
        if (featMetric == "Qn") {
            SCORE <- apply(dataMatrix, 1, Qn)
        }
        if (featMetric == "Sn") {
            SCORE <- apply(dataMatrix, 1, Sn)
        }
        if (featMetric == "F") {
            SCORE <- dgeSeg$Fvec
        }
        if (featMetric == "square") {
            SCORE <- apply(dataMatrix, 1, function(x) sum(x^2))
        }

        ## Get number of features
        modVec <- do.call(c, mclapply(seq(nBoots), function(x,
            clustFunc, SCORE, nFeats, dataMatrix, clustList) {
            SCOREboot <- sort(sample(SCORE, length(SCORE), replace=TRUE),
                decreasing=TRUE)[seq_len(nFeats)]
            Dboot <- dataMatrix[names(SCOREboot), ]
            mods <- clustFunc(Dboot, clustList)
            return(mods)
        }, clustFunc, SCORE, nFeats, dataMatrix, clustList,
            mc.cores=clustCors))

        ## Create table of results
        modTab <- sort(table(modVec), decreasing=TRUE)

        ## Get all clustering
        modList <- strsplit(modVec, "")

        ## Get stability statistics for each sample
        modataMatrix <- do.call(rbind, modList)
        colnames(modataMatrix) <- colnames(dataMatrix)

        ## For each cluster get mean cosine distance in each cluster
        mat <- matrix(NA, ncol(modataMatrix), ncol(modataMatrix))
        for (i in seq_len(ncol(modataMatrix))) {
            for (j in seq_len(ncol(modataMatrix))) {
                ## Convert to -1 and 1
                vec1 <- c(-1, 1)[as.numeric(modataMatrix[, i] ==
                    "2") + 1]
                vec2 <- c(-1, 1)[as.numeric(modataMatrix[, j] ==
                    "2") + 1]
                mat[i, j] <- t(vec1) %*% vec2/length(vec1)
            }
        }
        colnames(mat) <- rownames(mat) <- colnames(dataMatrix)

        ## Get mods
        mods <- cutree(hclust(as.dist(1 - mat), method=linkage),
            k=2)

        ## Get proportion of boots that had best split
        modString <- paste(mods, collapse="")
        if (modString %in% names(modTab)) {
            propBoots <- modTab[modString]/nBoots
            names(propBoots) <- NULL
        } else {
            propBoots <- 0
        }

        ## Get cluster stability (first eigenvalue in each)
        mod_stability <- vapply(seq_len(2), function(x, mods,
            mat) {

            ## Get matrix
            matMod <- as.matrix(mat[mods == x, mods == x, drop=FALSE])

            ## Calculate number of eigenvalues to 90% variance
            stabVals <- svd(matMod)$d
            CDstab <- cumsum(stabVals)/ncol(matMod)
            CDnull <- cumsum(rep(1, ncol(matMod))/ncol(matMod))
            max(CDstab - CDnull)

        }, mods, mat, FUN.VALUE=double(1))

        ## Get stability of the entire node
        stabVals <- svd(mat)$d
        CDstab <- cumsum(stabVals)/ncol(mat)
        CDnull <- cumsum(rep(1, ncol(mat))/ncol(mat))
        node_stability <- max(CDstab - CDnull)
        Dn <- node_stability * sqrt(ncol(mat))

    }

    ## Get mod and sample stability
    stabList <- list(node=node_stability, Dn=Dn, clusters=mod_stability,
        samples=as.dist(mat))

    return(list(mods=mods, propBoots=propBoots, stability=stabList))
}
