.doSplit <- function(dataM, K2r) {

    ## Order the chemicals in alphbetical order to generate
    ## predictable splits
    clabs <- sort(colnames(dataM))
    K2data(K2r) <- dataM[, clabs]

    ## If there are only two observation set known values for
    ## split
    if (ncol(dataM) == 2) {

        node_stability=0.5
        mod_stability=c(0, 0)
        Dn=0.5 * sqrt(2)
        mat <- matrix(c(1, -1, -1, 1), 2, 2)
        colnames(mat) <- rownames(mat) <- colnames(dataM)
        mods <- c("1", "2")
        names(mods) <- colnames(dataM)
        propBoots <- 1

    } else {

        if (K2meta(K2r)$recalcDataMatrix | K2meta(K2r)$featMetric == "F") {
            eKeep <- rownames(K2colData(K2r))[K2colData(K2r)[,
                K2meta(K2r)$cohorts] %in% colnames(dataM)]
            dgeSeg <- suppressWarnings(.dgeWrapper(K2r, eKeep, Fstat=TRUE))
            K2data(K2r) <- dgeSeg$dataM
        }

        ## Get consensus module
        if (K2meta(K2r)$featMetric == "mad") {
            SCORE <- apply(K2data(K2r), 1, mad)
        }
        if (K2meta(K2r)$featMetric == "sd") {
            SCORE <- apply(K2data(K2r), 1, sd)
        }
        if (K2meta(K2r)$featMetric == "Qn") {
            SCORE <- apply(K2data(K2r), 1, Qn)
        }
        if (K2meta(K2r)$featMetric == "Sn") {
            SCORE <- apply(K2data(K2r), 1, Sn)
        }
        if (K2meta(K2r)$featMetric == "F") {
            SCORE <- dgeSeg$Fvec
        }
        if (K2meta(K2r)$featMetric == "square") {
            SCORE <- apply(K2data(K2r), 1, function(x) sum(x^2))
        }
      
        ## Cluster data
        modVec <- do.call(c, bplapply(seq(K2meta(K2r)$nBoots), function(x, K2r) {
            SCOREboot <- sort(sample(SCORE, length(SCORE), replace=TRUE),
                decreasing=TRUE)[seq_len(K2meta(K2r)$nFeats)]
            mods <- K2meta(K2r)$clustFunc(clabs, names(SCOREboot), K2r)
            return(mods)
        }, K2r, BPPARAM = K2meta(K2r)$BPPARAM_useCors))

        ## Create table of results
        modTab <- sort(table(modVec), decreasing=TRUE)

        ## Get all clustering
        modList <- strsplit(modVec, "")

        ## Get stability statistics for each sample
        modataM <- do.call(rbind, modList)
        colnames(modataM) <- colnames(K2data(K2r))

        ## For each cluster get mean cosine distance in each cluster
        mat <- matrix(NA, ncol(modataM), ncol(modataM))
        for (i in seq_len(ncol(modataM))) {
            for (j in seq_len(ncol(modataM))) {
                ## Convert to -1 and 1
                vec1 <- c(-1, 1)[as.numeric(modataM[, i] ==
                    "2") + 1]
                vec2 <- c(-1, 1)[as.numeric(modataM[, j] ==
                    "2") + 1]
                mat[i, j] <- t(vec1) %*% vec2/length(vec1)
            }
        }
        colnames(mat) <- rownames(mat) <- colnames(K2data(K2r))

        ## Get mods
        mods <- cutree(hclust(as.dist(1 - mat), method=K2meta(K2r)$linkage),
            k=2)

        ## Get proportion of boots that had best split
        modString <- paste(mods, collapse="")
        if (modString %in% names(modTab)) {
            propBoots <- modTab[modString]/K2meta(K2r)$nBoots
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
