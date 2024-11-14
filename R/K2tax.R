#' Function to K2 Taxonomer partitioning
#'
#' This function performs K2 Taxonomer procedure only. Arguments values are
#' extracted from K2meta(K2res) unless specified otherwise.
#' @param K2res A K2 class object.
#' @return An object of class K2.
#' @references
#'    \insertRef{reed_2020}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2preproc
#' @export
#' @import parallel
#' @import robustbase

K2tax <- function(K2res, nFeats=NULL, featMetric=NULL, recalcDataMatrix=NULL,
    nBoots=NULL, clustFunc=NULL, useCors=NULL, clustList=NULL,
    linkage=NULL, oneoff=NULL, stabThresh=NULL) {

    ## Run checks
    .isK2(K2res)
  
    ## Change meta data if new value is specific
    K2meta(K2res)$nFeats <- .checkK2(K2res, "nFeats", nFeats)
    K2meta(K2res)$featMetric <- .checkK2(K2res, "featMetric", featMetric)
    K2meta(K2res)$recalcDataMatrix <- .checkK2(K2res, "recalcDataMatrix", 
                                               recalcDataMatrix)
    K2meta(K2res)$nBoots <- .checkK2(K2res, "nBoots", nBoots)
    K2meta(K2res)$clustFunc <- .checkK2(K2res, "clustFunc", clustFunc)
    K2meta(K2res)$useCors <- .checkK2(K2res, "useCors", useCors)
    K2meta(K2res)$clustList <- .checkK2(K2res, "clustList", clustList)
    K2meta(K2res)$linkage <- .checkK2(K2res, "linkage", linkage)
    K2meta(K2res)$oneoff <- .checkK2(K2res, "oneoff", oneoff)
    K2meta(K2res)$stabThresh <- .checkK2(K2res, "stabThresh", stabThresh)
    
    ## Set cluster function
    if(class(K2meta(K2res)$clustFunc) == "character") {
      clustFunc <- get(K2meta(K2res)$clustFunc)
    }
    
    ## Get parallel parameters
    if(K2meta(K2res)$useCors > 1) {
      bF <- get(class(bpparam())[[1]])
      K2meta(K2res)$BPPARAM_useCors <- bF(workers = K2meta(K2res)$useCors)
    } else {
      K2meta(K2res)$BPPARAM_useCors <- SerialParam()
    }

    ## Create Splits of the data
    taxList <- list(list(colnames(K2data(K2res))))
    resList <- stabList <- list(list(NULL))
    iter <- 1
    while (max(unlist(lapply(taxList[[iter]], function(x) {

        length(x[x != "Vehicle"])

    }))) > (2 - K2meta(K2res)$oneoff)) {
        outList <- lapply(taxList[[iter]], function(samps, K2res) {
            if (length(samps) > 1) {
                Dsub <- K2data(K2res)[, samps]
                outList <- .doSplit(Dsub, K2res)

                ## Get minimum size
                minSize <- min(table(outList$mods))

                ## Get determinant
                nodeStab <- outList$stability$node

                ## Check stopping criteria
                if ((!K2meta(K2res)$oneoff & minSize == 1 & iter != 1) | 
                    (nodeStab <= K2meta(K2res)$stabThresh & iter != 1)) {
                    outList <- list(mods=samps, propBoots=NULL,
                    clustStab=NULL)
                }

            } else {
                outList <- list(mods=samps, propBoots=NULL,
                    clustStab=NULL)
            }
            return(outList)
        }, K2res)
        iter <- iter + 1
        taxList[[iter]] <- list()
        outMods <- lapply(outList, function(x) x[[1]])
        slot <- 1
        for (i in seq_len(length(outMods))) {
            taxList[[iter]][[slot]] <- names(outMods[[i]])[outMods[[i]] ==
                1]
            slot <- slot + 1
            taxList[[iter]][[slot]] <- names(outMods[[i]])[outMods[[i]] ==
                2]
            slot <- slot + 1
        }

        if (length(taxList[[iter]]) == 0) {
            taxList[[iter]] <- list(NULL)
        }

        resList[[iter]] <- list()
        stabList[[iter]] <- list()
        outRes <- lapply(outList, function(x) x[[2]])
        outStab <- lapply(outList, function(x) x[[3]])
        slot <- 1
        for (i in seq_len(length(outRes))) {
            resList[[iter]][[slot]] <- outRes[[i]]
            stabList[[iter]][[slot]] <- outStab[[i]]
            slot <- slot + 1
        }
    }
    resList <- resList[-1]
    stabList <- stabList[-1]

    if (is.null(unlist(taxList[[iter]])))
        taxList <- taxList[-iter]

    ## Get instances where the split had > 1 samples in cluster
    modList <- lapply(taxList[-1], function(x) {
        combs <- seq_len(length(x)/2)
        combList <- list()
        for (i in combs) {
            combList[[i]] <- list()
            combList[[i]][[1]] <- x[[i * 2 - 1]]
            combList[[i]][[2]] <- x[[i * 2]]
            if (length(combList[[i]]) > 0) {
                if (length(combList[[i]][[1]]) < (2 - K2meta(K2res)$oneoff) |
                    length(combList[[i]][[2]]) < (2 - K2meta(K2res)$oneoff))
                    combList[[i]] <- list()
            }
        }
        combList
    })

    ## Create list of modules and bootstrap stats
    K2c <- lapply(seq_len(length(modList)), function(x) {
        lapply(seq_len(length(modList[[x]])), function(y) {
            modSub <- modList[[x]][[y]]
            if (length(modSub) > 0) {
                resSub <- resList[[x]][[y]]
                stabSub <- stabList[[x]][[y]]
                return(list(obs=modSub, bootP=resSub, stability=stabSub))
            } else {
                return(NULL)
            }
        })
    })
    K2c <- unlist(K2c, recursive=FALSE)
    K2c <- K2c[!unlist(lapply(K2c, is.null))]

    ## Add names
    len <- length(K2c)
    nAlphabets <- ceiling(len/length(LETTERS))
    ALPHABETS <- unlist(lapply(seq_len(nAlphabets), function(x) vapply(LETTERS,
        function(y) paste(rep(y, x), collapse=""), FUN.VALUE=character(1))))
    names(K2c) <- ALPHABETS[seq_len(length(K2c))]

    ## Add to K2res and return
    K2results(K2res) <- K2c
    return(K2res)

}
