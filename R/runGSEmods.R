#' Perform hyperenrichment of split specific gene signatures
#'
#' Adds hyperenrichment analysis results to the output of runDGEmods().
#' @param K2res An object of class K2. The output of runDGEmods().
#' @param genesets A named list of feature IDs
#' @param qthresh A numeric value between 0 and 1 of the FDR cuttoff to define
#' feature sets.
#' @param cthresh A positive value for the coefficient cuttoff to define
#' feature sets.
#' @return An object of class K2.
#' @keywords clustering
#' @export
#' @import Biobase
#' @examples
#' runGSEmods(K2res)

runGSEmods <- function(K2res, genesets = NULL, qthresh = NULL, cthresh = NULL, ntotal = NULL) {


    ## Change meta data if new value is specific
    K2meta(K2res)$qthresh <- qthresh <- .checkK2(K2res, "qthresh", qthresh)
    K2meta(K2res)$cthresh <- cthresh <- .checkK2(K2res, "cthresh", cthresh)
    K2meta(K2res)$ntotal <- ntotal <- .checkK2(K2res, "ntotal", ntotal)
    
    ## Check K2 object
    k2Check <- .checkK2(K2res)

    ## Set genesets if not null
    if (!is.null(genesets)) {
        K2genesets(K2res) <- genesets
    }

    ## If genesets is empty then stop
    if (length(K2genesets(K2res)) == 0) {
        stop("No named list of genesets provided")
    }

    ## Run hyperenrichment
    K2results(K2res) <- lapply(K2results(K2res), function(x) {

        ## Get DGE results
        res <- x$dge

        ## Assign genes to each group
        one <- res[res$mod == "1", ]
        two <- res[res$mod == "2", ]

        ## For each group get a set of up and down-regulated genes
        oneUp <- one$gene[one$fdr < K2meta(K2res)$qthresh & one$coef > K2meta(K2res)$cthresh]
        oneDown <- one$gene[one$fdr < K2meta(K2res)$qthresh & one$coef < (-K2meta(K2res)$cthresh)]
        twoUp <- two$gene[two$fdr < K2meta(K2res)$qthresh & two$coef > K2meta(K2res)$cthresh]
        twoDown <- two$gene[two$fdr < K2meta(K2res)$qthresh & two$coef < (-K2meta(K2res)$cthresh)]
        sigList <- list(oneUp, oneDown, twoUp, twoDown)

        ## Run hyperenrichment
        x$gse <- lapply(sigList, function(sig) {
            enrichFram <- NULL
            if (length(sig) > 0) {
                hits <- vapply(K2genesets(K2res), function(x, y) paste(intersect(x,
                  y), collapse = ","), sig, FUN.VALUE = character(1))
                nhits <- vapply(K2genesets(K2res), function(x, y) length(intersect(x,
                  y)), sig, FUN.VALUE = integer(1))
                ndrawn <- length(sig)
                ncats <- vapply(K2genesets(K2res), length, FUN.VALUE = integer(1))
                nleft <- K2meta(K2res)$ntotal - ncats
                pval <- phyper(q = nhits - 1, m = ncats, n = nleft, k = ndrawn, lower.tail = F)
                enrichFram <- data.frame(category = names(K2genesets(K2res)), pval = pval,
                  nhits = nhits, ndrawn = ndrawn, ncats = ncats, ntot = K2meta(K2res)$ntotal,
                  hits = hits, stringsAsFactors = F)
            }
            return(enrichFram)
        })
        names(x$gse) <- c("g1_up", "g1_down", "g2_up", "g2_down")
        return(x)
    })

    ## Calculate and merge FDR values
    pValueDF <- data.frame(pval = unlist(lapply(K2results(K2res), function(x) lapply(x$gse,
        function(y) y$pval))))
    pValueDF$fdr <- p.adjust(pValueDF$pval, method = "BH")
    pValueDF <- unique(pValueDF)
    K2results(K2res) <- lapply(K2results(K2res), function(x, pValueDF) {
        x$gse <- lapply(x$gse, function(y, pValueDF) {
            y <- merge(y, pValueDF)
            if (nrow(y) > 0) {
                y <- y[, c("category", "pval", "fdr", "nhits", "ndrawn", "ncats",
                  "ntot", "hits")]
            }
            return(y)
        }, pValueDF)
        return(x)
    }, pValueDF)

    ## Set gene2pathway and make K2gSet empty
    K2gene2Pathway(K2res) <- getGenePathways(K2genesets(K2res))
    K2gSet(K2res) <- ExpressionSet()

    return(K2res)
}
