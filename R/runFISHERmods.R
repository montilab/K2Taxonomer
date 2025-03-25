#' Perform Fisher test overrepresentation analysis of each subgroup
#'
#' Adds overrepresentation testing results using the output of runDGEmods().
#' @return An object of class K2.
#' @references
#'    \insertRef{reed_2020}{K2Taxonomer}
#'    \insertRef{bh}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2preproc
#' @export
#' @import Biobase

runFISHERmods <- function(K2res, genesets=NULL, qthresh=NULL,
    cthresh=NULL, ntotal=NULL) {

    ## Run checks
    .isK2(K2res)

    ## K2 algorithm
    if (length(K2results(K2res)) == 0) {
        stop("No results found. Please run K2tax() or runK2Taxonomer().\n")
    }

    ## DGE
    if (is.null(K2results(K2res)[[1]]$dge)) {
        stop("No differential analysis results found. Please run runDGEmods().\n")
    }

    ## Change meta data if new value is specific
    K2meta(K2res)$qthresh <- qthresh <- .checkK2(K2res, "qthresh",
        qthresh)
    K2meta(K2res)$cthresh <- cthresh <- .checkK2(K2res, "cthresh",
        cthresh)
    K2meta(K2res)$ntotal <- ntotal <- .checkK2(K2res, "ntotal",
        ntotal)

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
    
    ## Gene names not found in expression set
    if(nrow(K2eMatDS(K2res)) != 0) {
      gs <- rownames(K2eMatDS(K2res))
    } else {
      gs <- rownames(K2eMat(K2res))
    }
    if (length(K2genesets(K2res)) > 0 &&
        sum(unique(unlist(K2genesets(K2res))) %in%
            gs) == 0) {
      stop("No features in argument, genesets, found in data set.\n")
    }

    ## Run hyperenrichment
    K2results(K2res) <- lapply(K2results(K2res), function(x) {

        ## Get DGE results
        res <- x$dge

        ## Assign genes to each group
        one <- res[res$edge == "1", ]
        two <- res[res$edge == "2", ]

        ## For each group get a set of up and down-regulated genes
        oneUp <- one$gene[one$fdr < K2meta(K2res)$qthresh & one$coef >
            K2meta(K2res)$cthresh]
        oneDown <- one$gene[one$fdr < K2meta(K2res)$qthresh &
            one$coef < (-K2meta(K2res)$cthresh)]
        twoUp <- two$gene[two$fdr < K2meta(K2res)$qthresh & two$coef >
            K2meta(K2res)$cthresh]
        twoDown <- two$gene[two$fdr < K2meta(K2res)$qthresh &
            two$coef < (-K2meta(K2res)$cthresh)]
        sigList <- list(oneUp, oneDown, twoUp, twoDown)

        ## Run hyperenrichment
        x$gse <- lapply(sigList, function(sig) {
            enrichFram <- NULL
            if (length(sig) > 0) {
                hits <- vapply(K2genesets(K2res), function(x,
                    y) paste(intersect(x, y), collapse=","),
                    sig, FUN.VALUE=character(1))
                nhits <- vapply(K2genesets(K2res), function(x,
                    y) length(intersect(x, y)), sig, FUN.VALUE=integer(1))
                ndrawn <- length(sig)
                ncats <- vapply(K2genesets(K2res), length, FUN.VALUE=integer(1))
                nleft <- K2meta(K2res)$ntotal - ncats
                pval <- mapply(function(nh, ns, ng, nt) {
                  tr <- ng - nh
                  bl <- ns - nh
                  tl <- nt - tr - bl - nh
                  fisher.test(matrix(c(tl, bl, tr, nh), ncol = 2), 
                              alternative = "greater")$p.value
                }, nhits,  ncats, ndrawn, K2meta(K2res)$ntotal)
                enrichFram <- data.frame(category=names(K2genesets(K2res)),
                    pval=pval, fdr = NA, nhits=nhits, ndrawn=ndrawn,
                    ncats=ncats, ntot=K2meta(K2res)$ntotal,
                    hits=hits, stringsAsFactors=FALSE)
                enrichFram <- enrichFram[order(enrichFram$pval),]
            }
            return(enrichFram)
        })
        names(x$gse) <- c("g1_up", "g1_down", "g2_up", "g2_down")
        return(x)
    })
    
    ## Fix FDR values
    K2res <- .fixFDR(K2res, "gse")

    ## Set gene2pathway and make K2gSet empty
    K2gene2Pathway(K2res) <- getGenePathways(K2genesets(K2res))
    K2gMat(K2res) <- new("matrix")

    return(K2res)
}
