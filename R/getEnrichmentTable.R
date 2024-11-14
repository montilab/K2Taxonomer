#' Extract table of enrichment results from 'K2' object
#'
#' Create table fisher- and gene set scoring enrichment results from 'K2' object.
#' @return A data.frame object with the following columns:
#' \itemize{
#'    \item{category: }{The user-specified gene set names}
#'    \item{node: }{The identifier of the partition}
#'    \item{edge: }{Indication of which subgroup the gene was assigned at a
#'    given partition}
#'    \item{direction: }{The direction of coefficient for the assigned gene set}
#'    \item{pval_fisher: }{The p-value of the enrichment test comparing the
#'    subgroup-assigned gene set to the user-specified gene set}
#'    \item{fdr_fisher: }{The multiple hypothesis corrected FDR
#'    (Benjamini-Hochberg) p-value of enrichment, adjusted across all
#'    partitions}
#'    \item{nhits: }{The intersection of subgroup-assigned genes and the
#'    user-specified gene set}
#'    \item{ndrawn: }{The number of subgroup-assigned genes}
#'    \item{ncats: }{The number of genes in the user-specified gene set}
#'    \item{ntot: }{The background population of possible genes}
#'    \item{pval_limma: }{The p-value estimated by the `limma` R package}
#'    \item{fdr_limma: }{The multiple hypothesis corrected FDR
#'    (Benjamini-Hochberg) p-value of differential analysis, adjusted across all
#'    partitions}
#'    \item{coef: }{The difference between the means of each subgroup at a given
#'    partition}
#'    \item{mean: }{The mean across all observations at the given partition}
#'    \item{t: }{The test statistic estimated by the `limma` R package}
#'    \item{hits: }{A comma deliminated string of intersection gene identifers}
#' }
#' @references
#'    \insertRef{reed_2020}{K2Taxonomer}
#'    \insertRef{limma}{K2Taxonomer}
#'    \insertRef{bh}{K2Taxonomer}
#'    \insertRef{gsva}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2tax
#' @export

getEnrichmentTable <- function(K2res) {

    # Run checks
    .isK2(K2res)

    K2resList <- K2results(K2res)

    ## Format enrichment table
    EnrTable <- do.call(rbind, lapply(names(K2resList), function(x,
        K2resList) {

        ## Get GSE tables
        GSEtab <- NULL
        if (!is.null(K2resList[[x]]$dge)) {
            GSEtabList <- K2resList[[x]]$gse

            ## Add group informaiton and extract tables
            GSEtab <- do.call(rbind, lapply(names(GSEtabList),
                function(y) {
                    GSEsub <- GSEtabList[[y]]
                    if (!is.null(GSEsub) && nrow(GSEsub) > 0) {
                      GSEsub$node <- x    ## Add node ID
                      GSEsub$edge <- gsub("g|_up|_down", "", y)    ## Add group ID
                      GSEsub$direction <- gsub("g1_|g2_", "", y)    ## Add    dir
                    } else {
                      GSEsub <- NULL
                    }
                    return(GSEsub)
              }))

            ## Make pval and fdr column names unique
            colnames(GSEtab)[colnames(GSEtab) %in% c("pval",
                "fdr")] <- paste(colnames(GSEtab)[colnames(GSEtab) %in%
                c("pval", "fdr")], "fisher", sep="_")
        }
        return(GSEtab)

    }, K2resList))

    if (nrow(EnrTable) == 0) {
        EnrTable <- data.frame(category=NA, pval_fisher=NA,
            fdr_fisher=NA, nhits=NA, ndrawn=NA, ncats=NA,
            ntot=NA, hits=NA, node=NA, edge=NA, direction=NA)
    }

    ## Format single-sample enrichment results
    ssEnrTable <- do.call(rbind, lapply(names(K2resList), function(x,
        K2resList) {

        SSGSEAtab <- NULL
        if (!is.null(K2resList[[x]]$dge)) {
            ## Get SSGSEA tables
            SSGSEAtab <- K2resList[[x]]$dsse

            ## Add group information and extract tables
            SSGSEAtab$node <- x

            ## Add direction information
            SSGSEAtab$direction <- c("down", "up")[as.numeric(SSGSEAtab$t >
                0) + 1]

            ## Make pval and fdr column names unique
            colnames(SSGSEAtab)[colnames(SSGSEAtab) %in% c("pval",
                "fdr")] <- paste(colnames(SSGSEAtab)[colnames(SSGSEAtab) %in%
                c("pval", "fdr")], "limma", sep="_")
        }

        return(SSGSEAtab)
    }, K2resList))

    ## Merge the two and sort by fisher p-value
    EnrTable <- merge(EnrTable, ssEnrTable, all=TRUE)
    EnrTable <- EnrTable[!is.na(EnrTable$category), ]

    ## Sort columns
    EnrTable <- EnrTable[, c("category", "node", "edge", "direction",
        "pval_fisher", "fdr_fisher", "nhits", "ndrawn", "ncats",
        "ntot", "pval_limma", "fdr_limma", "coef", "mean", "t",
        "hits")]

    ## Sort by p-value
    EnrTable <- EnrTable[order(EnrTable$pval_limma, partial = EnrTable$pval_fisher),]

    rownames(EnrTable) <- NULL

    return(EnrTable)
}
