#' Extract table of differential results from 'K2' object
#'
#' Create table differential analysis results from 'K2' object.
#' @return A data.frame object with the following columns:
#' \itemize{
#'  \item{gene: }{The identifier of the gene}
#'  \item{coef: }{The difference between the means of each subgroup at
#'  a given partition}
#'  \item{mean: }{The mean across all observations at the given
#'  partition}
#'  \item{t: }{The test statistic estimated by differential analysis}
#'  \item{pval: }{The p-value estimated by differential analysis}
#'  \item{fdr: }{The multiple hypothesis corrected fdr p-value, adjusted across
#'  all partitions}
#'  \item{edge: }{Indication of which subgroup the gene was assigned at a given
#'  partition}
#'  \item{node: }{The identifier of the partition}
#'  \item{direction: }{The direction of coefficient for the assigned gene}
#' }
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{limma}{K2Taxonomer}
#'  \insertRef{bh}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2tax
#' @export

getDGETable <- function(K2res) {

    ## Run checks
    .isK2(K2res)

    K2resList <- K2results(K2res)

    ## Format single-sample enrichment results
    dgeTable <- do.call(rbind, lapply(names(K2resList), function(x,
        K2resList) {

        ## Get tables
        CompTab <- K2resList[[x]]$dge

        if (!is.null(CompTab)) {
            ## Add group information and extract tables
            CompTab$node <- x

            ## Add direction information
            CompTab$direction <- c("down", "up")[as.numeric(CompTab$coef >
                0) + 1]
        }

        return(CompTab)
    }, K2resList))

    ## Sort by p-value
    dgeTable <- dgeTable[order(dgeTable$pval, partial = -abs(dgeTable$coef)), ]

    rownames(dgeTable) <- NULL

    return(dgeTable)
}
