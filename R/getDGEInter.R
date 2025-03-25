#' Create interactive table of differential gene expression analysis results from 'K2' object
#'
#' Create interactive table differential gene expression analysis results from 'K2' object.
#' @param maxFDR Numeric. A value between 0 and 1 indicating the FDR cutoff for differential gene expressio analysis.
#' @param minDiff Numeric. A value between 0 and 1 indicating the mean difference cutoff for differential gene expressio analysis.
#' @param maxPval Numeric. A value between 0 and 1 indicating the p-value cutoff for differential gene expression analysis.
#' @param genes Character. A vector of gene identifiers to display.
#' @param nodes Character. A vector of node identifiers to display.
#' @param pagelength Numeric. Number of rows to display in each page of output.
#' @return An interactive data frame with the following columns:
#' \itemize{
#'  \item{Gene: }{The identifier of the gene}
#'  \item{Node: }{The identifier of the partition}
#'  \item{Edge: }{The identifier of partition subgroup}
#'  \item{Direction: }{The direction of coefficient for the assigned gene}
#'  \item{P Value: }{The p-value estimated by differential analysis}
#'  \item{FDR: }{The multiple hypothesis corrected FDR p-value, adjusted across
#'  all partitions}
#'  \item{Diff: }{The difference between the means of each subgroup at
#'  a given partition}
#'  \item{Mean: }{The mean across all observations at the given
#'  partition}
#' }
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{limma}{K2Taxonomer}
#'  \insertRef{bh}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2tax
#' @export

## Function to format differential results
getDGEInter <- function(K2res,
                        maxFDR = 0.05,
                        minDiff = NULL,
                        maxPval = NULL,
                        genes = NULL,
                        nodes = NULL,
                        pagelength = 50) {
  
  # Get results tables
  dgetab <- getDGETable(K2res)
  dgetab <- dgetab[,c("gene", "node", "edge", "direction", "pval", "fdr",
                          "coef", "mean")]
  colnames(dgetab) <- c("Gene", "Node", "Edge", "Direction", "P Value", "FDR",
                          "Diff", "Mean")
  
  # Handle inputs
  
  ## Genes
  genesOut <- !genes %in% dgetab$Gene
  if(sum(genesOut) > 0) {
    if(mean(genesOut) == 1) {
      stop("0 input genes found in results table")
    } else {
      warning(paste0(
        sum(genesOut), " input gene(s) not found in results table (", signif(mean(genesOut)*100, 2), "% missing)\n",
        paste(genes[genesOut], collapse = ", ")
      ))
      genes <- genes[!genesOut]
    }
  }
  
  # Filter by gene input
  if(!is.null(genes)) {
    cat("Filtering for", length(genes), "genes\n")
    dgetab <- dgetab[dgetab$Gene %in% genes,]
  }
  
  ## Nodes
  nodesOut <- !nodes %in% dgetab$Node
  if(sum(nodesOut) > 0) {
    if(mean(nodesOut) == 1) {
      stop("0 input nodes found in results table")
    } else {
      warning(paste0(
        sum(nodesOut), " input nodes(s) not found in results table (", signif(mean(nodesOut)*100, 2), "% missing)\n",
        paste(nodes[nodesOut], collapse = ", ")
      ))
      nodes <- nodes[!nodesOut]
    }
  }
  
  # Filter by node
  if(!is.null(nodes)) {
    cat("Filtering for", length(nodes), "nodes\n")
    dgetab <- dgetab[dgetab$Node %in% nodes,]
  }
  
  # Filter by p-values (if maxFDR = NULL)
  if(is.null(maxFDR)) {
    if(!is.null(maxPval)) {
      if(!is.numeric(maxPval) | maxPval < 0) {
        stop("'maxPval' must be numeric value > 0")
      } else {
        cat("Filtering for P-value <", maxPval, "\n")
        dgetab <- dgetab[dgetab$`P Value` < maxPval,]
      }
    }
  }
  
  # Filter by FDR (if maxFDR == NULL)
  if(!is.null(maxFDR)) {
    if(!is.numeric(maxFDR) | maxFDR < 0) {
      stop("'maxFDR' must be numeric value > 0")
    } else {
      cat("Filtering for FDR <", maxFDR, "\n")
      dgetab <- dgetab[dgetab$FDR < maxFDR,]
    }
  }
  
  # Filter by minumum difference
  if(!is.null(minDiff)) {
    if(!is.numeric(minDiff) | minDiff < 0) {
      stop("'minDiff' must be numeric value > 0")
    } else {
      cat("Filtering for |Differences (LogFC)| >", minDiff, "\n")
      dgetab <- dgetab[abs(dgetab$Diff) > minDiff,]
    }
  }
  
  # Check and print number of results
  if(nrow(dgetab) == 0) {
    stop("No results after filtering for input criteria")
  }
  cat("Generating table for", nrow(dgetab), "results\n")
  
  ## Create data table obect
  datatable(dgetab, rownames=FALSE, escape=FALSE,
            filter=list(position="top", clear=FALSE),
            options=list(columnDefs=list(list(className="dt-center", 
                                              targets="_all")),
                         search=list(regex=TRUE), 
                         scrollX=TRUE, 
                         scrollY=TRUE, 
                         dom="Brtp", 
                         paging=TRUE,
                         pageLength=pagelength,
            selection="none")) %>%
    formatRound(c("Mean", "Diff"), digits=2) %>%
    formatSignif(c("P Value", "FDR"), digits=2) %>%
    formatStyle(c("Direction", "Mean"), `border-right`="solid 2px")
}