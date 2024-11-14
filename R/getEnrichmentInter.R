#' Create interactive table of gene set enrichment results from 'K2' object
#'
#' Create table gene set enrichment results from 'K2' object.
#' @param maxFDR_score Numeric. A value between 0 and 1 indicating the FDR cutoff for differential analysis of enrichment scores.
#' @param minDiff_score Numeric. A value between 0 and 1 indicating the mean difference cutoff for differential analysis of enrichment scores.
#' @param maxPval_score Numeric. A value between 0 and 1 indicating the p-value cutoff for differential analysis of enrichment scores.
#' @param maxFDR_fisher Numeric. A value between 0 and 1 indicating the FDR cutoff for Fisher-based overrepresentation analysis.
#' @param maxPval_fisher Numeric. A value between 0 and 1 indicating the p-value cutoff for Fisher-based overrepresentation analysis.
#' @param gsNames Character. A vector of gene sets identifiers to display.
#' @param nodes Character. A vector of node identifiers to display.
#' @return An interactive data frame with the following columns:
#' \itemize{
#'  \item{Gene Set: }{The identifier of the gene set.}
#'  \item{Node: }{The identifier of the partition}
#'  \item{Edge: }{The identifier of partition subgroup}
#'  \item{Direction: }{The direction of coefficient for the assigned gene sets}
#'  \item{P Value Fisher: }{The p-value estimated by Fisher-based overrepresentation analysis}
#'  \item{FDR Fisher: }{The multiple hypothesis corrected FDR Fisher p-value, adjusted across
#'  all partitions}
#'  \item{N Overlap: }{The intersection of differentially expressed genes and gene set}
#'  \item{N Sig. Genes: }{The number differentially expressed genes}
#'  \item{N Gene Set: }{The number genes comprising the gene set}
#'  \item{P Value Score: }{The p-value estimated by differential analysis of enrichment scores}
#'  \item{FDR Score: }{The multiple hypothesis corrected FDR p-value estimated by differential analysis of enrichment scores, adjusted across
#'  all partitions}
#'  \item{Diff Score: }{The difference between the mean enrichment score of each subgroup at
#'  a given partition}
#' }
#' @inheritParams K2tax
#' @export

getEnrichmentInter <- function(K2res, 
                               maxFDR_score = 0.01,
                               minDiff_score = NULL,
                               maxPval_score = NULL,
                               maxFDR_fisher = NULL,
                               maxPval_fisher = NULL,
                               gsNames = NULL,
                               nodes = NULL) {
  
  # Get results tables
  enrtab <- getEnrichmentTable(K2res)
  enrtab <- enrtab[, !colnames(enrtab) %in% c("ntot", "t", "hits", "mean")]
  colnames(enrtab) <- c("Gene Set", "Node", "Edge", "Direction",
                        "P Value_Fisher", "FDR_Fisher", "N_Overlap", 
                        "N_Sig. Genes", "N_Gene Set", "P Value_Score", 
                        "FDR_Score", "Diff_Score")
  
  # Handle inputs
  
  ## Gene sets
  gsNamesOut <- !gsNames %in% enrtab$`Gene Set`
  if(sum(gsNamesOut) > 0) {
    if(mean(gsNamesOut) == 1) {
      stop("0 input gene set names found in results table")
    } else {
      warning(paste0(
        sum(gsNamesOut), " input gene set(s) not found in results table (", signif(mean(gsNamesOut)*100, 2), "% missing)\n",
        paste(gsNames[gsNamesOut], collapse = ", ")
      ))
      gsNames <- gsNames[!gsNamesOut]
    }
  }
  
  # Filter by gene input
  if(!is.null(gsNames)) {
    cat("Filtering for", length(gsNames), "gene sets\n")
    enrtab <- enrtab[enrtab$`Gene Set` %in% gsNames,]
  }
  
  ## Nodes
  nodesOut <- !nodes %in% enrtab$Node
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
    enrtab <- enrtab[enrtab$Node %in% nodes,]
  }
  
  # Filter by Enr. Score p-values (if maxFDR_score = NULL)
  if(is.null(maxFDR_score)) {
    if(!is.null(maxPval_score)) {
      if(!is.numeric(maxPval_score) | maxPval_score < 0) {
        stop("'maxPval_score' must be numeric value > 0")
      } else {
        cat("Filtering for Enr. Score P-value <", maxPval_score, "\n")
        enrtab <- enrtab[!is.na(enrtab$`P Value_Score`) & enrtab$`P Value_Score` < maxPval_score,]
      }
    }
  }
  
  # Filter by Enr. Score FDR (if maxFDR_score == NULL)
  if(!is.null(maxFDR_score)) {
    if(!is.numeric(maxFDR_score) | maxFDR_score < 0) {
      stop("'maxFDR_score' must be numeric value > 0")
    } else {
      cat("Filtering for Enr. Score FDR <", maxFDR_score, "\n")
      enrtab <- enrtab[!is.na(enrtab$FDR_Score) & enrtab$FDR_Score < maxFDR_score,]
    }
  }
  
  # Filter by Enr. Score minumum difference
  if(!is.null(minDiff_score)) {
    if(!is.numeric(minDiff_score) | minDiff_score < 0) {
      stop("'minDiff_score' must be numeric value > 0")
    } else {
      cat("Filtering for Enr. Score |Differences| >", minDiff_score, "\n")
      enrtab <- enrtab[!is.na(enrtab$Diff_Score) & abs(enrtab$Diff_Score) > minDiff_score,]
    }
  }
  
  # Filter by fisher p-values (if maxFDR_fisher = NULL)
  if(is.null(maxFDR_fisher)) {
    if(!is.null(maxPval_fisher)) {
      if(!is.numeric(maxPval_fisher) | maxPval_fisher < 0) {
        stop("'maxPval_fisher' must be numeric value > 0")
      } else {
        cat("Filtering for Fisher P-value <", maxPval_fisher, "\n")
        enrtab <- enrtab[enrtab$`P Value_Fisher` < maxPval_fisher,]
      }
    }
  }
  
  # Filter by fisher FDR (if maxFDR_fisher == NULL)
  if(!is.null(maxFDR_fisher)) {
    if(!is.numeric(maxFDR_fisher) | maxFDR_fisher < 0) {
      stop("'maxFDR_fisher' must be numeric value > 0")
    } else {
      cat("Filtering for Fisher FDR <", maxFDR_fisher, "\n")
      enrtab <- enrtab[enrtab$FDR_Fisher < maxFDR_fisher,]
    }
  }
  
  # Check and print number of results
  if(nrow(enrtab) == 0) {
    stop("No results after filtering for input criteria")
  }
  cat("Generating table for", nrow(enrtab), "results\n")
  
  ## Add line breaks
  colnames(enrtab) <- gsub("_", "<br>", colnames(enrtab))
  
  ## Create data table obect
  datatable(enrtab, rownames=FALSE, escape=FALSE,
            filter=list(position="top", clear=FALSE),
            options=list(columnDefs=list(list(className="dt-center", 
                                              targets="_all")),
                         search=list(regex=TRUE), 
                         scrollX=TRUE, 
                         scrollY=TRUE, 
                         dom="Brtp", 
                         paging=TRUE,
                         pageLength=50,
            selection="none")) %>%
    formatRound(c("Diff<br>Score"), digits=2) %>%
    formatSignif(
      c("P Value<br>Fisher", "FDR<br>Fisher", "P Value<br>Score",
        "FDR<br>Score"), digits=2) %>%
    formatStyle(c("Direction", "N<br>Gene Set", "Diff<br>Score"),
                `border-right`="solid 2px")
}