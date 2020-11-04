#' Extract table of phenotypic variable tests from 'K2' object
#'
#' Create table phenotypic variable results from 'K2' object.
#' @param K2res An object of class K2 or K2results().
#' @return A data.frame object with the following columns: 
#' \itemize{
#'  \item{value: }{The variable being tested}
#'  \item{node: }{The partition label}
#'  \item{edge: }{The subgroup for the given partition}
#'  \item{pval: }{Nominal p-value of test}
#'  \item{fdr: }{Benjamini-Hochberg FDR corrected p-value}
#'  \item{df: }{Degrees of freedom of test}
#'  \item{stat: }{Test statistic}
#'  \item{obsMean: }{Mean value across partition members}
#'  \item{altMean: }{Mean value for all other observations}
#'  \item{diffMean: }{Difference is mean}
#'  \item{nhits: }{The number of second label values in subgroup}
#'  \item{ncase: }{The total second-level label value}
#'  \item{nalt: }{The total first-level label value}
#'  \item{ndrawn: }{The total members in the subgroup}
#'  \item{hits: }{Members of subgroup with second-level label value}
#'  \item{test: }{The statistical test that produced this result}
#' }
#' @keywords clustering
#' @export
#' @examples
#' ## Read in ExpressionSet object
#' library(Biobase)
#' data(sample.ExpressionSet)
#' 
#' ## Pre-process and create K2 object
#' K2res <- K2preproc(sample.ExpressionSet)
#' 
#' ## Run K2 Taxonomer algorithm
#' K2res <- K2tax(K2res,
#'                stabThresh = 0.5)
#' 
#' K2res <- infoClassVector <- c(
#' sex = "factor",
#' score = "numeric1"
#' )
#'
#' runTestsMods <- function(K2res, infoClass = infoClassVector)
#'
#' head(getTestsModTable(K2res))
#' 

getTestsModTable <- function(K2res) {
    
    # Run checks
    .isK2(K2res)
    
    K2resList <- K2results(K2res)
    
    ## Format table
    K2modTestTable <- do.call(rbind, lapply(names(K2resList), function(x) {
        
        # Get phenotype results
        modTests <- K2resList[[x]]$modTests
        
        # Add edge ID to each subgroup
        modTests[[1]]$edge <- "1"
        modTests[[2]]$edge <- "2"
        
        # Concatenate
        modTests <- do.call(rbind, modTests)
        
        # Add node ID
        modTests$node <- x
        
        return(modTests)
    }))
    rownames(K2modTestTable) <- NULL
    
    ## Sort columns
    K2modTestTable <- K2modTestTable[, c("value", "node", "edge", "pval", "fdr", 
                                         "stat", "df", "obsMean", "altMean", 
                                         "diffMean", "nhits", "ncase", 
                                         "nalt", "ndrawn", "hits", "test")]
    
    ## Sort by p-value
    K2modTestTable <- K2modTestTable[order(K2modTestTable$pval),]
    
    return(K2modTestTable)
}
