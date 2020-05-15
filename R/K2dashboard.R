#' Generate interactive dashboard of K2 Taxonomer results
#'
#' This function will generate an interactive dashboard of the annotated K=2
#' clustering results in a specified directory.
#' @param K2res Object of class K2.
#' @param analysis_name Character string of the name of analysis to write files
#' and generate title of report.
#' @param output_dir Base directory to put report.
#' @return Nothing.  Writes files to a specified directory.
#' @keywords clustering
#' @export
#' @examples
#' K2dashboard(K2results, analysis_name = 'K2Taxonomer', output_dir = '.')
#'

K2dashboard <- function(K2res, analysis_name = "K2Taxonomer", output_dir = ".") {
    
    ## Run checks
    
    ## Check K2 object
    k2Check <- .checkK2(K2res)
    
    ## K2 algorithm
    if (length(K2results(K2res)) == 0) {
        "No results found. Please run K2tax() or runK2Taxonomer().\n"
    }
    
    ## DGE
    if (is.null(K2results(K2res)[[1]]$dge)) {
        "No differential analysis results found. Please run runDGEmods().\n"
    }
    
    ## GSE
    if (is.null(K2results(K2res)[[1]]$gse)) {
        "No enrichment results found. Please run runDGEmods().\n"
    }
    
    ## GSVA
    if (ncol(K2gSet(K2res)) == 0) {
        "No ssGSEA data found. Please run runGSVAmods().\n"
    }
    
    ## DSSE
    if (is.null(K2results(K2res)[[1]]$dsse)) {
        "No differential enrichment results found. Please run runDSSEmods().\n"
    }

    ## DSSE
    if (is.null(K2results(K2res)[[1]]$dsse)) {
        "No differential enrichment results found. Please run runDSSEmods().\n"
    }
    
    ## Create file paths
    analysis_name_nospace <- gsub(" ", "_", analysis_name)

    ## Directory to print
    dirPath <- file.path(output_dir, paste(analysis_name_nospace, gsub("-| |:", "_",
        Sys.time()), sample(1e+07, 1), sep = "_"))
    ## RMD file
    RMDpath <- file.path(dirPath, paste0(analysis_name_nospace, ".Rmd"))

    ## Create folder to put report and results file
    dir.create(dirPath)

    ## Save K2results to this folder
    saveRDS(K2results, file.path(dirPath, "K2results.rds"))

    ## Get rmd file location
    K2rmd <- system.file("dashboard", "K2Dashboard.Rmd", package = "K2Taxonomer")

    ## Write RMD file
    K2rmdLines <- readLines(K2rmd)
    K2rmdLines[2] <- gsub("K2 Taxonomic Clustering", analysis_name, K2rmdLines[2])
    writeLines(K2rmdLines, RMDpath)

    ## Print success message
    cat("Dashboard created in '", dirPath, "'.\n", sep = "")

    invisible(NULL)
}
