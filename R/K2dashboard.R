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
