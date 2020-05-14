## Function to format covariates string in formula
.formatCov <- function(covariates) {
    if (is.null(covariates)) 
        "" else paste0("+", paste(covariates, collapse = "+"))
}
