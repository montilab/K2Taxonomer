## Function to format variables string in formula
.formatCov <- function(variables) {
    if (is.null(variables))
        "" else paste0("+", paste(variables, collapse="+"))
}
