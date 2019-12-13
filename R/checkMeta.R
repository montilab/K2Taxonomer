# Function to change meta data if new value is entered.
.checkMeta <- function(K2res, arg, argValue) {
  
  if (is.null(argValue)) {
    argValue <- K2meta(K2res)[[arg]]
  }
  if (is.null(argValue)) {
    stop(cat("No value of ", arg, " specified.\n"))
  }
  
  return(argValue)
  
}
