## Check that object is class K2
.isK2 <- function(K2res) {
    
    if (!is(K2res, "K2")) {
        stop("Argument, K2res is not a K2 class object.\n")
    }
    
}