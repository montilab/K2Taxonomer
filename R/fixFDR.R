.fixFDR <- function(K2res, analysis) {
    
    ## Create data.frame of pvalues and calculate FDR
    pValueDF <- data.frame(pval = unlist(lapply(K2results(K2res), function(x) {
        x[[analysis]]$pval
    })))
    pValueDF$fdr <- p.adjust(pValueDF$pval, method = "BH")
    pValueDF <- unique(pValueDF)
    
    ## Merge with original results
    K2results(K2res) <- lapply(K2results(K2res), function(x, pValueDF, analysis) {
        
        ## Remove original FDR
        xSub <- x[[analysis]][, colnames(x[[analysis]]) != "fdr"]
        
        ## Merge new fdr
        xSub <- merge(xSub, pValueDF)
        
        ## Sort by p-value
        xSub <- xSub[order(xSub$pval), ]
        
        ## Reorder columns
        if(nrow(xSub) > 0) {
            if (analysis == "dge") {
                xSub <- xSub[, c("gene", "coef", "mean", "t", "pval", "fdr", "B", "edge")]
            } else {
                xSub <- xSub[, c("category", "coef", "mean", "t", "pval", "fdr", "B", 
                    "edge")]
            }
            ## Add back to K2results()
            x[[analysis]] <- xSub
        }
    
        return(x)
    }, pValueDF, analysis)
    
    return(K2res)
}
