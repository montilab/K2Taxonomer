.fixFDR <- function(K2r, a) {
  
    K2names <- names(K2results(K2r))
    
    if(a != "gse") {
  
      ## Concatenate resultsand calculate FDR
      Dconc <- do.call(rbind, lapply(K2names, function(i) {
        resi <- K2results(K2r)[[i]][[a]]
        if(!is.null(resi)) {
          resi$i <- i
        }
        return(resi)
      }))
      Dconc$fdr <- p.adjust(Dconc$pval)
      
      K2results(K2r) <- lapply(K2names, function(i) {
        resi <- Dconc[Dconc$i == i, -ncol(Dconc)]
        K2ri <- K2results(K2r)[[i]]
        if(nrow(resi) > 0) {
          K2ri[[a]] <- resi
        }
        return(K2ri)
      })
    
    } else {
      
      ## Concatenate resultsand calculate FDR
      resiNames <- names(K2results(K2r)[[1]][[a]])
      Dconc <- do.call(rbind, lapply(K2names, function(i) {
        resi <- K2results(K2r)[[i]][[a]]
        resiNames <- names(resi)
        resi <- do.call(rbind, lapply(resiNames, function(j) {
          resj <- resi[[j]]
          if(!is.null(resj)) {
            resj$j <- j
          }
          return(resj)
        }))
        if(!is.null(resi)) {
          resi$i <- i
        }
        return(resi)
      }))
      Dconc$fdr <- p.adjust(Dconc$pval)
      
      K2results(K2r) <- lapply(K2names, function(i) {
        resi <- Dconc[Dconc$i == i, -ncol(Dconc)]
        resi <- lapply(resiNames, function(j) {
          resj <- resi[resi$j == j, -ncol(resi)]
          if(nrow(resj) == 0) {
            resj <- NULL
          }
          return(resj)
        })
        names(resi) <- resiNames
        K2ri <- K2results(K2r)[[i]]
        K2ri[[a]] <- resi
        return(K2ri)
      })
    }
    
    names(K2results(K2r)) <- K2names
    
    return(K2r)
}
