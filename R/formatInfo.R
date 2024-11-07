.formatInfo <- function(K2r) {

  ## Format info summary
  if (is.null(K2meta(K2r)$cohorts)) {
    if(!is.null(K2meta(K2r)$info)) {
      infoDF <- data.frame(ID=colnames(K2data(K2r)),
                           K2colData(K2r)[colnames(K2data(K2r)), K2meta(K2r)$info],
                           row.names=colnames(K2data(K2r)),
                           stringsAsFactors=FALSE)
    } else {
      infoDF <- data.frame(ID=colnames(K2data(K2r)),
                           row.names=colnames(K2data(K2r)), 
                           stringsAsFactors=FALSE)
    }
      
  } else {
  
    tabDF <- as.data.frame(table(K2colData(K2r)[,K2meta(K2r)$cohort]))
    infoDF <- data.frame(ID = as.factor(tabDF$Var1),
                         N = tabDF$Freq,
                         row.names = as.character(tabDF$Var1))
    
    ## Format info summary
    if (!is.null(K2meta(K2r)$info)) {
      cdSub <- K2colData(K2r)[, c(K2meta(K2r)$cohort, K2meta(K2r)$info)]
      colnames(cdSub)[1] <- "ID"
      cdClasses <- sapply(cdSub,class)[-1]
      
      cdSum <- do.call(cbind, lapply(names(cdClasses), function(nam) {
        cl <- cdClasses[nam]
        cld <- cdSub[, c("ID", nam)]
        if(cl %in% c("character", "factor")) {
          cltab <- table(cld)
          sv <- as.factor(apply(cltab, 1, function(x) {
            paste(colnames(cltab)[x > 0], collapse = ";")
          })[rownames(infoDF)])
          sv <- data.frame(sv)
          colnames(sv) <- paste0(nam, "_labs")
        } else {
          sv <- sapply(rownames(infoDF), function(x) {
            median(cld[cld$ID == x, 2])
          })
          sv <- data.frame(sv)
          colnames(sv) <- paste0(nam, "_med")
        }
        return(data.frame(sv))
      }))
      
      infoDF <- cbind(infoDF, cdSum)
    }
  }
  
  return(infoDF)
  
}
