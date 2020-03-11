#' Function K2 Taxonomer
#'
#' This function performs K2 Taxonomer procedure only. Arguments values are extracted from K2meta(K2res) unless othewise specified.
#' 
#' @param K2
#' @param nFeats A numeric value <= P of subsets of the data to use.
#' @param featMatric Metric to use to assign variance/signal score. Options are
#' "square" (default), "mad" to use MAD scores, "sd" to use standard deviation
#' @param nBoots A numeric value of the number of bootstraps to run at each split.
#' @param clustFunc Wrapper function to cluster a P x N (See details).
#' @param clustCors Number of cores to use for clustering.
#' @param clustList List of objects to use for clustering procedure.
#' @param linkage Linkage criteria for splitting cosine matrix ("method" in hclust).
#' @param oneoff Logical. Allow 1 member clusters?
#' @param stabThresh A numeric value < 1, to set stopping threshold (use any negative value for no threshold).
#' @return An object of class K2.
#' @keywords clustering
#' @export
#' @import parallel
#' @examples
#' K2tax(K2res)
#'

K2tax <- function(K2res,
                    nFeats = NULL,
                    featMetric = NULL,
                    nBoots = NULL,
                    clustFunc =  NULL,
                    clustCors = NULL,
                    clustList = NULL,
                    linkage = NULL,
                    oneoff = NULL,
                    stabThresh = NULL){
  
  # Change meta data if new value is specific
  K2meta(K2res)$nFeats <- nFeats <- .checkMeta(K2res, "nFeats", nFeats)
  K2meta(K2res)$featMetric <- featMetric <- .checkMeta(K2res, "featMetric", featMetric)
  K2meta(K2res)$nBoots <- nBoots <- .checkMeta(K2res, "nBoots", nBoots)
  K2meta(K2res)$clustFunc <- clustFunc <- .checkMeta(K2res, "clustFunc", clustFunc)
  K2meta(K2res)$clustCors <- clustCors <- .checkMeta(K2res, "clustCors", clustCors)
  K2meta(K2res)$clustList <- clustList <- .checkMeta(K2res, "clustList", clustList)
  K2meta(K2res)$linkage <- linkage <- .checkMeta(K2res, "linkage", linkage)
  K2meta(K2res)$oneoff <- oneoff <- .checkMeta(K2res, "oneoff", oneoff)
  K2meta(K2res)$stabThresh <- stabThresh <- .checkMeta(K2res, "stabThresh", stabThresh)
  
  # Create Splits of the data
  taxList <- list(list(colnames(K2data(K2res))))
  resList <- stabList <- list(list(NULL))
  iter <- 1
  while(max(unlist(lapply(taxList[[iter]], 
                          function(x) length(x[x!="Vehicle"])
  ))) > (2 - oneoff)){
    outList <- lapply(taxList[[iter]], function(samps, nFeats, nBoots, clustFunc){
      if( length(samps) > 1 ){
        Dsub <- K2data(K2res)[,samps]
        outList <- .do_split(Dsub,
                             nFeats = nFeats,
                             featMetric = featMetric,
                             nBoots = nBoots,
                             clustFunc = clustFunc,
                             clustCors = clustCors,
                             clustList = clustList,
                             linkage = linkage)
        
        # Get minimum size
        minSize <- min(table(outList$mods))
        
        # Get determinant
        nodeStab <- outList$stability$node
        
        # Check stopping criteria 
        if ( (!oneoff & minSize == 1 & iter != 1) | 
             (nodeStab <= stabThresh & iter != 1) ) {
          outList <- list(mods = samps, propBoots = NULL, clustStab = NULL)
        }
        
        # If first split was below threshold, print warning.
        if(nodeStab <= stabThresh & iter == 1) {
          warning("First split below stability threshold.")
        }
        
      } else {
        outList <- list(mods = samps, propBoots = NULL, clustStab = NULL)
      }
      return(outList)
    }, nFeats, nBoots, clustFunc)
    iter <- iter + 1
    taxList[[iter]] <- list()
    outMods <- lapply(outList, function(x) x[[1]])
    slot <- 1
    for(i in 1:length(outMods)){
      taxList[[iter]][[slot]] <- names(outMods[[i]])[outMods[[i]] == 1]; slot <- slot + 1
      taxList[[iter]][[slot]] <- names(outMods[[i]])[outMods[[i]] == 2]; slot <- slot + 1
    }
    
    if(length(taxList[[iter]])==0){
      taxList[[iter]] <- list(NULL)
    }
    
    resList[[iter]] <- list()
    stabList[[iter]] <- list()
    outRes <- lapply(outList, function(x) x[[2]])
    outStab <- lapply(outList, function(x) x[[3]])
    slot <- 1
    for(i in 1:length(outRes)){
      resList[[iter]][[slot]] <- outRes[[i]]
      stabList[[iter]][[slot]] <- outStab[[i]]
      slot <- slot + 1
    }
  }
  resList <- resList[-1]
  stabList <- stabList[-1]
  
  if(is.null(unlist(taxList[[iter]]))) taxList <- taxList[-iter]
  
  # Get instances where the split had > 1 samples in cluster
  modList <- lapply(taxList[-1], function(x){
    combs <- 1:(length(x)/2)
    combList <- list()
    for(i in combs){
      combList[[i]] <- list()
      combList[[i]][[1]] <- x[[i*2-1]]
      combList[[i]][[2]] <- x[[i*2]]
      if(length(combList[[i]]) > 0){
        if(length(combList[[i]][[1]]) < (2 - oneoff) | length(combList[[i]][[2]]) < (2 - oneoff)) combList[[i]] <- list()
      }
    }
    combList
  })
  
  # Create list of modules and bootstrap stats
  K2c <- lapply(1:length(modList), function(x){
    lapply(1:length(modList[[x]]), function(y){
      modSub <- modList[[x]][[y]]
      if(length(modSub) > 0){
        resSub <- resList[[x]][[y]]
        stabSub <- stabList[[x]][[y]]
        return(list(obs = modSub, bootP = resSub, stability = stabSub))
      } else {
        return(NULL)
      }
    })
  })
  K2c <- unlist(K2c, recursive = F)
  K2c <- K2c[!unlist(lapply(K2c, is.null))]
  
  # Add names
  len <- length(K2c)
  nAlphabets <- ceiling(len/length(LETTERS))
  ALPHABETS <- unlist(lapply(1:nAlphabets, function(x) 
    sapply(LETTERS, function(y) paste(rep(y, x), collapse = ""))))
  names(K2c) <- ALPHABETS[1:length(K2c)]
  
  # Add to K2res and return
  K2results(K2res) <- K2c
  return(K2res)
  
}