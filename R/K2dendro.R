#' Reformat K2Taxonomer results to dendrogram object
#'
#' Reformats the output of K2tax() to an object of class, dendrogram,
#' which can be easily plotted.
#' @param K2res An object of class K2. The output of K2tax().
#' @return An object of class dendrogram.
#' @keywords clustering
#' @export
#' @import dendextend
#' @examples
#' K2dendro(K2res)

K2dendro <- function(K2res) {
    
    ## Get labels order
    K2labs <- colnames(K2data(K2res))
    
    ## Pull out K2 results
    K2res <- K2results(K2res)
    
    ## Get split identifiers
    splitNames = names(K2res)
    
    ## Create Matric of results
    mat <- matrix(0, nrow = length(K2res), ncol = length(K2res[[1]]$obs[[1]]) + length(K2res[[1]]$obs[[2]]))
    colnames(mat) <- c(K2res[[1]]$obs[[1]], K2res[[1]]$obs[[2]])
    for (i in 1:length(K2res)) {
        mat[i, K2res[[i]]$obs[[1]]] <- 1
        mat[i, K2res[[i]]$obs[[2]]] <- 2
    }
    
    ### Collapse mat
    matCollapse <- sort(apply(mat, 2, function(x) paste(x[x != 0], collapse = "")))
    matUnique <- unique(matCollapse)
    
    ## Get branchlist
    bList <- c()
    j <- 1
    for (i in 1:max(nchar(matUnique))) {
        sLength <- matUnique[nchar(matUnique) >= i]
        sLength <- unique(substr(sLength, 1, i))
        for (k in sLength) {
            bList[j] <- k
            j <- j + 1
        }
    }
    
    if (nrow(mat) > 1) {
        ## Add edgenames
        for (i in 2:nrow(mat)) {
            matNow <- mat[i, ]
            wn0 <- which(matNow != 0)[1]
            matPrev <- mat[1:(i - 1), wn0]
            matPrev <- paste(matPrev[matPrev != 0], collapse = "")
            names(splitNames)[i] <- matPrev
        }
    }
    
    ## Max height = number of observations
    mHeight <- ncol(mat)
    
    ## Add list
    aList <- list()
    newLabel <- splitNames[1]
    names(newLabel) <- NULL
    stability <- K2res[[splitNames[1]]]$stability$node
    log_stab_cum <- log(stability)
    attributes(aList) <- list(members = length(matCollapse), height = NULL, midpoint = (length(matCollapse) - 
        1)/2, label = newLabel, stability = stability, log_stab_cum = log_stab_cum)
    
    for (i in bList) {
        
        ## Add element
        iSplit <- unlist(strsplit(i, ""))
        iPaste <- paste0("aList", paste(paste0("[[", iSplit, "]]"), collapse = ""))
        eval(parse(text = paste0(iPaste, "<-list()")))
        
        ## Get parent node height
        pSplit <- iSplit[-length(iSplit)]
        if (length(pSplit) == 0) {
            pPaste <- "aList"
        } else {
            pPaste <- paste0("aList", paste(paste0("[[", pSplit, "]]"), collapse = ""))
        }
        pNode <- eval(parse(text = pPaste))
        
        ## Add attributes
        members <- sum(substr(matCollapse, 1, nchar(i)) == i)
        newLabel <- splitNames[i]
        names(newLabel) <- NULL
        
        ## Set height = 1 if terminal split, otherwise get stability and cumulative
        ## stability
        if (!any(substr(bList, 1, nchar(bList) - 1) == i)) {
            height <- 1
            stability <- log_stab_cum <- NULL
        } else {
            height = NULL
            stability <- K2res[[splitNames[i]]]$stability$node
            log_stab_cum = log(stability) + attributes(pNode)$log_stab_cum
        }
        att <- list(members = members, height = height, midpoint = (members - 1)/2, 
            label = newLabel, stability = stability, log_stab_cum = log_stab_cum, 
            index = iSplit)
        
        eval(parse(text = paste0("attributes(", iPaste, ") <- att")))
        
        ## Add leaves
        leaves <- matCollapse[matCollapse == i]
        if (length(leaves) > 0) {
            for (l in 1:length(leaves)) {
                
                ## Get attributes
                members <- 1
                leaf <- names(leaves)[l]
                names(leaf) <- NULL
                height <- 0
                
                ## Add element
                lPaste <- paste0(iPaste, "[[", l, "]]")
                eval(parse(text = paste0(lPaste, "<-", which(K2labs == leaf), "L")))
                
                att <- list(members = members, height = height, label = leaf, leaf = TRUE, 
                  index = c(iSplit, l))
                eval(parse(text = paste0("attributes(", lPaste, ") <- att")))
            }
        }
    }
    class(aList) <- "dendrogram"
    
    ## Get heights scaling factor
    logStabCums <- get_nodes_attr(aList, "log_stab_cum")
    members <- get_nodes_attr(aList, "members")
    minHeight <- min(log(members) + logStabCums, na.rm = t)
    
    ## Add heights to upper nodes
    aList <- dendrapply(aList, function(x, minHeight) {
        attr <- attributes(x)
        if (is.null(attr$height)) {
            attr$height <- log(attr$members) + attr$log_stab_cum - minHeight + 1
        }
        attributes(x) <- attr
        return(x)
    }, minHeight)
    
    return(aList)
}
