#' Interactive K2 dendrogram
#'
#' Create an interactive dendrogram of the K2 Taxonomer results
#' @param K2res A list object. The output of runGSE_mods().
#' @keywords clustering
#' @export
#' @examples
#' K2visNetwork(K2res)
#' 
K2visNetwork <- function(K2res){
  
  # Get results list
  K2res <- K2results(K2res)
  
  # Generate Matrix from visNetwork
  mat <- matrix(0, nrow = length(K2res), ncol = length(K2res[[1]]$obs[[1]]) + length(K2res[[1]]$obs[[2]]))
  colnames(mat) <- c(K2res[[1]]$obs[[1]], K2res[[1]]$obs[[2]])
  for(i in 1:length(K2res)){
    mat[i,K2res[[i]]$obs[[1]]] <- 1
    mat[i,K2res[[i]]$obs[[2]]] <- 2
  }
  rownames(mat) <- names(K2res)
  
  # Calculate sizes
  sizes <- apply(mat, 1, function(x) sum(x!=0))
  
  # Add Labels
  titles <- unlist(lapply(K2res, function(x) paste("Probability:", signif(x$bootP, 2),
                                            "<br>", "Members(Group 1):", length(x$obs[[1]]),
                                            "<br>", "Members(Group 2):", length(x$obs[[2]]),
                                            "<br>", "Stability(Group 1):", signif(x$stability$clusters[[1]], 2),
                                            "<br>", "Stability(Group 2):", signif(x$stability$clusters[[2]], 2)
                                            )))
  names(titles) <- names(K2res)
  
  # initialize leafe names
  len <- length(K2res)
  nalphabets <- ceiling(ncol(mat)/length(letters))
  nAlphabets <- ceiling(len/length(letters))
  alphabets <- unlist(lapply(1:nalphabets, function(x) sapply(letters, 
                                                              function(y) paste(rep(y, x), collapse = ""))))
  ALPHABETS <- unlist(lapply(1:nAlphabets, function(x) sapply(LETTERS, 
                                                              function(y) paste(rep(y, x), collapse = ""))))
  
  source <- c()
  target <- c()
  k <- 1
  
  # Get edges
  for(i in 1:nrow(mat)){
    
    source <- c(source, rep(rownames(mat)[i], 2))
    matRow <- mat[i,]
    sub1 <- which(matRow == 1)[1]
    sub2 <- which(matRow == 2)[1]
    matSub1 <- mat[-(1:i), sub1];names(matSub1) <- rownames(mat)[-(1:i)]
    matSub2 <- mat[-(1:i), sub2];names(matSub2) <- rownames(mat)[-(1:i)]
    
    target1 <- names(matSub1)[which(matSub1!=0)[1]]
    if(is.na(target1)){
      target1 <- alphabets[k]
      sizes[target1] <- sum(matRow == 1)
      titles[target1] <- paste(colnames(mat)[matRow == 1], collapse = "<br>")
      k <- k+1}
    target2 <- names(matSub1)[which(matSub2!=0)[1]]
    if(is.na(target2)){
      target2 <- alphabets[k]
      sizes[target2] <- sum(matRow == 2)
      titles[target2] <- paste(colnames(mat)[matRow == 2], collapse = "<br>")
      k <- k+1}
    target <- c(target, target1, target2)
  }
  
  # Set terminal nodes to 0 sizes
  sizes[names(sizes) %in% alphabets] <- 0
  
  # Add Labels
  labs <- titles
  labs <- gsub("<br>", "\n", labs)
  labs[names(labs) %in% ALPHABETS] <- names(labs)[names(labs) %in% ALPHABETS]
  
  # Add shapes
  shapes <- rep("diamond", length(sizes))
  names(shapes) <- names(sizes)
  shapes[names(shapes) %in% alphabets] <- "box"
  
  # Get terminal node level
  levs <- sapply(2:nrow(mat), function(x){
    matSub <- mat[1:x,]
    matSub <- matSub[,matSub[x,]!=0]
    matSub <- matSub[,1]
    sum(matSub!=0)
  })
  levs <- c(1, levs)
  levs <- c(levs, rep(max(levs) + 1, length(sizes) - length(levs)))
  names(levs) <- names(sizes)
  
  # Generate network
  nodes <- data.frame(id = names(levs), title = titles, level = levs, label = labs, shape = shapes, stringsAsFactors = FALSE)
  edges <- data.frame(from = source, to = target, stringsAsFactors = FALSE)
  
  # Generate plot
  p <- visNetwork(nodes, edges) %>%
    visEdges(arrows = "to") %>% 
    visHierarchicalLayout(direction = "LR")
  
  return(p)
}