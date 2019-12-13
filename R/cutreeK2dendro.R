#' Cut K2dendro output based on specified number of clusters
#'
#' Returns cluster for each observation based on input k value
#' @param K2dendrogram The output of K2dendr
#' @param k Specified number of clusters
#' @return A vector of cluster assignments
#' @keywords clustering
#' @export
#' @import dendextend
#' @examples
#' cutreeK2dendro(K2dendrogram, k = 2)


cutreeK2dendro <- function(K2dendrogram, k = 2) {
  
  # Get node heights
  k2Heights <- get_nodes_attr(K2dendrogram, "height")
  k2Index <- get_nodes_attr(K2dendrogram, "index")
  
  # Assign rank
  k2Order <- order(k2Heights, decreasing = TRUE)
  
  # Start at first split
  leaves <- get_leaves_attr(K2dendrogram, "label")
  child1 <- get_leaves_attr(K2dendrogram[[1]], "label")
  clusters <- as.numeric(leaves %in% child1) + 1
  
  # Initialize i
  i <- 2
  
  while (i < k) {
    
    # Get parent node to split
    splitIndex <- k2Index[[k2Order[i]]]
    dendSplit <- eval(parse(
      text = paste0("K2dendrogram[[", paste(splitIndex, collapse = "]][["), "]]")
    ))
    
    # Get new cluster assignment
    i <- i + 1
    
    # Get children
    child1 <- get_leaves_attr(dendSplit[[1]], "label")
    clusters[leaves %in% child1] <- i
    
  }
  
  # Set order from left to right
  clustOrds <- unique(clusters)
  clusters <- as.character(as.numeric(factor(clusters, levels = clustOrds)))
  
  # Add leaves
  names(clusters) <- leaves
  return(clusters)
  
}