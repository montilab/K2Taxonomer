#' Interactive K2 dendrogram
#'
#' Create an interactive dendrogram of the K2 Taxonomer results
#' @param labelsize Size of labels displayed in window.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#' @return An interactive dendrogram created by visNetwork::visNetwork().
#' @inheritParams K2tax
#' @export

K2visNetwork <- function(K2res, labelsize = 50) {

    ## Run checks
    .isK2(K2res)

    ## K2 algorithm
    if (length(K2results(K2res)) == 0) {
        stop("No results found. Please run K2tax() or runK2Taxonomer().\n")
    }

    ## Get results list
    K2r <- K2results(K2res)

    ## Generate Matrix from visNetwork
    mat <- matrix(0, nrow=length(K2r), ncol=length(K2r[[1]]$obs[[1]]) +
        length(K2r[[1]]$obs[[2]]))
    colnames(mat) <- c(K2r[[1]]$obs[[1]], K2r[[1]]$obs[[2]])
    for (i in seq_len(length(K2r))) {
        mat[i, K2r[[i]]$obs[[1]]] <- 1
        mat[i, K2r[[i]]$obs[[2]]] <- 2
    }
    rownames(mat) <- names(K2r)

    ## Calculate sizes
    sizes <- apply(mat, 1, function(x) sum(x != 0))

    ## Add Labels
    titles <- unlist(lapply(names(K2r), function(nam) {
      x <- K2r[[nam]]
      paste(
        "Node:", nam, "<br>",
        "Probability:", signif(x$bootP, 2), "<br>", 
        "Members(Edge:1):", length(x$obs[[1]]), "<br>", 
        "Members(Edge:2):", length(x$obs[[2]]), "<br>",
        "Stability(Edge:1):", signif(x$stability$clusters[[1]], 2), "<br>", 
        "Stability(Edge:2):", signif(x$stability$clusters[[2]], 2))
    }))
    names(titles) <- names(K2r)

    ## initialize leafe names
    len <- length(K2r)
    nalphabets <- ceiling(ncol(mat)/length(letters))
    nAlphabets <- ceiling(len/length(letters))
    alphabets <- unlist(lapply(seq_len(nalphabets), function(x) {

        vapply(letters, function(y) paste(rep(y, x), collapse=""),
            character(1))

    }))
    ALPHABETS <- unlist(lapply(seq_len(nAlphabets), function(x) vapply(LETTERS,
        function(y) paste(rep(y, x), collapse=""), character(1))))
    source <- c()
    target <- c()
    k <- 1

    ## Get edges
    for (i in seq_len(nrow(mat))) {

        source <- c(source, rep(rownames(mat)[i], 2))
        matRow <- mat[i, ]
        sub1 <- which(matRow == 1)[1]
        sub2 <- which(matRow == 2)[1]
        matSub1 <- mat[-seq_len(i), sub1]
        names(matSub1) <- rownames(mat)[-seq_len(i)]
        matSub2 <- mat[-seq_len(i), sub2]
        names(matSub2) <- rownames(mat)[-seq_len(i)]

        target1 <- names(matSub1)[which(matSub1 != 0)[1]]
        if (is.na(target1)) {
            target1 <- alphabets[k]
            sizes[target1] <- sum(matRow == 1)
            titles[target1] <- paste(colnames(mat)[matRow ==
                1], collapse="<br>")
            k <- k + 1
        }
        target2 <- names(matSub1)[which(matSub2 != 0)[1]]
        if (is.na(target2)) {
            target2 <- alphabets[k]
            sizes[target2] <- sum(matRow == 2)
            titles[target2] <- paste(colnames(mat)[matRow ==
                2], collapse="<br>")
            k <- k + 1
        }
        target <- c(target, target1, target2)
    }

    ## Set terminal nodes to 0 sizes
    sizes[names(sizes) %in% alphabets] <- 0

    ## Add Labels
    labs <- titles
    labs <- gsub("<br>", "\n", labs)
    labs[names(labs) %in% ALPHABETS] <- names(labs)[names(labs) %in%
        ALPHABETS]

    ## Add shapes
    shapes <- rep("diamond", length(sizes))
    names(shapes) <- names(sizes)
    shapes[names(shapes) %in% alphabets] <- "box"

    ## Get terminal node level
    levs <- vapply(2:nrow(mat), function(x) {
        matSub <- mat[seq_len(x), ]
        matSub <- matSub[, matSub[x, ] != 0]
        matSub <- matSub[, 1]
        sum(matSub != 0)
    }, numeric(1))
    levs <- c(1, levs)
    levs <- c(levs, rep(max(levs) + 1, length(sizes) - length(levs)))
    names(levs) <- names(sizes)

    ## Generate network
    nodes <- data.frame(id=names(levs), title=titles, level=levs,
        label=labs, shape=shapes, font.size=labelsize, 
        stringsAsFactors=FALSE)
    edges <- data.frame(from=source, to=target, stringsAsFactors=FALSE)

    ## Generate plot
    p <- visNetwork(nodes, edges) %>% visEdges(arrows="to") %>%
        visHierarchicalLayout(direction="LR") %>%
        visInteraction(zoomSpeed = 0.1)

    return(p)
}
