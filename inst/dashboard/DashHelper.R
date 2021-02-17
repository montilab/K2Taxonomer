## Function to format covariates string in formula
.formatCov <- function(covariates) if (is.null(covariates)) "" else paste0("+", paste(covariates,
    collapse = "+"))

## Function to generate differential signature
.signatureWrapper <- function(eSet, cohorts, mods, vehicle = NULL, covariates = NULL,
    logCounts = FALSE) {

    ## Remove vehicle from mods and make a data frame
    mods <- mods[names(mods) != "Vehicle"]
    mods <- data.frame(mods = as.character(mods), GROUP = names(mods), stringsAsFactors = FALSE)

    modStats <- NULL
    if (length(unique(mods$mods)) > 1 && nrow(mods) > 2) {

        ## If replicates in data get unique cohorts
        if (is.null(cohorts)) {
            cohorts <- "GROUP"
            pData(eSet)[, cohorts] <- colnames(eSet)
        }

        ## Get unique groups
        gUnique <- unique(pData(eSet)[, cohorts])

        if (!is.null(vehicle)) {
            gUnique <- gUnique[gUnique != vehicle]
        }

        ## Subset data for mods
        eSub <- eSet[, pData(eSet)[, cohorts] %in% c(vehicle, mods$GROUP)]

        ## Drop levels
        pData(eSub) <- droplevels(pData(eSub))

        ## Create new variable in pData by merging with mods data.frame
        pData(eSub)$GROUP <- factor(pData(eSub)[, cohorts], levels = c(vehicle, mods$GROUP))
        pData(eSub)$Rownames <- rownames(pData(eSub))
        pD <- merge(pData(eSub), mods, all.x = TRUE)
        rownames(pD) <- pD$Rownames
        pD <- pD[colnames(eSub), , drop = FALSE]

        ## Add vehicle mod
        pD$mods[is.na(pD$mods)] <- "0"
        pD$mods <- as.factor(pD$mods)

        ## Add back to eSet
        pData(eSub) <- pD

        ## Need to run different analysis if there are cohorts or not
        if (!is.null(cohorts)) {

            ## Create design matrix
            design <- model.matrix(as.formula(paste0("~ 0 +", "mods", .formatCov(covariates))),
                data = pData(eSub))
            colnames(design) <- sub("mods", "X", colnames(design))

            ## Fit model
            fit <- lmFit(eSub, design)

            ## Fit contrasts
            modFit <- lapply(paste0("X", unique(mods$mods)), function(x, design,
                fit) {
                conString <- paste0(x, " - (", paste(colnames(design)[colnames(design) !=
                  x & colnames(design) %in% paste0("X", unique(mods$mods))], collapse = "+"),
                  ")/", sum(colnames(design) != x & colnames(design) %in% paste0("X",
                    unique(mods$mods))))
                contrasts <- makeContrasts(contrasts = conString, levels = design)
                contFit <- suppressWarnings(topTable(eBayes(contrasts.fit(fit, contrasts),
                  trend = logCounts, robust = logCounts), number = Inf, sort.by = "none"))
                contFit <- contFit[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
                  "B")]
                contFit$edge <- sub("X", "", x)
                return(contFit)
            }, design, fit)

        } else {

            design <- model.matrix(as.formula(paste0("~ 0 + ", "GROUP", .formatCov(covariates))),
                data = pData(eSub))
            colnames(design) <- sub("GROUP", "X", colnames(design))

            ## Fit model
            fit <- lmFit(eSub, design)

            ## Create contrasts strings
            modsFull <- unique(pD[, c("GROUP", "mods")])
            modsTable <- table(modsFull$mods)
            cVec <- vapply(names(modsTable), function(x) paste0("(", paste(paste0("X",
                modsFull$GROUP[modsFull$mods == x], collapse = "+")), ")/", sum(modsFull$mods ==
                x)), FUN.VALUE = integer(1))

            ## Run each contrast
            modFit <- lapply(as.character(unique(mods$mods)), function(x, cVec, design,
                fit) {
                conString <- paste0(cVec[x], " - (", paste(cVec[names(cVec) != x],
                  collapse = "+"), ")/", sum(names(cVec) != x))
                contrasts <- makeContrasts(contrasts = conString, levels = design)
                contFit <- topTable(eBayes(contrasts.fit(fit, contrasts), trend = logCounts,
                  robust = logCounts), number = Inf, sort.by = "none")
                contFit <- contFit[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
                  "B")]
                contFit$edge <- x
                return(contFit)
            }, cVec, design, fit)
        }

        ## Create vector of where to assign result
        if (is.null(vehicle)) {
            if (length(modFit) == 2) {
                one2 <- c(1, 2)[as.numeric(modFit[[1]]$t < 0) + 1]
            } else {
                one2 <- vapply(seq(nrow(modFit[[1]])), function(row, modFit) {
                  which.max(vapply(seq(length(modFit)), function(g, modFit, row) {
                    modFit[[g]][row, "t"]
                  }, modFit, row, FUN.VALUE = numeric(1)))
                }, modFit, FUN.VALUE = integer(1))
            }
        } else {
            one2 <- c(1, 2)[as.numeric(vapply(seq(nrow(modFit[[1]])), function(row,
                modFit) {
                modFit[[1]]$P.Value[row] > modFit[[2]]$P.Value[row]
            }, modFit, FUN.VALUE = logical(1))) + 1]
        }

        ## Create just one data.frame
        modStats <- as.data.frame(t(vapply(1:nrow(modFit[[1]]), function(row, one2,
            modFit) {
            unlist(as.numeric(modFit[[one2[row]]][row, ]))
        }, one2, modFit, FUN.VALUE = double(7))))
        colnames(modStats) <- colnames(modFit[[1]])
        rownames(modStats) <- rownames(modFit[[1]])

        ## Order by p-value
        modStats <- modStats[order(modStats$P.Value), ]
        modStats$adj.P.Val <- p.adjust(modStats$P.Value, method = "BH")

    }

    ## Save mods as character
    modStats$edge <- as.character(modStats$edge)

    ## Change column names
    colnames(modStats) <- c("coef", "mean", "t", "pval", "fdr", "B", "edge")

    ## Return
    return(modStats)
}

## Function to format differential results
geneTable <- function(DGETABLE, nodeID = NULL, geneList = NULL) {

    ## Get exact match for nodeID
    if (!is.null(nodeID) && nodeID == "") {
        nodeID <- NULL
    }
    if (!is.null(nodeID))
        nodeID <- paste0("^", nodeID, "$")

    ## Create data table obect
    datatable(DGETABLE, rownames = FALSE, extensions = "Buttons", escape = FALSE, filter = list(position = "top",
        clear = FALSE), options = list(columnDefs = list(list(searchable = FALSE,
        orderable = FALSE, width = "3px", targets = c(8, 9, 10)), list(className = "dt-center",
        targets = "_all")), search = list(regex = TRUE), searchCols = list(list(search = geneList),
        list(search = nodeID), NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL),
        scrollX = TRUE, scrollY = "325px", dom = "Brtp", paging = TRUE, pageLength = 50,
        buttons = list(list(extend = "collection", text = "Help", action = DT::JS("function ( e, dt, node, config ) {
                                  Shiny.setInputValue('geneHelp', true, {priority: 'event'});
                                  }")),
            list(extend = "collection", text = "Download All Results", action = DT::JS("function ( e, dt, node, config ) {
                                  Shiny.setInputValue('geneDL', true, {priority: 'event'});
                                  }")))),
        selection = "none") %>% formatRound(c("Mean", "Diff"), digits = 2) %>% formatSignif(c("P Value",
        "FDR"), digits = 2) %>% formatStyle(c("Direction", "Mean"), `border-right` = "solid 2px")
}

## Function to format hyperenrichment results
genesetTable <- function(ENRTABLE, nodeID = NULL, edgeID = NULL, dgeHits = NULL) {

    ## Get exact match for nodeID and edgeID
    if (!is.null(nodeID) && nodeID == "") {
        nodeID <- NULL
    }
    if (!is.null(edgeID) && edgeID == "") {
        edgeID <- NULL
    }
    if (!is.null(nodeID))
        nodeID <- paste0("^", nodeID, "$")

    ## Add line breaks
    colnames(ENRTABLE) <- gsub("_", "<br>", colnames(ENRTABLE))

    ## Create DT object
    outDT <- datatable(ENRTABLE, rownames = FALSE, extensions = "Buttons", escape = FALSE,
        filter = list(position = "top", clear = FALSE), options = list(columnDefs = list(list(searchable = FALSE,
            orderable = FALSE, width = "3px", targets = c(14, 15, 16)), list(visible = FALSE,
            targets = c(12, 13)), list(className = "dt-center", targets = "_all")),
            search = list(regex = TRUE), searchCols = list(list(search = dgeHits),
                list(search = nodeID), list(search = edgeID), NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL),
            scrollX = TRUE, scrollY = "325px", dom = "Brtp", paging = TRUE, pageLength = 50,
            buttons = list(list(extend = "collection", text = "Help", action = DT::JS("function ( e, dt, node, config ) {
                                           Shiny.setInputValue('hyperHelp', true, {priority: 'event'});
}")),
                list(extend = "collection", text = "Download All Results", action = DT::JS("function ( e, dt, node, config ) {
                                    Shiny.setInputValue('hyperDL', true, {priority: 'event'});
                                    }")))),
        selection = "none") %>% formatRound(c("Mean<br>ssGSEA", "Diff<br>ssGSEA"),
        digits = 2) %>% formatSignif(c("P Value<br>Hyper", "FDR<br>Hyper", "P Value<br>ssGSEA",
        "FDR<br>ssGSEA"), digits = 2) %>% formatStyle(c("Direction", "N<br>Gene Set",
        "Diff<br>ssGSEA"), `border-right` = "solid 2px")

}

## Function to plot gene expression
plotGenePathway <- function(eSet, gene, obs1, obs2, cohorts, vehicle, yaxis = "Expression") {

    if (gene %in% rownames(eSet)) {

        ## Format subgroup names
        if (is.null(cohorts)) {
            nams <- colnames(eSet)
        } else {
            nams <- pData(eSet)[, cohorts]
        }
        nams[nams == vehicle] <- "Vehicle"

        ## Create data.frame of expression values
        e <- Biobase::exprs(eSet)[gene, ]
        df <- data.frame(e = e, ch = nams, stringsAsFactors = FALSE)

        ## Subset for obs in groups
        df <- df[df$ch %in% c(obs1, obs2, "Vehicle"), ]
        df$edge <- "Edge:1"
        df$edge[df$ch %in% obs2] <- "Edge:2"
        df$edge[df$ch == "Vehicle"] <- "Vehicle"

        ## Get per Observation mean
        dfMeans <- df %>% group_by(ch) %>% summarise(me = mean(e))
        dfMeans$ch <- as.character(dfMeans$ch)
        dfMeans <- dfMeans[order(dfMeans$me, decreasing = TRUE), ]

        ## Sort levels by mean expression
        df$ch <- factor(df$ch, levels = dfMeans$ch)
        df <- merge(df, dfMeans)

        ## Add levels for boxplots
        df$edge2 <- df$edge

        ## Add rows for boxplots
        df2 <- df
        df2$ch <- df$edge
        df2$edge2 <- "Comparison"
        df2$e2 <- df2$e
        df2$e <- NA

        ## Concatenate
        df$e2 <- NA
        df <- df[df$ch != "Vehicle", ]
        df <- rbind(df, df2)

        ## Fix levels
        df$ch <- factor(df$ch, levels = c(dfMeans$ch[dfMeans$ch != "Vehicle"], "Edge:1",
            "Vehicle", "Edge:2"))
        df$edge <- factor(df$edge, levels = c("Edge:1", "Vehicle", "Edge:2"))
        df$edge2 <- factor(df$edge2, levels = c("Edge:1", "Comparison", "Edge:2"))

        ## Remove Means from comparison
        df$me[df$edge2 == "Comparison"] <- NA

        ## Add column names
        colnames(df) <- c("Observation", "Expression", "Subgroup", "Mean", "Subgroup2",
            "Expression2")

        ## Plot
        p <- ggplot(data = df, aes(x = Observation, y = Expression)) +
                geom_boxplot(aes(y = Expression2, fill = Subgroup)) +
                geom_line(aes(group = Observation)) +
                geom_point(aes(colour = Subgroup), size = 3) +
                geom_point(aes(y = Mean), shape = 3, size = 3) +
                facet_grid(~Subgroup2, scales = "free_x") +
                scale_colour_manual(values = c(`Edge:1` = "darkorange", Vehicle = "grey", `Edge:2` = "darkorchid1")) +
                scale_fill_manual(values = c(`Edge:1` = "darkorange", Vehicle = "grey", `Edge:2` = "darkorchid1")) +
                scale_x_discrete() +
                scale_y_continuous(name = yaxis) +
                theme_bw() + ggtitle(gene) +
                theme(
                    plot.margin = margin(0, 10, 0, 10),
                    legend.position = "none",
                    axis.text.x = element_text(angle = 45, hjust = 0, size = 15),
                    axis.text.y = element_text(size = 15),
                    axis.title.x = element_blank()
                    )

        p <- ggplotly(p)

        # Remove duplicates in data and hover text for individual observations
        p$x$data <- lapply(p$x$data, function(dat, cohorts) {
            if("text" %in% names(dat)) {

                uni <- !duplicated(dat$text) | is.na(dat$text)
                dat$text <- dat$text[uni]
                dat$x <- dat$x[uni]
                dat$y <- dat$y[uni]

                if(!is.null(cohorts)) {
                    dat$text <- gsub("Observation", "Group", dat$text)
                }

                if(sum(grepl("Mean:", dat$text) & !is.na(dat$text)) > 0) {
                    dat$text <- gsub("<br />Expression(.*)", "", dat$text)
                } else {
                    dat$text <- NULL
                }
            }
            return(dat)
        }, cohorts)

        ## Fix xaxis due to a bug in plotly
        whXaxis <- which(grepl("xaxis", names(p$x$layout)))
        for (i in whXaxis) {
            ticktext <- p$x$layout[[i]]$ticktext
            if (length(ticktext) == 1) {
                p$x$layout[[i]]$tickvals <- seq(2)
                p$x$layout[[i]]$ticktext <- c(ticktext, "")
            }
        }
        return(p)
    }

}

## Function to plot gene expression
plotGenePathwayClusters <- function(eSet, gene, groupList, cohorts, vehicle, yaxis = "Expression") {

    if (gene %in% rownames(eSet)) {

        ## Format group names
        if (is.null(cohorts)) {
            nams <- colnames(eSet)
        } else {
            nams <- pData(eSet)[, cohorts]
        }
        nams[nams == vehicle] <- "Vehicle"

        ## Create data.frame of expression values
        e <- Biobase::exprs(eSet)[gene, ]
        df <- data.frame(e = e, ch = nams, stringsAsFactors = FALSE)

        ## Get clusters
        obs <- unlist(groupList)

        ## Subset for obs in groups
        df <- df[df$ch %in% c(obs, "Vehicle"), ]
        df$edge <- "Vehicle"
        for (i in names(groupList)) {
            df$edge[df$ch %in% groupList[[i]]] <- i
        }

        ## Get per Observation mean
        dfMeans <- df %>% group_by(ch) %>% summarise(me = mean(e))
        dfMeans$ch <- as.character(dfMeans$ch)
        dfMeans <- dfMeans[order(dfMeans$me, decreasing = FALSE), ]

        ## Sort levels by mean expression
        df$ch <- factor(df$ch, levels = dfMeans$ch)
        df <- merge(df, dfMeans)

        ## Add levels for boxplots
        df$edge2 <- df$edge

        ## Add rows for boxplots
        df2 <- df
        df2$ch <- df$edge
        df2$edge2 <- "Comparison"
        df2$e2 <- df2$e
        df2$e <- NA

        ## Concatenate
        df$e2 <- NA
        df <- df[df$ch != "Vehicle", ]
        df <- rbind(df, df2)

        ## Fix levels
        df$ch <- factor(df$ch, levels = c(dfMeans$ch[dfMeans$ch != "Vehicle"], "Vehicle",
            names(groupList)))
        df$edge <- factor(df$edge, levels = c("Vehicle", names(groupList)))
        df$edge2 <- factor(df$edge2, levels = c(names(groupList), "Comparison"))

        ## Remove Means from comparison
        df$me[df$edge2 == "Comparison"] <- NA

        ## Add column names
        colnames(df) <- c("Observation", "Expression", "Group", "Mean", "Group2",
            "Expression2")

        ## Create color manual
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual", ]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        colMan <- c("grey", col_vector[seq(length(groupList))])
        names(colMan) <- c("Vehicle", names(groupList))

        ## Plot
        p <- ggplot(data = df, aes(x = Observation, y = Expression)) +
            geom_boxplot(aes(y = Expression2, fill = Group)) +
            geom_line(aes(group = Observation)) +
            geom_point(aes(colour = Group), size = 3) +
            geom_point(aes(y = Mean), shape = 3, size = 3) +
            facet_grid(~Group2, scales = "free_x") +
            scale_colour_manual(values = colMan) +
            scale_fill_manual(values = colMan) +
            scale_x_discrete() +
            scale_y_continuous(name = yaxis) +
            theme_bw() +
            ggtitle(gene) +
            theme(
                plot.margin = margin(0,10, 0, 10),
                legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 0, size = 15),
                axis.text.y = element_text(size = 15),
                axis.title.x = element_blank()
            )

        p <- ggplotly(p)

        # Remove duplicates in data and hover text for individual observations
        p$x$data <- lapply(p$x$data, function(dat, cohorts) {
            if("text" %in% names(dat)) {

                uni <- !duplicated(dat$text) | is.na(dat$text)
                dat$text <- dat$text[uni]
                dat$x <- dat$x[uni]
                dat$y <- dat$y[uni]

                if(!is.null(cohorts)) {
                    dat$text <- gsub("Observation", "Group", dat$text)
                }

                if(sum(grepl("Mean:", dat$text) & !is.na(dat$text)) > 0) {
                    dat$text <- gsub("<br />Expression(.*)", "", dat$text)
                } else {
                    dat$text <- NULL
                }
            }
            return(dat)
        }, cohorts)


        ## Fix xaxis due to a bug in plotly
        whXaxis <- which(grepl("xaxis", names(p$x$layout)))
        for (i in whXaxis) {
            ticktext <- p$x$layout[[i]]$ticktext
            if (length(ticktext) == 1) {
                p$x$layout[[i]]$tickvals <- seq(2)
                p$x$layout[[i]]$ticktext <- c(ticktext, "")
            }
        }
        return(p)
    }
}

## Function to format differential results
geneTableClusters <- function(clusterRes, subgroupID = NULL, geneList = NULL) {

    ## Get exact match for nodeID
    if (!is.null(subgroupID) && subgroupID == "") {
        subgroupID <- NULL
    }
    if (!is.null(subgroupID)) {
        nodeID <- paste0("^", subgroupID, "$")
    }

    ## Create data table obect
    datatable(clusterRes, rownames = FALSE, extensions = "Buttons", escape = FALSE, filter = list(position = "top",
        clear = FALSE), options = list(columnDefs = list(list(searchable = FALSE,
        orderable = FALSE, width = "3px", targets = c(6, 7, 8)), list(className = "dt-center",
        targets = "_all")), search = list(regex = TRUE), searchCols = list(list(search = geneList),
        list(search = subgroupID), NULL, NULL, NULL, NULL, NULL), scrollX = TRUE,
        scrollY = "325px", dom = "Brtp", paging = TRUE, pageLength = 50, buttons = list(list(extend = "collection",
            text = "Help", action = DT::JS("function ( e, dt, node, config ) {
                                  Shiny.setInputValue('geneHelp', true, {priority: 'event'});
                                  }")),
            list(extend = "collection", text = "Download All Results", action = DT::JS("function ( e, dt, node, config ) {
                                  Shiny.setInputValue('geneDLMulti', true, {priority: 'event'});
                                  }")))),
        selection = "none") %>% formatRound(c("Mean", "Diff"), digits = 2) %>% formatSignif(c("P Value",
        "FDR"), digits = 2) %>% formatStyle(c("Subgroup", "Mean"), `border-right` = "solid 2px")
}

## Function to format hyperenrichment results from multiple comparisons
genesetTableClusters <- function(ENRTABLE, subgroupID = NULL, dgeHits = NULL) {

    ## Get exact match for nodeID
    if (!is.null(subgroupID) && subgroupID == "") {
        subgroupID <- NULL
    }
    if (!is.null(subgroupID))
        nodeID <- paste0("^", subgroupID, "$")

    ## Add line breaks
    colnames(ENRTABLE) <- gsub("_", "<br>", colnames(ENRTABLE))

    ## Create DT object
    outDT <- datatable(ENRTABLE, rownames = FALSE, extensions = "Buttons", escape = FALSE,
        filter = list(position = "top", clear = FALSE), options = list(columnDefs = list(list(searchable = FALSE,
            orderable = FALSE, width = "3px", targets = c(12, 13, 14)), list(visible = FALSE,
            targets = 11), list(className = "dt-center", targets = "_all")), search = list(regex = TRUE),
            searchCols = list(list(search = dgeHits), list(search = subgroupID),
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL), scrollX = TRUE, scrollY = "325px", dom = "Brtp", paging = TRUE,
            pageLength = 50, buttons = list(list(extend = "collection", text = "Help",
                action = DT::JS("function ( e, dt, node, config ) {
                                           Shiny.setInputValue('hyperHelpMulti', true, {priority: 'event'});
}")),
                list(extend = "collection", text = "Download All Results", action = DT::JS("function ( e, dt, node, config ) {
                                    Shiny.setInputValue('hyperDLMulti', true, {priority: 'event'});
                    }")))),
        selection = "none") %>% formatRound(c("Mean<br>ssGSEA", "Diff<br>ssGSEA"),
        digits = 2) %>% formatSignif(c("P Value<br>Hyper", "FDR<br>Hyper", "P Value<br>ssGSEA",
        "FDR<br>ssGSEA"), digits = 2) %>% formatStyle(c("Subgroup", "N<br>Gene Set",
        "Mean<br>ssGSEA"), `border-right` = "solid 2px")

}

## Generete hyperenrichment results
hyperenrichmentClusters <- function(clusterRes, groupList, genesets, qthresh, cthresh,
    ntotal) {

    ## Create list of gene signatures
    sigList <- lapply(seq(length(groupList)), function(mod, clusterRes, qthresh) {

        ## Get subset of the clusters
        cSub <- clusterRes[clusterRes$edge == mod, ]

        ## Get genes with sig pvalues
        genes <- cSub$gene[cSub$fdr < qthresh & cSub$coef > cthresh]

        return(genes)

    }, clusterRes, qthresh)
    names(sigList) <- names(groupList)

    ## Run hyperenrichment
    gseList <- lapply(sigList, function(sig, genesets, ntotal) {
        enrichFram <- NULL
        if (length(sig) > 0) {
            hits <- vapply(genesets, function(x, y) paste(intersect(x, y), collapse = ","),
                sig, FUN.VALUE = character(1))
            nhits <- vapply(genesets, function(x, y) length(intersect(x, y)), sig,
                FUN.VALUE = integer(1))
            ndrawn <- length(sig)
            ncats <- vapply(genesets, length, FUN.VALUE = integer(1))
            nleft <- ntotal - ncats
            pval <- phyper(q = nhits - 1, m = ncats, n = nleft, k = ndrawn, lower.tail = FALSE)
            enrichFram <- data.frame(category = names(genesets), pval = pval, nhits = nhits,
                ndrawn = ndrawn, ncats = ncats, ntot = ntotal, hits = hits, stringsAsFactors = FALSE)
        }
        return(enrichFram)
    }, genesets, ntotal)

    ## Calculate and merge FDR values
    pValueDF <- data.frame(pval = unlist(lapply(gseList, function(y) y$pval)))
    pValueDF$fdr <- p.adjust(pValueDF$pval, method = "BH")
    pValueDF <- unique(pValueDF)
    gseList <- lapply(gseList, function(y, pValueDF) {
        y <- merge(y, pValueDF)
        if (nrow(y) > 0) {
            y <- y[, c("category", "pval", "fdr", "nhits", "ndrawn", "ncats", "ntot",
                "hits")]
        }
        return(y)
    }, pValueDF)

    return(gseList)
}
