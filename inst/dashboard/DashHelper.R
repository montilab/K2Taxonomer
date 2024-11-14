## Function to generate differential signature
DashSignatureWrapper <- function(eM, cD, cohorts, mods, vehicle=NULL,
                                 variables=NULL, logCounts=FALSE,
                                 DGEexpThreshold = NULL,
                                 GENE = TRUE) {
  
  
  # Use internal functions from R package (or environment when developing)
  if(!exists(".signatureWrapper")) {
    .signatureWrapper <- K2Taxonomer:::.signatureWrapper
  }
  if(!exists(".K2filterGenes")) {
    .K2filterGenes <- K2Taxonomer:::.K2filterGenes
  }
  if(!exists(".formatCov")) {
    .formatCov <- K2Taxonomer:::.formatCov
  }
  if(!exists(".limmaTable")) {
    .limmaTable <- K2Taxonomer:::.limmaTable
  }
  
  # Create K2 object
  if(GENE) {
    K2obj <- new("K2", 
                 eMat=eM,
                 colData=cD)
  } else {
    K2obj <- new("K2", 
                 gMat=eM,
                 colData=cD)
  }
  K2meta(K2obj)$cohorts <- cohorts
  K2meta(K2obj)$vehicle <- vehicle
  K2meta(K2obj)$variables <- cohorts
  K2meta(K2obj)$logCounts <- logCounts
  K2meta(K2obj)$DGEexpThreshold <- DGEexpThreshold
  
  modStats <- .signatureWrapper(K2obj, mods, GENE)$modStats
  modStats$fdr <- p.adjust(modStats$pval, method = "BH")
  
  ## Return
  return(modStats)
}

## Function to format differential results
geneTable <- function(DGETABLE, nodeID=NULL, geneList=NULL) {
  
  ## Get exact match for nodeID
  if (!is.null(nodeID) && nodeID == "") {
    nodeID <- NULL
  }
  if (!is.null(nodeID))
    nodeID <- paste0("^", nodeID, "$")
  
  ## Create data table obect
  datatable(DGETABLE, rownames=FALSE, extensions="Buttons", escape=FALSE,
            filter=list(position="top", clear=FALSE),
            options=list(columnDefs=list(list(searchable=FALSE,
                                              orderable=FALSE, width="3px", targets=c(8, 9, 10)),
                                         list(className="dt-center", targets="_all")),
                         search=list(regex=TRUE), searchCols=list(list(search=geneList),
                                                                  list(search=nodeID), NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                                                  NULL), scrollX=TRUE, scrollY="325px", dom="Brtp", paging=TRUE,
                         pageLength=50, buttons=list(list(extend="collection", text="Help",
                                                          action=DT::JS(
                                                            "function ( e, dt, node, config ) {
                    Shiny.setInputValue('geneHelp', true, {priority: 'event'});
                    }")),
                                                     list(extend="collection", text="Download All Results",
                                                          action=DT::JS(
                                                            "function ( e, dt, node, config ) {
                    Shiny.setInputValue('geneDL', true, {priority: 'event'});
                    }")))),
            selection="none") %>%
    formatRound(c("Mean", "Diff"), digits=2) %>%
    formatSignif(c("P Value", "FDR"), digits=2) %>%
    formatStyle(c("Direction", "Mean"), `border-right`="solid 2px")
}

## Function to format enrichment results
genesetTable <- function(ENRTABLE, nodeID=NULL, edgeID=NULL, dgeHits=NULL) {
  
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
  outDT <- datatable(ENRTABLE, rownames=FALSE, extensions="Buttons",
                     escape=FALSE, filter=list(position="top", clear=FALSE),
                     options=list(columnDefs=list(list(searchable=FALSE,
                                                       orderable=FALSE, width="3px", targets=c(14, 15, 16)),
                                                  list(visible=FALSE, targets=c(12, 13)),
                                                  list(className="dt-center", targets="_all")),
                                  search=list(regex=TRUE), searchCols=list(list(search=dgeHits),
                                                                           list(search=nodeID),
                                                                           list(search=edgeID), NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                                                           NULL, NULL, NULL, NULL, NULL, NULL, NULL), scrollX=TRUE,
                                  scrollY="325px", dom="Brtp", paging=TRUE, pageLength=50,
                                  buttons=list(list(extend="collection",
                                                    text="Help",
                                                    action=DT::JS(
                                                      "function ( e, dt, node, config ) {
                Shiny.setInputValue('enrHelp', true, {priority: 'event'});
                }")),
                                               list(extend="collection", text="Download All Results",
                                                    action=DT::JS(
                                                      "function ( e, dt, node, config ) {
                Shiny.setInputValue('enrDL', true, {priority: 'event'});
                }")))),
                     selection="none") %>%
    formatRound(c("Mean<br>Score", "Diff<br>Score"), digits=2) %>%
    formatSignif(
      c("P Value<br>Fisher", "FDR<br>Fisher", "P Value<br>Score",
        "FDR<br>Score"), digits=2) %>%
    formatStyle(c("Direction", "N<br>Gene Set", "Diff<br>Score"),
                `border-right`="solid 2px")
  
}

plotGenePathwayDash <- function(eM, cD, gene, obs1, obs2, cohorts, vehicle,
                                yaxis="Expression") {
  
  if (gene %in% rownames(eM)) {
    
    ## Format subgroup names
    if (is.null(cohorts)) {
      nams <- colnames(eM)
    } else {
      nams <- as.character(cD[, cohorts])
    }
    nams[nams == vehicle] <- "Vehicle"
    
    ## Create data.frame of expression values
    e <- eM[gene, ]
    df <- data.frame(e=e, ch=nams, stringsAsFactors=FALSE)
    
    ## Subset for obs in groups
    df <- df[df$ch %in% c(obs1, obs2, "Vehicle"), ]
    df$edge <- "Edge:1"
    df$edge[df$ch %in% obs2] <- "Edge:2"
    df$edge[df$ch == "Vehicle"] <- "Vehicle"
    df$edge <- factor(df$edge, levels =c("Combined","Edge:1", "Edge:2", "Vehicle"))
    
    ## Get per Observation means for each cohort
    if(!is.null(cohorts)) {
      dfMeans <- do.call(rbind, lapply(as.character(unique(df$ch)), function(coh) {
        data.frame(ch = coh, me = mean(df[df$ch == coh, "e"]))
      }))
    } else {
      dfMeans <- df
      colnames(dfMeans)[2] <- "me"
    }
    dfMeans <- dfMeans[order(dfMeans$me, decreasing=TRUE), ]
    
    ## Sort levels by mean expression
    df$ch <- factor(df$ch, levels=dfMeans$ch)
    df <- merge(df, dfMeans)
    df[duplicated(df$ch), "me"] <- NA
    
    ## Add levels for boxplots
    df$edge2 <- df$edge
    
    ## Add rows for boxplots
    df2 <- df
    df2$ch <- df$edge
    df2$edge2 <- "Combined"
    df2$e2 <- df2$e
    df2$e <- NA
    
    df2 <- df2[, colnames(df2) != "me"]
    
    ## Get per Observation means for each edge
    dfMeans2 <- do.call(rbind, lapply(sort(as.character(unique(df$edge))), function(ed) {
      data.frame(edge = ed, me = mean(df[df$edge == ed, "e"]))
    }))
    df2 <- merge(df2, dfMeans2)
    df2[duplicated(df2$edge), "me"] <- NA
    
    ## Concatenate
    df$e2 <- NA
    df <- df[df$ch != "Vehicle", ]
    df <- rbind(df, df2[, colnames(df)])
    
    ## Fix levels
    df$ch <- factor(df$ch, levels=c(dfMeans$ch[dfMeans$ch != "Vehicle"],
                                    "Edge:1", "Vehicle", "Edge:2"))
    df$edge <- factor(df$edge, levels=c("Edge:1", "Vehicle", "Edge:2"))
    df$edge2 <- factor(df$edge2, levels=c("Edge:1", "Combined", "Edge:2"))
    
    ## Add column names
    colnames(df) <- c("Observation", "Expression", "Subgroup", "Mean",
                      "Subgroup2", "Expression2")
    
    ## Plot
    if(!is.null(cohorts)) {
      
      df <- df[order(df$Subgroup2, df$Observation),]
      
      df$Expression <- unlist(lapply(as.character(unique(df$Observation)), function(coh) {
        evec <- df[df$Observation == coh, "Expression"]
        if(!coh %in% c("Edge:1", "Vehicle", "Edge:2")) {
          if(length(evec) > 501) {
            quants <- quantile(evec, seq(0, 1, by = 0.002))
            q <- rep(NA, length(evec))
            q[seq(length(quants))] <- quants
          } else {
            q <- evec
          }
        } else {
          q <- rep(NA, length(evec))
        }
        return(q)
      }))
      
      df$Expression2 <- unlist(lapply(as.character(unique(df$Subgroup)), function(s1) {
        evec <- df[df$Subgroup == s1, "Expression2"]
        s2 <- df[df$Subgroup == s1, "Subgroup2"]
        sl <- s2 == "Combined"
        evecSub <- evec[sl]
        if(sum(sl) > 0) {
          if(length(evecSub) > 501) {
            quants <- quantile(evecSub, seq(0, 1, by = 0.002))
            q <- rep(NA, length(evecSub))
            q[seq(length(quants))] <- quants
            if(!sl[1]) {
              q <- c(rep(NA, sum(!sl)), q)
            } else {
              q <- c(q, rep(NA, sum(!sl)))
            }
          } else {
            q <- evec
          }
        } else {
          q <- rep(NA, length(evec))
        }
        return(q)
      }))
      
      df <- df[!(is.na(df$Expression) & is.na(df$Expression2) & is.na(df$Mean)),]
      
      p <- ggplot(data=df, aes(x=Observation, y=Expression)) +
        geom_boxplot(aes(y=Expression2), fill = "white", outliers = FALSE) +
        geom_violin(aes(y=Expression2, fill=Subgroup), alpha = 0.5) +
        geom_boxplot(fill = "white", linewidth = 0.25, outliers = FALSE) +
        geom_violin(aes(fill = Subgroup), adjust = 10, linewidth = 0.25, alpha = 0.5) +
        geom_point(aes(y=Mean, color=Subgroup), shape=3, size=3) +
        facet_grid(~Subgroup2, scales="free_x") +
        scale_colour_manual(values=c(`Edge:1`="darkorange3", Vehicle="grey",
                                     `Edge:2`="darkorchid3")) +
        scale_fill_manual(values=c(`Edge:1`="darkorange", Vehicle="grey",
                                   `Edge:2`="darkorchid1")) +
        scale_x_discrete() +
        scale_y_continuous(name=yaxis) +
        theme_bw() +
        ggtitle(gene) +
        theme(plot.margin=margin(0, 10, 0, 10), legend.position="none",
              axis.text.x=element_text(angle=45, hjust=0, size=15),
              axis.text.y=element_text(size=15), axis.title.x=element_blank(),
              strip.background = element_rect(fill = "white"))
      
      p <- suppressWarnings(
        style(ggplotly(p), hoverinfo = 'none', traces = c(2, 3, 6, 7))
      )
      
      # Remove and edit hover text
      p$x$data <- lapply(p$x$data, function(dat) {
        if ("hoverinfo" %in% names(dat)) {
          
          if(dat$hoverinfo == "none") {
            dat$text <- NULL
          } else {
            dat$text <- gsub("Observation: Edge:1<br />|Observation: Edge:2<br />", "", dat$text)
            dat$text <- gsub("Observation: ", "", dat$text)
          }
          
        }
        return(dat)
      })
      
    } else {
      p <- ggplot(data=df, aes(x=Observation, y=Expression)) +
        geom_boxplot(aes(y=Expression2), fill = "white", outlier.color = "grey90") +
        geom_violin(aes(y=Expression2, fill=Subgroup), alpha = 0.5) +
        geom_point(data = subset(df, Subgroup2 != "Combined"), aes(color=Subgroup), shape=3, size=3) +
        geom_point(data = subset(df, Subgroup2 == "Combined"), aes(y=Mean, color=Subgroup), shape=3, size=3) +
        facet_grid(~Subgroup2, scales="free_x", space = "free_x") +
        scale_colour_manual(values=c(`Edge:1`="darkorange3", Vehicle="grey",
                                     `Edge:2`="darkorchid3")) +
        scale_fill_manual(values=c(`Edge:1`="darkorange", Vehicle="grey",
                                   `Edge:2`="darkorchid1")) +
        scale_x_discrete() +
        scale_y_continuous(name=yaxis) +
        theme_bw() +
        ggtitle(gene) +
        theme(plot.margin=margin(0, 10, 0, 10), legend.position="none",
              axis.text.x=element_text(angle=45, hjust=0, size=15),
              axis.text.y=element_text(size=15), axis.title.x=element_blank(),
              strip.background = element_rect(fill = "white"))
      
      if(length(obs1) > 10 | length(obs2) > 10) {
        p <- p + theme(
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()
        )
      }
      
      
      p <- suppressWarnings(
        style(ggplotly(p), hoverinfo = 'none', traces = c(2, 3))
      )
      
      # Remove and edit hover text
      p$x$data <- lapply(p$x$data, function(dat) {
        if ("hoverinfo" %in% names(dat)) {
          
          if(dat$hoverinfo == "none") {
            dat$text <- NULL
          } else {
            dat$text <- gsub("Observation: Edge:1<br />|Observation: Edge:2<br />", "", dat$text)
            dat$text <- gsub("Observation: ", "", dat$text)
          }
          
        }
        return(dat)
      })
      
    }
    
    rm(df); invisible(gc(verbose = FALSE))
    
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
plotGenePathwayClustersDash <- function(eM, cD, gene, groupList, cohorts, vehicle,
                                        yaxis="Expression") {
  
  if (gene %in% rownames(eM)) {
    
    ## Format group names
    if (is.null(cohorts)) {
      nams <- colnames(eM)
    } else {
      nams <- as.character(cD[, cohorts])
    }
    nams[nams == vehicle] <- "Vehicle"
    
    ## Create data.frame of expression values
    e <- eM[gene, ]
    df <- data.frame(e=e, ch=nams, stringsAsFactors=FALSE)
    
    ## Get clusters
    obs <- unlist(groupList)
    
    ## Subset for obs in groups
    df <- df[df$ch %in% c(obs, "Vehicle"), ]
    df$edge <- "Vehicle"
    for (i in names(groupList)) {
      df$edge[df$ch %in% groupList[[i]]] <- i
    }
    
    ## Get per Observation means for each cohort
    if(!is.null(cohorts)) {
      dfMeans <- do.call(rbind, lapply(as.character(unique(df$ch)), function(coh) {
        data.frame(ch = coh, me = mean(df[df$ch == coh, "e"]))
      }))
    } else {
      dfMeans <- df[, c("ch", "e")]
      colnames(dfMeans)[2] <- "me"
    }
    dfMeans <- dfMeans[order(dfMeans$me, decreasing=TRUE), ]
    
    ## Sort levels by mean expression
    df$ch <- factor(df$ch, levels=dfMeans$ch)
    df <- merge(df, dfMeans)
    df[duplicated(df$ch), "me"] <- NA
    
    ## Add levels for boxplots
    df$edge2 <- df$edge
    
    ## Add rows for boxplots
    df2 <- df
    df2$ch <- df$edge
    df2$edge2 <- "Combined"
    df2$e2 <- df2$e
    df2$e <- NA
    
    df2 <- df2[, colnames(df2) != "me"]
    
    ## Get per Observation means for each edge
    dfMeans2 <- do.call(rbind, lapply(sort(as.character(unique(df$edge))), function(ed) {
      data.frame(edge = ed, me = mean(df[df$edge == ed, "e"]))
    }))
    df2 <- merge(df2, dfMeans2)
    df2[duplicated(df2$edge), "me"] <- NA
    
    ## Concatenate
    df$e2 <- NA
    df <- df[df$ch != "Vehicle", ]
    df <- rbind(df, df2[, colnames(df)])
    
    ## Fix levels
    df$ch <- factor(df$ch, levels=c(dfMeans$ch[dfMeans$ch != "Vehicle"],
                                    "Vehicle", names(groupList)))
    df$edge <- factor(df$edge, levels=c("Vehicle", names(groupList)))
    df$edge2 <- factor(df$edge2, levels=c(names(groupList), "Combined"))
    
    ## Add column names
    colnames(df) <- c("Observation", "Expression", "Group", "Mean",
                      "Group2", "Expression2")
    
    ## Create color manual
    qual_col_pals=brewer.pal.info[brewer.pal.info$category == "qual", ]
    col_vector=unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                             rownames(qual_col_pals)))
    colMan <- c("grey", col_vector[seq(length(groupList))])
    names(colMan) <- c("Vehicle", names(groupList))
    
    ## Create darker colors
    .darkcol = function(col, x){
      colorRampPalette(c(col, "black"))(100)[x]
    }
    colMan2 <- sapply(colMan, .darkcol, 50)
    
    ## Plot
    if(!is.null(cohorts)) {
      
      df <- df[order(df$Group, df$Observation),]
      
      df$Expression <- unlist(lapply(as.character(unique(df$Observation)), function(coh) {
        evec <- df[df$Observation == coh, "Expression"]
        if(coh %in% nams) {
          if(length(evec) > 501) {
            quants <- quantile(evec, seq(0, 1, by = 0.002))
            q <- rep(NA, length(evec))
            q[seq(length(quants))] <- quants
          } else {
            q <- evec
          }
        } else {
          q <- rep(NA, length(evec))
        }
        return(q)
      }))
      
      df$Expression2 <- unlist(lapply(as.character(unique(df$Group)), function(s1) {
        evec <- df[df$Group == s1, "Expression2"]
        s2 <- df[df$Group == s1, "Group2"]
        sl <- s2 == "Combined"
        evecSub <- evec[sl]
        if(sum(sl) > 0) {
          if(length(evecSub) > 501) {
            quants <- quantile(evecSub, seq(0, 1, by = 0.002))
            q <- rep(NA, length(evecSub))
            q[seq(length(quants))] <- quants
            if(!sl[1]) {
              q <- c(rep(NA, sum(!sl)), q)
            } else {
              q <- c(q, rep(NA, sum(!sl)))
            }
          } else {
            q <- evec
          }
        } else {
          q <- rep(NA, length(evec))
        }
        return(q)
      }))
      
      df <- df[!(is.na(df$Expression) & is.na(df$Expression2) & is.na(df$Mean)),]
      
      p <- ggplot(data=df, aes(x=Observation, y=Expression)) +
        geom_boxplot(aes(y=Expression2), fill = "white", outliers = FALSE) +
        geom_violin(aes(y=Expression2, fill=Group), alpha = 0.5) +
        geom_boxplot(fill = "white", linewidth = 0.25, outliers = FALSE) +
        geom_violin(aes(fill = Group), adjust = 10, linewidth = 0.25, alpha = 0.5) +
        geom_point(aes(y=Mean, color=Group), shape=3, size=3) +
        facet_grid(~Group2, scales="free_x") +
        scale_colour_manual(values=colMan2) +
        scale_fill_manual(values=colMan) +
        scale_x_discrete() +
        scale_y_continuous(name=yaxis) +
        theme_bw() +
        ggtitle(gene) +
        theme(plot.margin=margin(0, 10, 0, 10), legend.position="none",
              axis.text.x=element_text(angle=45, hjust=0, size=15),
              axis.text.y=element_text(size=15), axis.title.x=element_blank(),
              strip.background = element_rect(fill = "white"))
      
      p <- suppressWarnings(
        style(ggplotly(p), hoverinfo = 'none', traces = c(2, 3, 6, 7))
      )
      
      # Remove and edit hover text
      p$x$data <- lapply(p$x$data, function(dat) {
        if ("hoverinfo" %in% names(dat)) {
          
          if(dat$hoverinfo == "none") {
            dat$text <- NULL
          } else {
            dat$text <- gsub("Observation: Edge:1<br />|Observation: Edge:2<br />", "", dat$text)
            dat$text <- gsub("Observation: ", "", dat$text)
          }
          
        }
        return(dat)
      })
      
    } else {
      
      p <- ggplot(data=df, aes(x=Observation, y=Expression)) +
        geom_boxplot(aes(y=Expression2), fill = "white", outlier.color = "grey90") +
        geom_violin(aes(y=Expression2, fill=Group), alpha = 0.5) +
        geom_point(data = subset(df, Group2 != "Combined"), aes(color=Group), shape=3, size=3) +
        geom_point(data = subset(df, Group2 == "Combined"), aes(y=Mean, color=Group), shape=3, size=3) +
        facet_grid(~Group2, scales="free_x", space = "free_x") +
        scale_colour_manual(values=colMan2) +
        scale_fill_manual(values=colMan) +
        scale_x_discrete() +
        scale_y_continuous(name=yaxis) +
        theme_bw() +
        ggtitle(gene) +
        theme(plot.margin=margin(0, 10, 0, 10), legend.position="none",
              axis.text.x=element_text(angle=45, hjust=0, size=15),
              axis.text.y=element_text(size=15), axis.title.x=element_blank(),
              strip.background = element_rect(fill = "white"))
      
      if(length(obs1) > 10 | length(obs2) > 10) {
        p <- p + theme(
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )
      }
      
      p <- suppressWarnings(
        style(ggplotly(p), hoverinfo = 'none', traces = c(2, 3))
      )
      
    }
    
    # Remove and edit hover text
    p$x$data <- lapply(p$x$data, function(dat) {
      if ("hoverinfo" %in% names(dat)) {
        
        if(dat$hoverinfo == "none") {
          dat$text <- NULL
        } else {
          dat$text <- gsub("Observation: Edge:1<br />|Observation: Edge:2<br />", "", dat$text)
          dat$text <- gsub("Observation: ", "", dat$text)
        }
        
      }
      return(dat)
    })
    
  }
  
  rm(df); invisible(gc(verbose = FALSE))
  
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

## Function to format differential results
geneTableClusters <- function(clusterRes, subgroupID=NULL, geneList=NULL) {
  
  ## Get exact match for nodeID
  if (!is.null(subgroupID) && subgroupID == "") {
    subgroupID <- NULL
  }
  if (!is.null(subgroupID)) {
    nodeID <- paste0("^", subgroupID, "$")
  }
  
  ## Create data table obect
  datatable(clusterRes, rownames=FALSE, extensions="Buttons", escape=FALSE,
            filter=list(position="top", clear=FALSE),
            options=list(columnDefs=list(list(searchable=FALSE,
                                              orderable=FALSE, width="3px", targets=c(6, 7, 8)),
                                         list(className="dt-center", targets="_all")),
                         search=list(regex=TRUE),
                         searchCols=list(list(search=geneList), list(search=subgroupID),
                                         NULL, NULL, NULL, NULL, NULL), scrollX=TRUE,
                         scrollY="325px", dom="Brtp", paging=TRUE, pageLength=50,
                         buttons=list(list(extend="collection",
                                           text="Help",
                                           action=DT::JS(
                                             "function ( e, dt, node, config ) {
                Shiny.setInputValue('geneHelp', true, {priority: 'event'});
                }")),
                                      list(extend="collection", text="Download All Results",
                                           action=DT::JS(
                                             "function ( e, dt, node, config ) {
                Shiny.setInputValue('geneDLMulti', true, {priority: 'event'});
                }")))),
            selection="none") %>%
    formatRound(c("Mean", "Diff"), digits=2) %>%
    formatSignif(c("P Value", "FDR"), digits=2) %>%
    formatStyle(c("Subgroup", "Mean"), `border-right`="solid 2px")
}

## Function to format enrichment results from multiple comparisons
genesetTableClusters <- function(ENRTABLE, subgroupID=NULL, dgeHits=NULL) {
  
  ## Get exact match for nodeID
  if (!is.null(subgroupID) && subgroupID == "") {
    subgroupID <- NULL
  }
  if (!is.null(subgroupID))
    nodeID <- paste0("^", subgroupID, "$")
  
  ## Add line breaks
  colnames(ENRTABLE) <- gsub("_", "<br>", colnames(ENRTABLE))
  
  ## Create DT object
  outDT <- datatable(ENRTABLE, rownames=FALSE, extensions="Buttons",
                     escape=FALSE, filter=list(position="top", clear=FALSE),
                     options=list(columnDefs=list(list(searchable=FALSE, orderable=FALSE,
                                                       width="3px", targets=c(12, 13, 14)), list(visible=FALSE,
                                                                                                 targets=11), list(className="dt-center", targets="_all")),
                                  search=list(regex=TRUE), searchCols=list(list(search=dgeHits),
                                                                           list(search=subgroupID), NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                                                           NULL, NULL, NULL, NULL, NULL, NULL),
                                  scrollX=TRUE, scrollY="325px", dom="Brtp", paging=TRUE,
                                  pageLength=50, buttons=list(list(extend="collection", text="Help",
                                                                   action=DT::JS(
                                                                     "function ( e, dt, node, config ) {
            Shiny.setInputValue('enrHelpMulti', true, {priority: 'event'});
            }")),
                                                              list(extend="collection", text="Download All Results",
                                                                   action=DT::JS("function ( e, dt, node, config ) {
                Shiny.setInputValue('enrDLMulti', true, {priority: 'event'});
                }")))),
                     selection="none") %>%
    formatRound(c("Mean<br>Score", "Diff<br>Score"), digits=2) %>%
    formatSignif(c("P Value<br>Fisher", "FDR<br>Fisher", "P Value<br>Score",
                   "FDR<br>Score"), digits=2) %>%
    formatStyle(c("Subgroup", "N<br>Gene Set", "Mean<br>Score"),
                `border-right`="solid 2px")
  
}

## Generete enrichment results
enrichmentClusters <- function(clusterRes, groupList, genesets, qthresh,
                               cthresh, ntotal) {
  
  ## Create list of gene signatures
  sigList <- lapply(seq(length(groupList)),
                    function(mod, clusterRes, qthresh) {
                      
                      ## Get subset of the clusters
                      cSub <- clusterRes[clusterRes$edge == mod, ]
                      
                      ## Get genes with sig pvalues
                      genes <- cSub$gene[cSub$fdr < qthresh & cSub$coef > cthresh]
                      
                      return(genes)
                      
                    }, clusterRes, qthresh)
  names(sigList) <- names(groupList)
  
  ## Run enrichment
  gseList <- lapply(sigList, function(sig, genesets, ntotal) {
    enrichFram <- NULL
    if (length(sig) > 0) {
      hits <- vapply(genesets, function(x,
                                        y) paste(intersect(x, y), collapse=","),
                     sig, FUN.VALUE=character(1))
      nhits <- vapply(genesets, function(x,
                                         y) length(intersect(x, y)), sig, FUN.VALUE=integer(1))
      ndrawn <- length(sig)
      ncats <- vapply(genesets, length, FUN.VALUE=integer(1))
      nleft <- ntotal - ncats
      pval <- mapply(function(nh, ns, ng, nt) {
        tr <- ng - nh
        bl <- ns - nh
        tl <- nt - tr - bl - nh
        fisher.test(matrix(c(tl, bl, tr, nh), ncol = 2), 
                    alternative = "greater")$p.value
      }, nhits,  ncats, ndrawn, ntotal)
      enrichFram <- data.frame(category=names(genesets),
                               pval=pval, nhits=nhits, ndrawn=ndrawn,
                               ncats=ncats, ntot=ntotal,
                               hits=hits, stringsAsFactors=FALSE)
    }
    return(enrichFram)
  }, genesets, ntotal)
  
  ## Calculate and merge FDR values
  pValueDF <- data.frame(pval=unlist(lapply(gseList, function(y) y$pval)))
  pValueDF$fdr <- p.adjust(pValueDF$pval, method="BH")
  pValueDF <- unique(pValueDF)
  gseList <- lapply(gseList, function(y, pValueDF) {
    y <- merge(y, pValueDF)
    if (nrow(y) > 0) {
      y <- y[, c("category", "pval", "fdr", "nhits", "ndrawn", "ncats",
                 "ntot", "hits")]
    }
    return(y)
  }, pValueDF)
  
  return(gseList)
}
