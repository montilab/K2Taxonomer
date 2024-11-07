plotGenePathway <- function(K2res,
                            feature, 
                            node, 
                            type = c("eMatDS", "eMat", "gMat"),
                            use_plotly = TRUE,
                            subsample = TRUE) {
  
  type <- match.arg(type)
  
  # Get expression/pathway data
  if(type == "eMatDS") {
    eM <- K2eMatDS(K2res)
    yaxis <- "Expression"
    if(nrow(eM) == 0) {
      cat("K2eMatDS() is empty, using K2eMat() instead\n")
      eM <- K2eMat(K2res)
    }
  }
  
  if(type == "eMat") {
    eM <- K2eMat(K2res)
    yaxis <- "Expression"
  }
  
  if(type == "gMat") {
    eM <- K2gMat(K2res)
    if(nrow(eM) == 0) {
      stop("K2gMat() is empty, first use runScoreGenesets().")
    }
    if(K2meta(K2res)$ScoreGeneSetMethod == "GSVA") {
      yaxis <- "GSVA Score"
    } else {
      yaxis <- "Log AUCell Score"
    }
  }
  
  if(!feature %in% rownames(eM)) {
    stop("Feature not in expression/pathway data.")
  }
    
  ## Format subgroup names
  if (is.null(K2meta(K2res)$cohorts)) {
    nams <- colnames(eM)
  } else {
    nams <- as.character(K2colData(K2res)[, K2meta(K2res)$cohorts])
  }
  nams[nams == K2meta(K2res)$vehicle] <- "Vehicle"
  
  ## Create data.frame of expression values
  e <- eM[feature, ]
  df <- data.frame(e=e, ch=nams, stringsAsFactors=FALSE)
  
  ## Get node edge members
  obsList <- K2results(K2res)[[node]]$obs
  obs1 <- obsList[[1]]
  obs2 <- obsList[[2]]
  
  ## Subset for obs in groups
  df <- df[df$ch %in% c(obs1, obs2, "Vehicle"), ]
  df$edge <- "Edge:1"
  df$edge[df$ch %in% obs2] <- "Edge:2"
  df$edge[df$ch == "Vehicle"] <- "Vehicle"
  df$edge <- factor(df$edge, levels =c("Combined","Edge:1", "Edge:2", "Vehicle"))
  
  ## Get per Observation means for each cohort
  if(!is.null(K2meta(K2res)$cohorts)) {
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
  if(!is.null(K2meta(K2res)$cohorts)) {
    
    if(subsample) {
    
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
    }
    
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
      ggtitle(feature) +
      theme(plot.margin=margin(0, 10, 0, 10), legend.position="none",
            axis.text.x=element_text(angle=45, hjust=0, size=15),
            axis.text.y=element_text(size=15), axis.title.x=element_blank(),
            strip.background = element_rect(fill = "white"))
    
    if(use_plotly) {
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
      
      whXaxis <- which(grepl("xaxis", names(p$x$layout)))
      for (i in whXaxis) {
        ticktext <- p$x$layout[[i]]$ticktext
        if (length(ticktext) == 1) {
          p$x$layout[[i]]$tickvals <- seq(2)
          p$x$layout[[i]]$ticktext <- c(ticktext, "")
        }
      }
    }
    
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
      ggtitle(feature) +
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
    
    if(use_plotly) {
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
      
      ## Fix xaxis due to a bug in plotly
      whXaxis <- which(grepl("xaxis", names(p$x$layout)))
      for (i in whXaxis) {
        ticktext <- p$x$layout[[i]]$ticktext
        if (length(ticktext) == 1) {
          p$x$layout[[i]]$tickvals <- seq(2)
          p$x$layout[[i]]$ticktext <- c(ticktext, "")
        }
      }
    }
  }
  
  rm(df); invisible(gc(verbose = FALSE))
  
  if(use_plotly) {
    return(p)
  } else {
    return(suppressWarnings(print(p)))
  }
  
}