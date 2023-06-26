plotEnrichment <- function(pathway, stats,
                           gseaParam=1,
                           ticksSize=0.2) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_point(color="green", size=0.1) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_line(color="green") + theme_bw() +
    
    # geom_hline(yintercept=0, colour="black") +
    # geom_segment(data=data.frame(x=pathway),
    #              mapping=aes(x=x, y=-diff/2,
    #                          xend=x, yend=diff/2),
    #              size=ticksSize) +
    
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    
    labs(x="rank", y="enrichment score")
  g
}
return_gsea_df <- function(pathway, stats, gseaParam=1){
  
}
print_gsea_plot <- function(pathway, stats, gseaParam=1, ticksSize=0.2, path){
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_point(color="#228c22", size=0.1) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_line(color="#228c22") + theme_bw() +
    
    # geom_hline(yintercept=0, colour="black") +
    # geom_segment(data=data.frame(x=pathway),
    #              mapping=aes(x=x, y=-diff/2,
    #                          xend=x, yend=diff/2),
    #              size=ticksSize) +
    ggtitle(names(pathway)) +
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    
    labs(x="rank", y="enrichment score")
  
  g_line <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_hline(yintercept=0, colour="black") +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,
                             xend=x, yend=diff/2),
                 size=ticksSize) +
    theme_bw() +
    theme(panel.border=element_blank(),
          panel.grid=element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())
  
  print_plot_to_ppt(g, path, width=5, height=4)
  print_plot_to_ppt(g_line, path, width=5, height=1)
}
