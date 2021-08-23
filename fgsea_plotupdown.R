# Code edited by Ian Dryg August 2021
# Original code is the plotEnrichment() function from fgsea. 
# source: https://github.com/ctlab/fgsea/blob/master/R/plot.R
# I am altering it to take in TWO pathways instead of one.
# One should be the UP pathway, and one the DOWN pathway. 
# The goal is to plot them both on the same enrichment plot. 


#' Function: plotTwoEnrichments
#' Plots GSEA enrichment plot.
#' @param pathway1 Gene set to plot (up)
#' @param pathway2 Gene set to plot (down).
#' @param stats Gene-level statistics.
#' @param gseaParam GSEA parameter.
#' @param ticksSize width of vertical line corresponding to a gene (default: 0.2)
#' @return ggplot object with the enrichment plot.
#' @export
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' \dontrun{
#' plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
#'                exampleRanks)
#' }
plotTwoEnrichments <- function(pathway1, pathway2,
                               stats, 
                           gseaParam=1,
                           ticksSize=0.2) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway1 <- unname(as.vector(na.omit(match(pathway1, names(statsAdj)))))
  pathway1 <- sort(pathway1)
  pathway2 <- unname(as.vector(na.omit(match(pathway2, names(statsAdj)))))
  pathway2 <- sort(pathway2)
  
  print('statsAdj - - - - - - - - -')
  print(statsAdj)
  print('pathways2 and 2 - - - - - - - -')
  print(pathway1)
  print(pathway2)
  
  gseaRes1 <- calcGseaStat(statsAdj, selectedStats = pathway1,
                          returnAllExtremes = TRUE)
  gseaRes2 <- calcGseaStat(statsAdj, selectedStats = pathway2,
                           returnAllExtremes = TRUE)
  
  print('gseaRes1 and 2 - - - - - - -')
  print(gseaRes1)
  print(gseaRes2)
  
  bottoms1 <- gseaRes1$bottoms
  tops1 <- gseaRes1$tops
  bottoms2 <- gseaRes2$bottoms
  tops2 <- gseaRes2$tops
  
  print('bottoms1, tops1, bottoms2, tops2 - - - - - - -')
  print(bottoms1)
  print(tops1)
  print(bottoms2)
  print(tops2)
  
  n <- length(statsAdj)
  xs1 <- as.vector(rbind(pathway1 - 1, pathway1))
  ys1 <- as.vector(rbind(bottoms1, tops1))
  toPlot1 <- data.frame(x=c(0, xs1, n + 1), y=c(0, ys1, 0))
  xs2 <- as.vector(rbind(pathway2 - 1, pathway2))
  ys2 <- as.vector(rbind(bottoms2, tops2))
  toPlot2 <- data.frame(x=c(0, xs2, n + 1), y=c(0, ys2, 0))
  
  diff1 <- (max(tops1) - min(bottoms1)) / 8
  diff2 <- (max(tops2) - min(bottoms2)) / 8
  
  # prepare text grob for custom legend
  legend_grob1 <- textGrob("Exhaustion: Up", x=0.07, y = 0.955,
                           gp=gpar(col="black", fontsize=12),
                           just="left")
  legend_grob2 <- textGrob("Exhaustion: Down", x=0.07, y = 0.925,
                           gp=gpar(col="black", fontsize=12),
                           just="left")
    
  # Getting rid of NOTEs
  x=y=NULL
  # initialize ggplot
  g <- ggplot() +
    
    # plot toPlot1, adding points then lines connecting them. In green
    geom_point(data=toPlot1, aes(x=x, y=y), color="green", size=0.1) +
    geom_line(data=toPlot1, aes(x=x, y=y), color="green", size=2) + theme_bw() +
    
    # plot toPlot2, adding points then lines connecting them. In red
    geom_point(data=toPlot2, aes(x=x, y=y), color="red", size=0.1) +
    geom_line(data=toPlot2, aes(x=x, y=y), color="red", size=2) + theme_bw() +
    
    # add the x-axis... a black line along y=0 
    geom_hline(yintercept=0, colour="black") +
    
    # add the tick-marks denoting the leading edge genes
    # I want to put these below the plot, above the x-label
    geom_segment(data=data.frame(x=pathway1),
                 #mapping=aes(x=x, y=-diff1/2 - 0.9,
                 aes(x=x, y=-diff1/2 - 0.9,
                             xend=x, yend=diff1/2 - 0.9),
                             #color="green",
                 size=ticksSize) +
    
    geom_segment(data=data.frame(x=pathway2),
                 #mapping=aes(x=x, y=-diff2/2 - 0.91 - diff1,
                 aes(x=x, y=-diff2/2 - 0.9 - diff1,
                             xend=x, yend=diff2/2 - 0.9 - diff1),
                             #color="red",
                 size=ticksSize) +
    
    # dots denoting the grouping for the tickmarks
    geom_point(aes(x=-500, y = -0.9), color="green", size=5) +
    geom_point(aes(x=-500, y = -0.9-diff1), color="red", size=5) +
    
    # custom legend on top
    geom_point(aes(x=-500, y = 1.15), color="green", size=5) +
    geom_point(aes(x=-500, y = 1.075), color="red", size=5) +
    annotation_custom(legend_grob1) +
    annotation_custom(legend_grob2) +
    
    # customize border and grid to blank
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    
    # add x and y labels
    labs(x="Rank in gene list", y="Running Enrichment Score")
  g
}
