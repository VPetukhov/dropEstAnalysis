#' @export
PlotClustering <- function(pgd, clusters) {
  PlotPagodaEmbeding(pgd, clusters=clusters, show.ticks=F, size=1, alpha=0.5, nudge_x=1, nudge_y=1, font.size=NULL) +
    ggplot2::scale_size_continuous(range=c(3, 6))
}

#' @export
AnnotateClusters <- function(clusters, type.clusters) {
  clusters.annotated <- as.integer(clusters) %>% setNames(names(clusters))
  for (n in names(type.clusters)) {
    clusters.annotated[clusters.annotated %in% type.clusters[[n]]] <- n
  }

  return(clusters.annotated)
}

#' @export
PlotExpressionHeatmap <- function(matrix, clusters, genes) {
  clusters <- data.frame(Cluster=as.factor(sort(clusters)))
  pheatmap::pheatmap(t(matrix[rownames(clusters), genes]), annotation_col=clusters,
                     cluster_rows=F, cluster_cols=F, show_colnames=F)
}
