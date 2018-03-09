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

#' @export
PlotGeneFraction <- function(gene, r, mtx, title.x=0.5, title.y=0.5, alpha=0.5, size=0.7, show.legend=T, plot.na=F,
                             show.ticks=F, legend.position=c(1, 0), legend.only=F, limits=NULL, class.label.layer=NULL, ...) {
  gg <- PlotPagodaEmbeding(r, colors=mtx[, gene], alpha=alpha, size=size, show.legend=show.legend,
                           plot.na=plot.na, show.ticks=show.ticks, ...) +
    ggrastr::theme_pdf(legend.pos=legend.position, show.ticks=show.ticks) +
    ggplot2::theme(plot.margin=ggplot2::margin(), axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank()) +
    ggplot2::guides(color=ggplot2::guide_colorbar(title='Expression', direction='horizontal',
                                                  title.position='top')) +
    ggplot2::scale_color_distiller(palette="Spectral")

  if (!is.null(limits)) {
    gg <- gg + ggplot2::xlim(limits[1,]) + ggplot2::ylim(limits[2,])
  }

  if (!is.null(class.label.layer)) {
    gg$layers <- c(gg$layers, list(class.label.layer))
    gg <- gg + ggplot2::scale_size_continuous(range=c(3, 3))
  }

  if (legend.only)
    return(cowplot::plot_grid(cowplot::get_legend(gg)))

  return(cowplot::ggdraw(gg) + cowplot::draw_label(gene, x=title.x, y=title.y, hjust=title.x, vjust=title.y))
}

#' @export
BordersOfClusterUnion <- function(cluster.names, clusters, type.ids, embedding) {
  type.ids <- unlist(type.ids[cluster.names])
  borders <- embedding[names(clusters)[clusters %in% type.ids],] %>%
    apply(2, quantile, c(0.01, 0.99)) %>% t()

  borders <- (borders - mean(borders)) * 0.1 + borders

  return(borders)
}
