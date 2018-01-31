#' @export
GetPagoda <- function(cm, n.cores=10, clustering.type='infomap', embeding.type='tSNE', verbose=TRUE) {
  r <- pagoda2::Pagoda2$new(cm, modelType='plain', trim=5, n.cores=n.cores, verbose=verbose)
  r$adjustVariance(plot=F, do.par=F, gam.k=10, verbose=verbose)

  r$calculatePcaReduction(nPcs=100, n.odgenes=1000, maxit=1000)
  r$makeKnnGraph(k=30,type='PCA', center=T,distance='cosine',weight.type='none', verbose=verbose)
  if (clustering.type == 'infomap') {
    r$getKnnClusters(method=igraph::infomap.community,type='PCA',name='infomap')
  } else if (clustering.type == 'multilevel') {
    r$getKnnClusters(method=igraph::multilevel.community,type='PCA',name='multilevel')
  } else stop("Unknown clustering  type")

  if ('largeVis' %in% embeding.type) {
    r$getEmbedding(type='PCA', embeddingType = 'largeVis')
  }

  if ('tSNE' %in% embeding.type) {
    r$getEmbedding(type='PCA', perplexity=30, embeddingType = 'tSNE')
  }

  return(r)
}

#' @export
PlotPagodaEmbeding <- function(r, embeding.type='tSNE', clusters=NULL, clustering.type=NULL, colors=NULL, plot.na=TRUE,
                               min.cluster.size=0, mark.clusters=TRUE, show.legend=FALSE, alpha=0.4, size=0.8, title=NULL,
                               font.size=5.5, show.ticks=TRUE, raster=FALSE, raster.width=NULL, raster.height=NULL, raster.dpi=300,
                               return.df=FALSE, ...) {
  labels <- ggplot2::labs(x='Component 1', y='Component 2')
  plot.df <- tibble::rownames_to_column(as.data.frame(r$embeddings$PCA[[embeding.type]]), var='CellName')
  if (raster) {
    geomp_point_w <- function(...) ggrastr::geom_point_rast(..., width=raster.width, height=raster.height, dpi=raster.dpi)
  } else {
    geomp_point_w <- ggplot2::geom_point
  }

  if (is.null(clusters) & !is.null(clustering.type)) {
    clusters <- r$clusters$PCA[[clustering.type]]
  }

  if (is.null(colors) & !is.null(clusters)) {
    plot.df <- plot.df %>% dplyr::mutate(Cluster=clusters[CellName])

    plot.df$Cluster <- as.character(plot.df$Cluster)

    big.clusts <- (plot.df %>% dplyr::group_by(Cluster) %>% dplyr::summarise(Size=n()) %>%
      dplyr::filter(Size >= min.cluster.size))$Cluster %>% as.vector()

    plot.df$Cluster[!(plot.df$Cluster %in% big.clusts)] <- NA
    na.plot.df <- plot.df %>% filter(is.na(Cluster))
    plot.df <- plot.df %>% filter(!is.na(Cluster))

    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=V1, y=V2)) +
      geomp_point_w(ggplot2::aes(col=Cluster), alpha=alpha, size=size) +
      labels

    if (plot.na) {
      gg <- gg + geomp_point_w(data=na.plot.df, alpha=alpha, size=size, color='black', shape=4)
    }

    if (mark.clusters) {
      labels.data <- plot.df %>% dplyr::group_by(Cluster) %>%
        dplyr::summarise(V1=mean(V1, tirm=0.4), V2=mean(V2, trim=0.4), Size=n())

      if (is.null(font.size)) {
        gg_repel <- ggrepel::geom_label_repel(data=labels.data, ggplot2::aes(label=Cluster, size=Size), color='black',
                                              fill=ggplot2::alpha('white', 0.7), label.size = NA,
                                              label.padding=ggplot2::unit(1, "pt"), ...)
      } else {
        gg_repel <- ggrepel::geom_label_repel(data=labels.data, ggplot2::aes(label=Cluster), color='black', size=font.size,
                                              fill=ggplot2::alpha('white', 0.7), label.size = NA,
                                              label.padding=ggplot2::unit(1, "pt"), ...)
      }
      gg <- gg + gg_repel +
        ggplot2::scale_size_continuous(range=c(3, 7), trans='identity', guide='none')
    }
  } else if (!is.null(colors)) {
    plot.df <- plot.df %>% dplyr::mutate(Color=colors[CellName])
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=V1, y=V2)) +
      geomp_point_w(ggplot2::aes(col=Color), alpha=alpha, size=size) +
      labels
  } else {
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=V1, y=V2)) +
      geomp_point_w(alpha=alpha, size=size) +
      labels
  }

  if (return.df)
    return(plot.df)

  if (!is.null(title)) {
    gg <- gg + ggplot2::ggtitle(title)
  }

  if (!show.legend) {
    gg <- gg + ggplot2::theme(legend.position="none")
  }

  if (!show.ticks) {
    gg <- gg + ggplot2::theme(axis.ticks=ggplot2::element_blank(), axis.text=ggplot2::element_blank())
  }

  return(gg)
}
