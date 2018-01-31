#' @export
BuildPanel4 <- function(gg.plots, ylabel, xlabel, show.legend=F, return.raw=F, show.ticks=T,
                        labels=c('A', 'B', 'C', 'D'), plot.theme=NULL, ...) {
  margin.theme <- ggplot2::theme(plot.margin=ggplot2::margin(l=0.03, r=0.03, b=0.03, t=0.06, "in"))

  gg.plots <- lapply(gg.plots, function(gg) gg + ggrastr::theme_pdf(show.ticks=show.ticks) +
                       margin.theme +  ggpubr::rremove('xylab') + ggpubr::rremove('legend'))

  gg.plots[[1]] <- gg.plots[[1]] + ggpubr::rremove("x.ticks") + ggpubr::rremove("x.text")
  gg.plots[[2]] <- gg.plots[[2]] + ggpubr::rremove("ticks") + ggpubr::rremove("xy.text")
  if (show.legend) {
    gg.plots[[3]] <- gg.plots[[3]] + ggplot2::theme(legend.position=c(0, 0), legend.justification=c(0, 0))
  }
  gg.plots[[4]] <- gg.plots[[4]] + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text")

  if (!is.null(plot.theme)) {
    gg.plots <- lapply(gg.plots, `+`, plot.theme)
  }

  if (return.raw)
    return(gg.plots)

  gg.res <- ggpubr::annotate_figure(cowplot::plot_grid(plotlist=gg.plots, ncol=2, nrow=2, labels=labels, ...),
                            left=ggpubr::text_grob(ylabel, size=14, rot=90),
                            bottom=ggpubr::text_grob(xlabel, size=14))

  return(gg.res)
}

#' @export
HeatmapAnnotGG <- function(df, umi.per.cell.limits=c(2, 4.5), raster.width=5, raster.height=5, annot.width=0.05,
                           raster.dpi=300, palette=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)[1:80]) {
  gg <- ggplot2::ggplot(df, ggplot2::aes(y=Barcode)) +
    ggrastr::theme_pdf() +
    ggpubr::rremove('xy.text') + ggpubr::rremove('ticks')

  heatmap.width <- 1 - 2 * annot.width

  ggs <- list(
    heatmap = gg +
      ggrastr::geom_tile_rast(ggplot2::aes(x=Gene, fill=Expression),
                              width=raster.width * heatmap.width,
                              height=raster.height, dpi=raster.dpi) +
      ggplot2::scale_fill_gradientn(colours=palette, values=),
    clust = gg + ggrastr::geom_tile_rast(ggplot2::aes(x=1, fill=Cluster),
                                         width=raster.width * annot.width,
                                         height=raster.height, dpi=raster.dpi) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_fill_discrete(drop=F) +
      ggplot2::theme(plot.margin=ggplot2::margin()) + ggpubr::rremove('xylab'),
    umis = gg + ggrastr::geom_tile_rast(ggplot2::aes(x=1, fill=UmisPerCb),
                                        width=raster.width * annot.width,
                                        height=raster.height, dpi=raster.dpi) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_fill_distiller(palette='OrRd', limits=umi.per.cell.limits, direction=1) +
      ggplot2::theme(plot.margin=ggplot2::margin()) + ggpubr::rremove('xylab')
  )

  return(ggs)
}

#' @export
HeatmapLegendGuide <- function(title, barwidth=1.3, guide=ggplot2::guide_colorbar, ...) {
  ggplot2::guides(fill = guide(title.position='top', direction='horizontal', title=title,
                               barwidth=ggplot2::unit(barwidth, 'in'), ...))
}

#' @export
PlotFiltrationResults <- function(pgd, clusters, rescued.clusters, filtered.cbs=NULL, raster.width=NULL, raster.height=NULL,
                                  rescued.alpha=0.9, rescued.size=2, rescued.stroke=0.4, unchanged.alpha=0.4,
                                  unchanged.size=0.5) {
  clusters <- as.factor(clusters)

  if (!is.null(filtered.cbs)) {
    filt.clusters <- rep("none", length(filtered.cbs)) %>% setNames(filtered.cbs)
    filt.df <- PlotPagodaEmbeding(pgd, clusters=filt.clusters, return.df=T, plot.na=F)
  }

  rescued.df <- PlotPagodaEmbeding(pgd, clusters=rescued.clusters, return.df=T)

  gg <- PlotPagodaEmbeding(pgd, clusters=clusters, alpha=unchanged.alpha, size=unchanged.size, font.size=NULL, plot.na=F,
                           raster=T, raster.width=raster.width, raster.height=raster.height,
                           point.padding=ggplot2::unit(0, 'pt'), nudge_y=2)

  gg$layers[[1]]$mapping$shape <- 'unchanged'

  if(is.null(filtered.cbs)) {
    gg.filt <- NULL
  } else {
    gg.filt <- ggrastr::geom_point_rast(data=filt.df, mapping=ggplot2::aes(x=V1, y=V2, shape='filtered'), size=1, alpha=0.5,
                                        width=raster.width, height=raster.height)
  }
  gg$layers <- c(gg$layers[[1]],
                 gg.filt,
                 ggrastr::geom_point_rast(data=rescued.df, mapping=ggplot2::aes(x=V1, y=V2, shape='rescued', fill=Cluster),
                                 size=rescued.size, alpha=rescued.alpha, stroke=rescued.stroke, color='black',
                                 width=raster.width, height=raster.height),
                 gg$layers[[2]])

  color.values <- scales::hue_pal()(length(levels(clusters))) %>% setNames(levels(clusters))
  gg.tsne <- gg + ggplot2::scale_color_discrete(guide="none") +
    ggplot2::scale_shape_manual(values=c(unchanged=4, rescued=24, filtered=19), name='Cell filtration') +
    ggplot2::scale_size_continuous(range=c(3, 5), trans='identity') +
    ggplot2::scale_fill_manual(values=color.values) +
    ggplot2::scale_color_manual(values=color.values) +
    ggplot2::guides(fill="none", color="none", size="none") +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank())

  return(gg.tsne)
}

#' @export
theme_base <- ggplot2::theme_bw(base_size=14, base_family='Helvetica') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
