#' @export
BuildCountMatrixFromReads <- function(filt.rpus, reads.per.umi.per.cb.info, collisions.info=NULL) {
  filt.umis.per.gene <- sapply(filt.rpus, length)
  if (!is.null(collisions.info)) {
    filt.umis.per.gene <- collisions.info[filt.umis.per.gene]
  }
  reads.per.umi.per.cb.info$umis_per_gene <- filt.umis.per.gene
  return(dropestr::BuildCountMatrix(reads.per.umi.per.cb.info))
}

#' @export
PlotCorrectionSize <- function(cms, correction.colors, xlim=NULL, ylim=NULL, min.umi.per.gene=10,
                               facet=T, mapping=NULL, ...) {
  if (is.null(mapping)) {
    mapping <- ggplot2::aes(x=`no correction`, y=`no correction`-value, color=Correction)
  }
  plot.df <- lapply(cms, function(cm) cm@x) %>% tibble::as_tibble() %>%
    dplyr::filter(`no correction` > min.umi.per.gene) %>%
    reshape2::melt(id.vars='no correction', variable.name='Correction')

  gg <- ggplot2::ggplot(plot.df) +
    ggrastr::geom_point_rast(mapping, size = 0.2, alpha=0.05, ...) +
    ggplot2::geom_abline(aes(intercept=0, slope=1)) +
    ggplot2::scale_x_log10(limits=xlim, expand=c(0, 0)) +
    ggplot2::scale_y_log10(limits=ylim, expand=c(0, 0)) +
    ggplot2::annotation_logticks(sides='b') +
    ggplot2::scale_color_manual(values=correction.colors) +
    ggrastr::theme_pdf(legend.pos = c(1, 0)) +
    ggplot2::theme(strip.text.x=ggplot2::element_blank(), panel.spacing=ggplot2::unit(3, 'pt')) +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(alpha=1.0, size=1.5), title='Correction'))

  if (facet) {
    gg <- gg + ggplot2::facet_wrap(~Correction)
  }

  return(gg)
}

#' @export
SampleNoReps <- function(size, ids, probs) {
  umis <- unique(sample(ids, size=size, prob=probs, replace=T))
  while (length(umis) < size) {
    umis <- unique(c(umis, sample(ids, size=size, prob=probs, replace=T)))
  }

  return(umis[1:size])
}

# UMI trimming
GetExprDf <- function(filt.umis.per.gene.vec, real.umis.per.gene.vec, raw.umis.per.gene.vec) {
  expressions <- as.data.frame(filt.umis.per.gene.vec)
  expressions$RealValue <- real.umis.per.gene.vec
  expressions$NoCorrection <- raw.umis.per.gene.vec
  return(expressions %>% dplyr::filter(NoCorrection > 1))
}

GetPlotExrDf <- function(trimmed.cur, raw, umi.length, return.all=FALSE, adjusted.raw=FALSE) {
  raw.name <- ifelse(adjusted.raw, 'umis.per.gene.adj', 'umis.per.gene')
  expr.df <- GetExprDf(trimmed.cur$filt_cells, raw$filt_umis_per_gene_simple, trimmed.cur[[raw.name]]) %>%
    reshape2::melt(id.vars=c('RealValue'), variable.name='Correction') %>% dplyr::mutate(UmiLen=umi.length)

  if (return.all)
    return(expr.df)

  return(expr.df %>% dplyr::group_by(Correction, RealValue) %>%
           dplyr::summarise(Max=quantile(value, 0.95), Min=quantile(value, 0.05), TMean=mean(value, trim=0,2)))
}

#' @export
PlotTrimmedCorrections <- function(trimmed.data, raw.data, trimed.length, log=T, adjusted.raw=FALSE, raster=T,
                                   rast.width=NULL, rast.height=NULL, rast.dpi=300, heights.ratio=c(1, 1)) {
  if (raster) {
    geom_point_w <- function(...) ggrastr::geom_point_rast(..., width=rast.width, height=rast.height * heights.ratio[1], dpi=rast.dpi)
  } else {
    geom_point_w <- ggplot2::geom_point
  }
  heights.ratio <- heights.ratio / sum(heights.ratio)

  plot.df <- GetPlotExrDf(trimmed.data, raw.data, trimed.length, adjusted.raw=adjusted.raw)
  plot.df.all <- GetPlotExrDf(trimmed.data, raw.data, trimed.length, return.all=T, adjusted.raw=adjusted.raw)

  plot.labs <- ggplot2::labs(x='Corrected #UMIs without trimming', y='Error on trimmed data, %',
                             title=paste0('Length of trimmed UMI: ', trimed.length))

  plot.df.subset <- plot.df.all %>% dplyr::filter(RealValue < 50) %>%
    dplyr::mutate(RealValueRound=ceiling(RealValue / 5) * 5)
  title.theme <- ggplot2::theme(plot.title=ggplot2::element_text(margin=ggplot2::margin(0, 0, 0.03, 0, 'in')))

  labels <- c('Bayesian', 'cluster', 'cluster-neq', 'directional', 'no correction')
  color.scale <- ggplot2::scale_color_manual(values=c("#017A5A", "#9B3BB8", "#E69F00", "#BD5500", '#757575'), labels=labels)
  # plot.df.subset$RealValueRound[plot.df.subset$RealValue == 1] <- 1

  trans <- if(log) 'log10' else 'identity'
  gg.large <- ggplot2::ggplot(plot.df,
                              ggplot2::aes(x=RealValue,
                                           y=100 * (TMean - RealValue) / RealValue,
                                           color=Correction, linetype=Correction)) +
    geom_point_w(size=0.3, alpha=0.3) +
    ggplot2::geom_smooth(alpha=0.1) +
    ggplot2::scale_linetype_manual(values=c(rep('solid', 4), 'dashed'), guide=F) +
    color.scale +
    ggplot2::scale_x_continuous(expand=c(0, 0), limits=c(1, 7999), trans=trans) +
    ggplot2::scale_y_continuous(expand=c(0, 0), limits=c(-101, 50)) +
    ggplot2::theme(legend.position=c(0, 0), legend.justification=c(0, 0),
                   strip.background=ggplot2::element_blank(),
                   strip.text=ggplot2::element_text(size=14)) +
    title.theme +
    ggplot2::guides(color=ggplot2::guide_legend(ncol=3)) +
    plot.labs

  gg.small <- ggplot2::ggplot(plot.df.subset,
                              ggplot2::aes(x=as.factor(RealValueRound), y=100 * (value - RealValue) / RealValue, color=Correction)) +
    ggrastr::geom_boxplot_jitter(outlier.jitter.width=NULL, outlier.jitter.height=1, outlier.size=0.3,
                                 outlier.alpha=0.7, raster=raster, raster.width=rast.width,
                                 raster.height=rast.height * heights.ratio[2], raster.dpi=rast.dpi) +
    color.scale +
    ggplot2::ylim(-50, 50) +
    ggplot2::theme(legend.position="none") + title.theme +
    plot.labs

  return(list(small=gg.small, large=gg.large))
}
