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
PlotTrimmedCorrections <- function(trimmed.data, raw.data, trimed.length, log=T, adjusted.raw=FALSE) {
  plot.df <- GetPlotExrDf(trimmed.data, raw.data, trimed.length, adjusted.raw=adjusted.raw)
  plot.df.all <- GetPlotExrDf(trimmed.data, raw.data, trimed.length, return.all=T, adjusted.raw=adjusted.raw)

  plot.labs <- ggplot2::labs(x='Corrected #UMIs without trimming', y='Error on trimmed data, %',
                             title=paste0('Length of trimmed UMI: ', trimed.length))

  plot.df.subset <- plot.df.all %>% dplyr::filter(RealValue < 50)
  gg.small <- ggplot2::ggplot(plot.df.subset,
                              ggplot2::aes(x=as.factor(ceiling(RealValue / 5) * 5) ,
                                           y=100 * (value - RealValue) / RealValue,
                                           color=Correction)) +
    ggrastr::geom_boxplot_jitter_outlier(outlier.jitter.width=NULL, outlier.jitter.height=1, outlier.size=0.3, outlier.alpha=0.7) +
    ggplot2::scale_linetype_manual(values=c('dashed', 'solid')) + ggsci::scale_color_jco() +
    ggplot2::ylim(-50, 50) +
    ggplot2::theme(legend.position="none") +
    plot.labs

  trans <- ifelse(log, 'log10', 'identity')
  gg.large <- ggplot2::ggplot(plot.df, ggplot2::aes(x=RealValue, y=100 * (TMean - RealValue) / RealValue, color=Correction)) +
    ggplot2::geom_point(size=0.3, alpha=0.3) +
    ggplot2::geom_smooth(alpha=0.1) +
    ggplot2::scale_linetype_manual(values=c('dashed', 'solid')) +
    ggsci::scale_color_jco(labels=c('Bayesian', 'cluster', 'cluster,\nno equals', 'directional', 'no correction')) +
    ggplot2::scale_x_continuous(expand=c(0, 0), limits=c(1, 7999), trans=trans) +
    ggplot2::scale_y_continuous(expand=c(0, 0), limits=c(-101, 70)) +
    ggplot2::theme(legend.position=c(0, 0), legend.justification=c(0, 0), strip.background=ggplot2::element_blank(),
                   strip.text=ggplot2::element_text(size=14)) +
    ggplot2::guides(color=ggplot2::guide_legend(ncol=3)) +
    plot.labs


  return(list(small=gg.small, large=gg.large))
}
