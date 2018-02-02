#' @export
ScoringFunction <- function(clf.funcs, data, real.cbs, background.cbs) {
  tr.data <- data[c(real.cbs, background.cbs),]
  tr.answers <- c(rep(1, length(real.cbs)), rep(0, length(background.cbs)))
  clf <- clf.funcs$train(tr.data, tr.answers)
  return(clf.funcs$predict(clf, data))
}

#' @export
CvResultsTable <- function(cv.res)  {
  df <- lapply(names(cv.res), function(n)
    cbind(cv.res[[n]], Classifier=n, Measure=rownames(cv.res[[n]]))) %>% dplyr::bind_rows() %>%
    dplyr::mutate(txt=paste0(round(mean * 100, 1), ' (±', round(sd * 100, 1), ')')) %>%
    reshape2::dcast(Classifier ~ Measure, value.var='txt')
  return(df)
}

#' @export
PlotTestedClassifierErrors <- function(scores.with.err, wrong.frac.vals, wrong.labels.vars, var.types, filt.subset = NULL, measure.names=NULL) {
  plot.df <- lapply(scores.with.err, lapply, GetTraialsSummary) %>%
    lapply(lapply, function(df) lapply(wrong.labels.vars, function(vars) cbind(df[vars,], Measure=var.types)) %>%
             dplyr::bind_rows(.id="Subset")) %>%
    lapply(function(df) lapply(1:length(wrong.frac.vals), function(i) cbind(df[[i]], Offset=wrong.frac.vals[i])) %>% dplyr::bind_rows()) %>%
    dplyr::bind_rows(.id="Classifier") %>% dplyr::filter(Subset != "all")

  if (!is.null(filt.subset)) {
    plot.df <- plot.df %>% dplyr::filter(Subset == filt.subset)
  }

  if (!is.null(measure.names)) {
    plot.df$Measure <- as.character(plot.df$Measure)
    for (n in names(measure.names)) {
      plot.df$Measure[plot.df$Measure == n] <- measure.names[n]
    }

    plot.df$Measure <- as.factor(plot.df$Measure)
  }

  res <- ggplot2::ggplot(plot.df, ggplot2::aes(x=Offset, ymin=pmax(mean-sd, 0), ymax=pmin(mean+sd, 1), y=mean,
                                               color=Classifier, shape=Measure)) +
    ggplot2::geom_pointrange(ggplot2::aes(linetype=Measure), fatten=2.5, alpha=0.8, size=0.8) +
    ggplot2::geom_line(ggplot2::aes(linetype=Measure), alpha=0.7, size=0.8) +
    ggplot2::scale_color_manual(values = c("#4E58E0", "#F08000", "#349147")) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::labs(x='Error fraction', y='Measure')

  if (is.null(filt.subset)) {
    res <- res + ggplot2::facet_wrap(~ Subset)
  }

  return(res)
}

#' @export
ScoreCellsWithWrongLabels <- function(data, answers, clf, wrong.good.frac, wrong.bad.frac, test.frac) {
  MSE <- function(arr) mean(arr^2)
  good.cbs <- names(answers)[answers > 0.5]
  bad.cbs <- names(answers)[answers < 0.5]

  test.good.cbs <- sample(good.cbs, round(length(good.cbs) * test.frac))
  test.bad.cbs <- sample(bad.cbs, round(length(bad.cbs) * test.frac))

  good.cbs <- setdiff(good.cbs, test.good.cbs)
  bad.cbs <- setdiff(bad.cbs, test.bad.cbs)

  wrong.bad.cbs <- sample(good.cbs, round(length(good.cbs) * wrong.bad.frac))
  wrong.good.cbs <- sample(bad.cbs, round(length(bad.cbs) * wrong.good.frac))

  train.good.cbs <- c(wrong.good.cbs, setdiff(good.cbs, wrong.bad.cbs))
  train.bad.cbs <- c(wrong.bad.cbs, setdiff(bad.cbs, wrong.good.cbs))

  score.with.err <- ScoringFunction(clf, data, train.good.cbs, train.bad.cbs)

  return(list(wrong.fnr=mean(score.with.err[wrong.bad.cbs] < 0.5),
              wrong.fpr=mean(score.with.err[wrong.good.cbs] > 0.5),
              wrong.mse=MSE(c(score.with.err[wrong.good.cbs], 1 - score.with.err[wrong.bad.cbs])),
              all.fnr=mean(score.with.err[good.cbs] < 0.5),
              all.fpr=mean(score.with.err[bad.cbs] > 0.5),
              all.mse=MSE(c(1 - score.with.err[good.cbs],  score.with.err[bad.cbs])),
              test.fnr=mean(score.with.err[test.good.cbs] < 0.5),
              test.fpr=mean(score.with.err[test.bad.cbs] > 0.5),
              test.mse=MSE(c(1 - score.with.err[test.good.cbs],  score.with.err[test.bad.cbs]))))
}

#' @export
ScoreDependsOnBorders <- function(data, answers, clf, good.cbs, bad.cbs, offset.l=0, offset.r=0, mc.cores=1) {
  answers <- round(answers)

  good.cbs.crop <- good.cbs[1:round(length(good.cbs) * (1 - offset.l))]
  bad.cbs.crop <- intersect(bad.cbs, names(answers))
  bad.cbs.crop <- bad.cbs.crop[round(length(bad.cbs.crop) * offset.r):length(bad.cbs.crop)]

  cv.data <- data[c(good.cbs.crop, bad.cbs.crop),]
  cv.answers <- answers[c(good.cbs.crop, bad.cbs.crop)]

  test.force <- list(x=data[intermediate_cbs,], y=answers[intermediate_cbs])

  r <- KFoldCV(cv.data, cv.answers, clf$train, clf$predict, k=10, stratify = F,
               measure = c('sensitivity', 'specifisity'), test.force = test.force, mc.cores=mc.cores)
  return(r)
}
