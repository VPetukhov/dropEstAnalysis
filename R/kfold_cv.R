# library(parallel)
# library(pROC)

GetTraialsSummary <- function(trials.data) {
  r <- t(sapply(trials.data, unlist))
  return(data.frame(mean=apply(r, 2, mean), sd=apply(r, 2, sd)))
}

#' @export
PlotTrialsSummary <- function(trials.data, x.vals, limits=c(0, 1), alpha=0.8) {
  vars.num <- nrow(trials.data[[1]])
  means <- cbind(data.table::rbindlist(lapply(trials.data, function(v) as.list(v$mean))), X=x.vals)
  stds <- cbind(data.table::rbindlist(lapply(trials.data, function(v) as.list(v$sd))), X=x.vals)
  colnames(means)[1:vars.num] <- rownames(trials.data[[1]])
  colnames(stds)[1:vars.num] <- rownames(trials.data[[1]])

  lower.ci <- means
  upper.ci <- means

  lower.ci[1:vars.num] <- lower.ci[1:vars.num] - 3 * stds[1:vars.num]
  upper.ci[1:vars.num] <- upper.ci[1:vars.num] + 3 * stds[1:vars.num]

  if (!is.null(limits)) {
    lower.ci[lower.ci < limits[1]] <- limits[1]
    upper.ci[upper.ci > limits[2]] <- limits[2]
  }

  upper.ci <- reshape2::melt(upper.ci, id.vars = 'X', variable.name = 'Measure', value.name = 'Upper')$Upper
  lower.ci <- reshape2::melt(lower.ci, id.vars = 'X', variable.name = 'Measure', value.name = 'Lower')$Lower

  means.long <- cbind(reshape2::melt(means, id.vars = 'X', variable.name = 'Measure'), upper=upper.ci, lower=lower.ci)
  means.long <- reshape2::melt(means.long, id.vars = c('X', 'Measure'), variable.name = 'Type')

  ggplot2::ggplot(means.long) +
    ggplot2::geom_line(ggplot2::aes(x=X, y=value, color=Measure, linetype=Type), alpha=alpha) +
    ggplot2::scale_linetype_manual(values = c('solid', 'dashed', 'dashed'), labels = c('Value', 'Upper CI', 'Lower CI'))
}

SplitClass <- function(cl.data, k) {
  marks <- cut(sample(seq(1,nrow(cl.data))),breaks=k,labels=FALSE)
  inds <- lapply(1:k, function(i) which(marks == i))
  return(lapply(inds, function(s.inds) cl.data[s.inds,]))
}

SplitToFolds <- function(data, answers, k, stratify=FALSE) {
  answers <- round(answers)
  if (stratify) {
    min.num <- min(table(answers))
    inds <- as.vector(sapply(0:1, function(i) {which(answers == i)[sample.int(n=sum(answers==i), size=min.num)]}))
    data <- data[inds,]
    answers <- answers[inds]
  }

  r <- lapply(0:1, function(cl) SplitClass(data[answers == cl,], k))
  cl1.folds <- r[[1]]
  cl2.folds <- r[[2]]
  res <- lapply(1:k, function(i) list(x=rbind(cl1.folds[[i]], cl2.folds[[i]]),
                                      y=rep(0:1, c(nrow(cl1.folds[[i]]), nrow(cl2.folds[[i]])))))
  names(res) <- paste(1:k)
  return(res)
}

GetTrainFolds <- function(test.fold, folds) {
  k <- length(folds)
  train.x <- do.call('rbind', lapply(folds[paste(setdiff(1:k, test.fold))], function(a) a$x))
  train.y <- as.vector(unlist(sapply(folds[paste(setdiff(1:k, test.fold))], function(a) a$y)))
  return(list(x=train.x, y=train.y))
}

PredictFold <- function(fold, folds, train.func, predict.func, test.force=NULL,
                        measure=c('accuracy', 'mse', 'specifisity', 'sensitivity', 'roc')) {
  train <- GetTrainFolds(fold, folds)
  test <- if (is.null(test.force)) folds[[paste(fold)]] else test.force

  clf <- train.func(train$x, train$y)
  ans <- predict.func(clf, test$x)
  res <- list()
  if ('accuracy' %in% measure) res$accuracy <- mean(round(ans) == round(test$y))
  if ('mse' %in% measure) res$mse <- sqrt(mean((test$y - ans)^2))
  if ('specifisity' %in% measure) res$specifisity <- sum((round(ans) == 0) & (round(test$y) == 0)) / sum(round(test$y) == 0)
  if ('sensitivity' %in% measure) res$sensitivity <- sum((round(ans) == 1) & (round(test$y) == 1)) / sum(round(test$y) == 1)
  if ('roc' %in% measure) res$roc <- pROC::auc(round(test$y), ans)
  return(res)
}

#' @export
KFoldCV <- function(data, answers, train.func, predict.func, k=5, stratify=FALSE, test.force=NULL,
                    measure=c('accuracy', 'mse', 'specifisity', 'sensitivity', 'roc'), mc.cores=NULL) {
  if (is.null(mc.cores)) {
    mc.cores <- if (is.null(options()$mc.cores)) 1 else options()$mc.cores
  }
  folds <- SplitToFolds(data, answers, k, stratify)
  res <- parallel::mclapply(1:k, function(fold)
    PredictFold(fold, folds, train.func, predict.func, measure = measure, test.force=test.force),
    mc.cores=mc.cores, mc.allow.recursive=F)

  return(GetTraialsSummary(res))
}
