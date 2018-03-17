#' @export
TimeCollisions <- function(holder) {
  t0 <- proc.time();
  umi_dist <- dropestr::GetUmisDistribution(holder$reads_per_umi_per_cell)
  coll_info <- dropestr::FillCollisionsAdjustmentInfo(umi_dist / sum(umi_dist), max(holder$cm_raw@x))
  r <- coll_info[holder$cm_raw@x]
  return((proc.time() - t0)[3])
}

#' @export
TimeUmiErrors <- function(rpu.per.cell, method = 'Bayesian', mc.cores=1) {
  t0 <- proc.time()
  if (method == 'directional') {
    corrected_mtx <- dropestr::CorrectUmiSequenceErrors(rpu.per.cell, method='Classic', mult=2, mc.cores=mc.cores)
  } else {
    corrected_mtx <- dropestr::CorrectUmiSequenceErrors(rpu.per.cell, mc.cores=mc.cores)
  }

  return((proc.time() - t0)[3])
}

#' @export
TimeQuality <- function(holder, ...) {
  t0 <- proc.time()
  dropestr::ScorePipelineCells(holder, ...)
  return((proc.time() - t0)[3])
}
