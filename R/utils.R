#' @useDynLib dropEstAnalysis
NULL

#' @importFrom dplyr %>%
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("dropEstAnalysis", libpath)
}

#' @export
ExpressionMatrixToDataFrame <- function(matrix, umis.per.cb, clusters, rescued.cbs = NULL, filtration.type=NULL,
                                        normalize=TRUE) {
  clusters <- as.factor(clusters)
  if (normalize) {
    matrix <- log10(1e-3 + matrix) %>% apply(2, function(vec) scales::rescale(rank(vec)))
  }
  h.clust.col <- hclust(dist(t(matrix)))
  gene.order <- h.clust.col$labels[h.clust.col$order]
  barcodes.order <- names(clusters)[order(-as.integer(clusters), umis.per.cb[names(clusters)])]

  res <- matrix %>% as.data.frame() %>% tibble::rownames_to_column('Barcode') %>%
    reshape2::melt(id.vars='Barcode', variable.name='Gene', value.name='Expression') %>%
    dplyr::mutate(Gene = factor(Gene, levels = gene.order),
                  Barcode = factor(Barcode, levels = barcodes.order),
                  UmisPerCb = umis.per.cb[as.character(Barcode)],
                  Cluster = clusters[as.character(Barcode)]) %>%
    dplyr::arrange(Gene, Barcode)

  if (!is.null(rescued.cbs)) {
    res$IsRescued <- (res$Barcode %in% rescued.cbs)
  }

  if (!is.null(filtration.type)) {
    res$FiltrationType <- filtration.type[as.character(res$Barcode)]
  }

  return(res)
}

#' @export
FindClusterMarkers <- function(clust2, clust1, srt.obj, max.pval=1e-5) {
  res <- Seurat::FindMarkers(object = srt, ident.1 = clust1, ident.2 = clust2, min.pct = 0.25, only.pos=T) %>%
    tibble::rownames_to_column('gene') %>% filter(p_val_adj < 1e-5) %>% .$gene

  return(res)
}

#' @export
GetOverexpressedGenes <- function(srt, compared.clusters, cluster.markers, genes.from.cluster=50, expression.threshold=0.6) {
  genes <- lapply(cluster.markers, function(x)
    (unlist(x) %>% table() %>% sort(decreasing=T) %>% names())[1:genes.from.cluster]) %>% unlist() %>% unique()
  genes <- genes[!is.na(genes)]

  gene.mask <- lapply(compared.clusters, function(cl)
    Matrix::rowMeans(srt@data[genes, names(srt@ident)[srt@ident == cl]] > 0) > expression.threshold)

  gene.mask <- Reduce(`|`, gene.mask)
  return(names(gene.mask)[gene.mask])
}

#' @export
GetCellsChull <- function(cbs, tsne, chull.quantile=0.95, offset.x=0.5, offset.y=0.5) {
  tsne <- tsne[cbs, ]
  rob.stats <- robustbase::covMcd(tsne)
  dists <- mahalanobis(tsne, center=rob.stats$center, cov=rob.stats$cov)
  cbs <- names(dists)[dists < quantile(dists, chull.quantile)]
  tsne <- tsne[cbs, ]
  res <- tsne[chull(tsne),]
  res <- res + c(offset.x, offset.y) * sign(res - rob.stats$center)
  return(list(chull=res[chull(res), ], cbs=cbs))
}

#' @export
AnnotateClustersByGraph <- function(graph, annotated.clusters, notannotated.cells, max.iter=20, mc.cores=1) {
  getNeighbourCorrs <- function(source.cb, graph) {
    n.cbs <- as.list(igraph::neighbors(graph, source.cb)) %>% names()
    dists <- sapply(n.cbs, function(target.cb) igraph::get.edge.ids(graph, c(source.cb, target.cb)) %>%
                      igraph::edge.attributes(graph=graph) %>% .$weight)
    return(1 - dists[dists < 1])
  }

  corrs <- parallel::mclapply(notannotated.cells, getNeighbourCorrs, graph=graph, mc.cores=mc.cores) %>%
    setNames(notannotated.cells)

  for (i in 1:max.iter) {
    notannotated.clusters <- lapply(corrs, function(corr) corr[names(corr) %in% names(annotated.clusters)]) %>%
      lapply(function(corr) split(corr, annotated.clusters[names(corr)]) %>% sapply(sum, na.rm=T) %>% which.max() %>% names()) %>%
      unlist()

    annotated.clusters[names(notannotated.clusters)] <- notannotated.clusters
    if (length(notannotated.clusters) == length(notannotated.cells))
      break()
  }

  return(annotated.clusters)
}

# Merge
#' @export
MergeComparisonSummary <- function(merge.targets, cell.species, dataset,
                                   filt.merge.types=c('real', 'merge_all', 'poisson')) {
  cbPair <- function(cbs) paste0(cbs, names(cbs))

  wrong.merge.fractions <- lapply(merge.targets, function(mt) mean(cell.species[mt] != cell.species[names(mt)])) %>%
    unlist()
  same.as.real <- merge.targets %>% sapply(function(mt) intersect(cbPair(mt), cbPair(merge.targets$real)) %>% length())

  res <- data.frame(WrongPercent=wrong.merge.fractions, MergesNum=sapply(merge.targets, length), SameAsReal=same.as.real) %>%
    tibble::rownames_to_column('Merge') %>% filter(Merge %in% filt.merge.types) %>%
    dplyr::mutate(SameAsReal = round(100 * SameAsReal / MergesNum, 2) %>% paste0("%"),
           WrongPercent = round(100 * WrongPercent, 2) %>% paste0("%"),
           Dataset=dataset) %>%
    dplyr::select(Dataset, Merge, MergesNum, WrongPercent, SameAsReal) %>%
    dplyr::rename(`Merge type`=Merge, `#Merges`=MergesNum, `Fraction of mixed merges`=WrongPercent,
                  `Similarity to merge with barcodes`=SameAsReal)

  return(res)
}

#' @export
FilterNUmis <- function(reads.per.umi) {
  return(lapply(reads.per.umi, function(rpus) rpus[grep("^[^N]+$", names(rpus))]))
}

#' @export
FillNa <- function(data, value=0) {
  data[is.na(data)] <- value
  return(data)
}

#' @export
Read10xMatrix <- function(path, use.gene.names=FALSE) {
  gene.var <- if (use.gene.names) 'V2' else 'V1'
  mtx <- as(Matrix::readMM(paste0(path, 'matrix.mtx')), 'dgCMatrix')
  colnames(mtx) <- read.table(paste0(path, 'barcodes.tsv'), stringsAsFactors=F)$V1
  rownames(mtx) <- read.table(paste0(path, 'genes.tsv'), stringsAsFactors=F)[[gene.var]]
  return(mtx)
}
