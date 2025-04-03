getCounts <- function(obj) {
  if ("SingleCellExperiment" %in% class(obj)) {
    return(as.matrix(SummarizedExperiment::assay(obj, "counts")))
  }
  if ("Seurat" %in% class(obj)) {
    return(as.matrix(obj@assays$RNA@counts))
  }
  if (is.data.frame(obj) &&
      all(c("bc", "type", "sample", "batch") %in% names(obj))) {
    w <- reshape2::dcast(obj, type ~ sample, value.var = "bc",
                         fun.aggregate = length)
    rownames(w) <- w$type
    w <- as.matrix(w[, -1])
    return(w)
  }
  stop("Unsupported object type.")
}