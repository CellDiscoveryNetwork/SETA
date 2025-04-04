setaCounts <- function(obj) {
  # SingleCellExperiment import
  if ("SingleCellExperiment" %in% class(obj)) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' must be installed.", call. = FALSE)
    }
    coldata <- SingleCellExperiment::colData(obj)
    if (!all(c("type", "sample") %in% colnames(coldata))) {
      stop("colData(obj) must contain 'type' and 'sample' columns.")
    }
    return(as.matrix(table(coldata$type, coldata$sample)))
  }

  # Seurat import
  if ("Seurat" %in% class(obj)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' must be installed.", call. = FALSE)
    }
    meta <- obj@meta.data
    if (!all(c("type", "sample") %in% colnames(meta))) {
      stop("obj@meta.data must contain 'type' and 'sample' columns.")
    }
    return(as.matrix(table(meta$type, meta$sample)))
  }

  # Long-form data.frame import
  if (is.data.frame(obj)) {
    requiredCols <- c("bc", "type", "sample")
    if (!all(requiredCols %in% colnames(obj))) {
      stop("data.frame must have columns: 'bc', 'type', 'sample' at minimum.")
    }
    # Deduplicate by barcode so each cell is counted once
    df <- unique(obj[, c("bc", "type", "sample")])
    return(as.matrix(table(df$type, df$sample)))
  }
  stop("Unsupported object type.")
}