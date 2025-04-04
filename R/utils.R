#' Extract Taxonomic Counts from Various Single Cell Objects
#'
#' Given a \code{SingleCellExperiment}, \code{Seurat}, or long-form \code{data.frame},
#' this function produces a type-by-sample matrix of cell counts. For \code{SingleCellExperiment}
#' or \code{Seurat}, it looks for \code{type} and \code{sample} columns in the object-level metadata.
#' For long-form data frames, it expects columns \code{bc}, \code{type}, and \code{sample}.
#'
#' @param obj Either a \code{SingleCellExperiment}, a \code{Seurat} object, or
#'   a \code{data.frame} with columns \code{bc}, \code{type}, and \code{sample}.
#'
#' @return A matrix whose rows are cell types and whose columns are samples,
#'   with entries giving the count of unique barcodes per type-sample combination.
#'
#' @details
#' \itemize{
#'   \item \strong{SingleCellExperiment}: Reads \code{colData(obj)} for \code{type} and \code{sample}.
#'   \item \strong{Seurat}: Uses \code{obj@meta.data} for \code{type} and \code{sample}.
#'   \item \strong{data.frame}: Duplicates (by \code{bc}) are removed so each cell is counted once.
#' }
#' If \code{type} or \code{sample} columns are missing, an error is thrown.
#'
#' @examples
#' \donttest{
#' # For a data.frame:
#' df <- data.frame(
#'   bc = paste0("cell", 1:10),
#'   type = sample(c("Tcell", "Bcell"), 10, TRUE),
#'   sample = sample(c("sample1","sample2"), 10, TRUE)
#' )
#' cmat <- setaCounts(df)
#' cmat
#' }
#'
#' @export
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