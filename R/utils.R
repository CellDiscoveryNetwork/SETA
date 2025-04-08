#' Extract Taxonomic Counts from Various Single Cell Objects
#'
#' Given a \code{SingleCellExperiment}, \code{Seurat}, or
#' long-form \code{data.frame}, creates a type-by-sample matrix of cell counts.
#' Users can specify the column names for cell types, samples, and barcodes.
#'
#' @param obj Either a \code{SingleCellExperiment}, a \code{Seurat} object, or
#'   a \code{data.frame}. This function requires specific metadata or columns
#'   to function correctly.
#' @param cell_type_col The name of the column representing cell types.
#'   Default is "type".
#' @param sample_col The name of the column representing sample identifiers.
#'   Default is "sample".
#' @param bc_col The name of the column representing barcodes (only needed for
#'   \code{data.frame} input). Default is "bc".
#'
#' @return A matrix whose rows are samples and whose columns are cell types,
#'   with entries of the count of unique barcodes per type-sample combination.
#'
#' @details
#' \itemize{
#'   \item \strong{SingleCellExperiment}: Reads \code{colData(obj)}
#'   for the specified \code{cell_type_col} and \code{sample_col}.
#'   \item \strong{Seurat}: Uses \code{obj@meta.data} for the specified
#'   \code{cell_type_col} and \code{sample_col}.
#'   \item \strong{data.frame}: Counts cells per specified \code{cell_type_col}
#'   and \code{sample_col}.
#' }
#' If the specified columns are missing, an error is thrown.
#'
#' @examples
#' # For a data.frame with custom column names:
#' set.seed(687)
#' df <- data.frame(
#'   barcode = paste0("cell", 1:10),
#'   cellType = sample(c("Tcell", "Bcell"), 10, TRUE),
#'   sampleID = sample(c("sample1","sample2"), 10, TRUE)
#' )
#' cmat <- setaCounts(df,
#'                    cell_type_col = "cellType",
#'                    sample_col = "sampleID",
#'                    bc_col = "barcode")
#' print(cmat)
#'
#'
#' @export
setaCounts <- function(obj,
                       cell_type_col = "type",
                       sample_col = "sample",
                       bc_col = "bc") {
  # SingleCellExperiment import
  if ("SingleCellExperiment" %in% class(obj)) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' must be installed.", call. = FALSE)
    }
    coldata <- SingleCellExperiment::colData(obj)
    if (!all(c(cell_type_col, sample_col) %in% colnames(coldata))) {
      stop(sprintf("colData(obj) must contain '%s' and '%s' columns.",
                   cell_type_col,
                   sample_col),
           call. = FALSE)
    }
    return(as.matrix(table(coldata[[sample_col]], coldata[[cell_type_col]])))
  }

  # Seurat import
  if ("Seurat" %in% class(obj)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' must be installed.", call. = FALSE)
    }
    meta <- obj@meta.data
    if (!all(c(cell_type_col, sample_col) %in% colnames(meta))) {
      stop(sprintf("obj@meta.data must contain '%s' and '%s' columns.",
                   cell_type_col,
                   sample_col),
           call. = FALSE)
    }
    return(as.matrix(table(meta[[sample_col]], meta[[cell_type_col]])))
  }

  # Long-form data.frame import
  if (is.data.frame(obj)) {
    requiredCols <- c(bc_col, cell_type_col, sample_col)
    if (!all(requiredCols %in% colnames(obj))) {
      stop(sprintf("data.frame must have columns: '%s', '%s', '%s' at minimum.",
                   bc_col,
                   cell_type_col,
                   sample_col),
           call. = FALSE)
    }
    # Deduplicate by barcode so each cell is counted once
    df <- unique(obj[, c(bc_col, cell_type_col, sample_col)])
    return(as.matrix(table(df[[sample_col]], df[[cell_type_col]])))
  }
  stop("Unsupported object type.")
}
