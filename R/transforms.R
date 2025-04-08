#' Centered Log-Ratio (CLR) Transform
#'
#' Applies a CLR transform to a matrix of counts. In this version, samples are assumed to be
#' in rows and taxa (cell types) in columns. For each sample, the transform computes
#' \eqn{\mathrm{CLR}(x)_i = \log\left(\frac{x_i + \text{pseudocount}}{g(x + \text{pseudocount})}\right)},
#' where \eqn{g(x + \text{pseudocount})} is the geometric mean of the row.
#'
#' @param counts A numeric matrix of counts with rows as samples and columns as taxa.
#' @param pseudocount Numeric. Added to all entries to avoid \code{log(0)}. Default is 1.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{A string indicating the transform ("CLR").}
#'   \item{counts}{A matrix of the same dimensions as the input after CLR transform.}
#' }
#'
#' @details
#' The CLR transform is defined sample-wise as:
#' \deqn{\mathrm{CLR}(x)_{ij} = \log\left(\frac{x_{ij} + \text{pseudocount}}{g_i}\right),}{
#' \log\left(\frac{x_{ij} + \text{pseudocount}}{g_i}\right)}
#' where \eqn{g_i = \exp\left(\frac{1}{p}\sum_{j=1}^{p}\log(x_{ij} + \text{pseudocount})\right)} for sample \(i\)
#' and \(p\) is the number of taxa.
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 44(2), 139-177.
#'
#' @examples
#' # Example matrix with 2 samples and 2 taxa:
#' mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("Taxon1", "Taxon2")
#' out <- setaCLR(mat, pseudocount = 0)
#' out$counts
#'
#' @export
setaCLR <- function(counts, pseudocount = 1) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  counts <- counts + pseudocount
  # Compute geometric mean for each sample (row)
  gm <- exp(apply(log(counts), 1, mean))
  # Subtract log(geometric mean) from each log-transformed element in each row
  clr_mat <- sweep(log(counts), 1, log(gm), FUN = "-")
  # Restore column names
  colnames(clr_mat) <- colnames(counts)
  list(method = "CLR", counts = clr_mat)
}

#' Isometric Log-Ratio (ILR) Transform
#'
#' Applies the ILR transform to a counts matrix using an orthonormal Helmert basis.
#' For each sample (row), the data are log-transformed (with an optional Box-Cox–like transformation),
#' and then projected onto an orthonormal Helmert basis to reduce dimensionality by one.
#'
#' @param counts A numeric matrix of counts with rows as samples and columns as taxa or cell types.
#' @param boxcox_p Numeric. If nonzero, a Box-Cox–type transform is applied to the log-values.
#'   Default is 0 (no Box-Cox transformation).
#' @param taxTree Currently unused. Reserved for future taxonomic-balance approaches.
#' @param pseudocount Numeric. Added to avoid \code{log(0)}. Default is 1.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{A string indicating the ILR transform. If \code{boxcox_p} is nonzero,
#'                  the value is indicated in the method string.}
#'   \item{counts}{A matrix of ILR-transformed values with \code{ncol(counts) - 1} columns
#'                 and the same number of rows (samples) as the input.}
#' }
#'
#' @details
#' The ILR transform is computed as follows:
#' \enumerate{
#'   \item Add a pseudocount and take the natural logarithm:
#'     \deqn{y = \log(x + \text{pseudocount})}
#'   \item If \code{boxcox_p != 0}, apply the Box-Cox–like transform:
#'     \deqn{y = \frac{\exp(p \, y) - 1}{p}}
#'   \item Project the log-transformed data onto an orthonormal Helmert basis computed via QR decomposition.
#' }
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 44(2), 139-177.
#'
#' @examples
#' # Example matrix: rows are samples, columns are cell types.
#' mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("A", "B")
#' # ILR transformation reduces the dimension by 1.
#' out <- setaILR(mat, boxcox_p = 0, pseudocount = 1)
#' out$counts
#'
#' @export
setaILR <- function(counts, boxcox_p = 0, taxTree = NULL, pseudocount = 1) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (!is.null(taxTree)) {
    message("A taxTree was provided but is not yet supported. Defaulting to Helmert basis.")
  }
  counts <- counts + pseudocount
  log_x <- log(counts)
  if (boxcox_p != 0) {
    log_x <- (exp(boxcox_p * log_x) - 1) / boxcox_p
  }
  n_taxa <- ncol(counts)
  if (n_taxa < 2) stop("ILR requires at least 2 taxa (columns).")

  # Compute orthonormal Helmert basis
  H <- stats::contr.helmert(n_taxa)
  H <- qr.Q(qr(H))  # orthonormalize via QR decomposition
  # Directly project log-values onto the orthonormal Helmert basis
  ilr_mat <- log_x %*% H
  list(
    method = paste0("ILR_Helmert", if (boxcox_p != 0) paste0(" (boxcox_p=", boxcox_p, ")") else ""),
    counts = ilr_mat
  )
}

#' Additive Log-Ratio (ALR) Transform
#'
#' Applies the ALR transform to a matrix of counts using a specified reference taxon.
#' Samples are in rows and taxa in columns. 
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#' @param ref Either the reference taxon name (a character string, which must appear in \code{colnames(counts)})
#'   or the column index of the reference.
#' @param pseudocount Numeric. Added to every count to avoid \code{log(0)}. Default is 1.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{A string indicating the ALR transform with the reference taxon.}
#'   \item{counts}{A matrix with the same number of rows and \eqn{(\text{n_taxa} - 1)} columns.}
#' }
#'
#' @details
#' Applies the ALR transform to a matrix of counts using a specified reference taxon.
#' Samples are in rows and taxa in columns. For each sample, the transform computes
#' \deqn{\mathrm{ALR}(x)_i = \log\left(\frac{x_i + \text{pseudocount}}{x_{ref} + \text{pseudocount}}\right)}{
#' \log\left(\frac{x_i + \text{pseudocount}}{x_{ref} + \text{pseudocount}}\right)}
#' for all taxa \(i\) except the reference.
#'
#' @examples
#' # Example with 2 samples and 2 taxa:
#' mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("TaxonA", "TaxonB")
#' # Using TaxonA as the reference.
#' out <- setaALR(mat, ref = "TaxonA", pseudocount = 0)
#' out$counts
#'
#' @export
setaALR <- function(counts, ref, pseudocount = 1) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (missing(ref)) stop("Please specify a reference taxon (by name or index).")
  if (is.character(ref)) {
    if (!(ref %in% colnames(counts))) {
      stop("Reference taxon not found in colnames(counts).")
    }
    refCol <- which(colnames(counts) == ref)
  } else if (is.numeric(ref)) {
    if (ref < 1 || ref > ncol(counts)) stop("Reference taxon index out of range.")
    refCol <- ref
  } else {
    stop("'ref' must be a character or numeric.")
  }
  counts <- counts + pseudocount
  # For each sample, subtract log(value at reference taxon) from log(counts)
  alr_mat <- sweep(log(counts), 1, log(counts[, refCol]), FUN = "-")
  # Remove the reference taxon column from output
  alr_mat <- alr_mat[, -refCol, drop = FALSE]
  list(method = paste0("ALR_ref=", if (is.character(ref)) ref else colnames(counts)[refCol]),
       counts = alr_mat)
}

#' Percentage Transform
#'
#' Converts each row (sample) of a counts matrix to percentages of its row sum.
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{The string \code{"percent"}.}
#'   \item{counts}{A matrix of the same dimensions as \code{counts}, where each row sums to 100.}
#' }
#'
#' @details
#' Useful for simplified comparisons and as an input to non-parametric tests.
#'
#' @examples
#' mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
#' out <- setaPercent(mat)
#' out$counts
#'
#' @export
setaPercent <- function(counts) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  pct_counts <- sweep(counts, 1, rowSums(counts), FUN = "/") * 100
  list(method = "percent", counts = pct_counts)
}

#' log2(CPM) Transform
#'
#' Computes the log2 counts-per-million (CPM) for each sample.
#' Samples are in rows and taxa in columns.
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#' @param pseudocount Numeric. Added to counts to avoid \code{log2(0)}. Default is 1.
#' @param size_factors Optional numeric vector of library sizes for each sample.
#'   If \code{NULL}, the row sums are used.
#' @param scale_factor Numeric. The scaling factor, typically 1e6 for CPM. Default is 1e6.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{The string \code{"logCPM"}.}
#'   \item{counts}{A matrix of the same dimensions with log2-transformed CPM values.}
#' }
#'
#' @details
#' The transform is defined as:
#' \deqn{\log_2\left(\frac{(x + \text{pseudocount})}{\text{library size}} \times \text{scale_factor}\right),}{
#' \log_2\left(\frac{(x + \text{pseudocount})}{\text{library size}} \times \text{scale_factor}\right)}
#' where the library size is computed per sample.
#' @examples
#' mat <- matrix(c(10, 20, 100, 200), nrow = 2, byrow = TRUE)
#' out <- setaLogCPM(mat, pseudocount = 1)
#' out$counts
#'
#' @export
setaLogCPM <- function(counts, pseudocount = 1, size_factors = NULL, scale_factor = 1e6) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (is.null(size_factors)) size_factors <- rowSums(counts)
  cpm <- sweep(counts + pseudocount, 1, size_factors + pseudocount, FUN = "/") * scale_factor
  log_cpm <- log2(cpm)
  list(method = "logCPM", counts = log_cpm)
}

#' Wrapper for Compositional Transforms
#'
#' A convenience function that dispatches to one of the transforms:
#' CLR, ALR, ILR, percent, or logCPM. Note that the input \code{counts} matrix should have rows as samples
#' and columns as taxa.
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#' @param method A character string specifying which transform to apply. One of
#'   \code{"CLR"}, \code{"ALR"}, \code{"ILR"}, \code{"percent"}, or \code{"logCPM"}.
#' @param ref Reference taxon (only used if \code{method="ALR"}). Can be a taxon name or column index.
#' @param taxTree Optional tree for ILR (not yet implemented).
#' @param pseudocount Numeric, used by CLR, ALR, ILR, and logCPM. Default is 1.
#' @param size_factors For logCPM scaling. If \code{NULL}, uses row sums.
#'
#' @return A list as returned by the corresponding transform function.
#'
#' @examples
#' mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
#' # Apply CLR transform:
#' setaTransform(mat, method = "CLR", pseudocount = 1)
#' # Apply percent transform:
#' setaTransform(mat, method = "percent")
#'
#' @export
setaTransform <- function(counts,
                          method = c("CLR", "ALR", "ILR", "percent", "logCPM"),
                          ref = NULL,
                          taxTree = NULL,
                          pseudocount = 1,
                          size_factors = NULL) {
  method <- match.arg(method)
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  switch(method,
    "CLR"     = setaCLR(counts,
                        pseudocount = pseudocount),
    "ALR"     = setaALR(counts,
                        ref = ref,
                        pseudocount = pseudocount),
    "ILR"     = setaILR(counts,
                        taxTree = taxTree,
                        pseudocount = pseudocount),
    "percent" = setaPercent(counts),
    "logCPM"  = setaLogCPM(counts,
                           pseudocount = pseudocount,
                           size_factors = size_factors)
  )
}