#' Centered Log-Ratio Transform
#'
#' Applies a CLR transform to a matrix of counts.
#' This transform is commonly used in compositional data analysis (CoDA)
#' to project counts into a log-ratio space.
#'
#' @param counts A matrix of counts (rows = features, columns = samples).
#' @param pseudocount Numeric.
#'        Added to all entries to avoid taking log(0). Default is 1.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{method}{A string indicating the transform ("CLR").}
#'   \item{counts}{A matrix of the same dimensions
#'                 as the input after CLR transform.}
#' }
#'
#' @details
#' The CLR transform is defined as \eqn{\log(x / g(x))}
#' where \eqn{g(x)} is the geometric mean
#' of each row (sample) in the log scale.
#' A pseudocount helps avoid log(0) - default is 1, as scRNA data can be sparse.
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society.
#'       Series B (Methodological)}, 44(2), 139-177.
#'
#' @examples
#' mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
#' out <- setaCLR(mat, pseudocount=0)
#' out$counts
#'
#' @export
setaCLR <- function(counts, pseudocount = 1) {
  # Basic sanity check
  if(!is.matrix(counts)) stop("'counts' must be a matrix.")

  # Centered log-ratio transform
  # Pseudocount = 1 to avoid issues with zeros
  counts <- counts + pseudocount
  gm <- exp(rowMeans(log(counts)))
  clr_mat <- log(counts / gm)
  
  # Restore column names
  colnames(clr_mat) <- colnames(counts)
  # Return a list with the transformed matrix and metadata
  list(
    method = "CLR",
    counts = clr_mat
      )
}

#' Isometric Log-Ratio Transform
#'
#' Applies an ILR transform to a matrix of counts, using a Helmert basis by default.
#' Optionally includes a Box-Cox–like step on the log scale.
#'
#' @param counts A matrix of counts.
#' @param boxcox_p Numeric. If nonzero, applies a Box-Cox–type transform to the log-values.
#' @param taxTree Currently unused. Reserved for future taxonomic-balance approaches.
#' @param pseudocount Numeric. Pseudocount to avoid log(0). Default is 1.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{A string indicating ILR with a Helmert basis (potentially noting \code{boxcox_p}).}
#'   \item{counts}{A matrix of ILR-transformed values with \code{ncol(counts) - 1} columns.}
#' }
#'
#' @details
#' The ILR transform is a key tool in compositional data analysis. By default, it
#' uses a Helmert contrast matrix. The parameter \code{boxcox_p} allows an
#' additional transform on the log-values, as described by whuber on
#' \url{https://stats.stackexchange.com/questions/259208/how-to-perform-isometric-log-ratio-transformation}.
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 44(2), 139-177.
#'
#' @examples
#' mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
#' out <- setaILR(mat, boxcox_p=0)
#' out$counts
#'
#' @export
setaILR <- function(counts, boxcox_p = 0, taxTree = NULL, pseudocount = 1) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (!is.null(taxTree)) {
    message("A taxTree was provided but is not yet used. Defaulting to Helmert basis.")
  }
  # Add pseudocount to avoid log(0)
  counts <- counts + pseudocount
  # Take the log
  y <- log(counts)
  # Optional Box-Cox–like step on the log-values if p != 0
  # This is suggested by whuber:
  # https://stats.stackexchange.com/questions/259208/how-to-perform-isometric-log-ratio-transformation
  if (boxcox_p != 0) {
    y <- (exp(boxcox_p * y) - 1) / boxcox_p
  }
  
  # Subtract row means to recenter
  y <- y - rowMeans(y, na.rm = TRUE)
  
  # Build the Helmert basis
  k <- ncol(y)
  if (k < 2) stop("ILR requires at least 2 columns (cell types).")
  
  # contr.helmert(k) yields a k x (k-1) matrix; 
  # we transpose and normalize each row by sqrt((2:k)*(2:k-1)).
  # This ensures an orthonormal basis for ILR rotation.
  H <- stats::contr.helmert(k)          # k x (k-1)
  H <- t(H) / sqrt((2:k)*(2:k-1))       # (k-1) x k
  
  # For interpretability, we rename the resulting columns.
  if (!is.null(colnames(counts)) && k > 1) {
    colnames(y)[-1] <- paste0(colnames(y)[-1], ".ILR")
  }
  
  # Apply the Helmert rotation
  ilr_vals <- y %*% t(H)
  
  # Return a list with metadata
  list(
    method = paste0("ILR_Helmert (boxcox_p=", boxcox_p, ")"),
    counts = ilr_vals
  )
}

#' Additive Log-Ratio Transform
#'
#' Applies an ALR transform to a matrix of counts using a specified reference column.
#'
#' @param counts A matrix of counts.
#' @param ref The reference column, specified either by name or by index.
#' @param pseudocount Numeric. Added to avoid log(0). Default is 1.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{A string noting \code{ALR_ref=<ref>}.}
#'   \item{counts}{A matrix of dimension \code{nrow(counts) x (ncol(counts) - 1)}.}
#' }
#'
#' @details
#' ALR transforms the data by taking \eqn{\log(x_i / x_{ref})} for each column
#' \eqn{i} except the reference column. This is another way to map data from
#' the simplex to a Euclidean space in compositional data analysis.
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 44(2), 139-177.
#'
#' @examples
#' mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
#' colnames(mat) <- c("A", "B")
#' out <- setaALR(mat, ref="A", pseudocount=0)
#' out$counts
#'
#' @export
setaALR <- function(counts, ref, pseudocount = 1) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (missing(ref)) stop("Please specify a reference cell type or column.")
  if (is.character(ref)) {
    if (!(ref %in% colnames(counts))) {
      stop("Reference celltype not found in column names.")
    }
    refCol <- which(colnames(counts) == ref)
  } else if (is.numeric(ref)) {
    if (ref < 1 || ref > ncol(counts)) {
        stop("Reference celltype index out of range.")
    }
    refCol <- ref
  } else {
    stop("'ref' must be a character or numeric.")
  }
  counts <- counts + pseudocount
  refVec <- counts[, refCol]
  alrCounts <- log(counts / refVec)
  alrCounts <- alrCounts[, -refCol, drop = FALSE]
  refName <- if (is.character(ref)) ref else colnames(counts)[ref]
  list(method = paste0("ALR_ref=", refName), counts = alrCounts)
}

#' Percentage Transform
#'
#' Converts columns (samples) of a counts matrix to percentages of their respective column sums.
#'
#' @param counts A matrix of counts.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{\code{"percent"}.}
#'   \item{counts}{A matrix of the same dimension, with each column summing to 100.}
#' }
#'
#' @details
#' Simple re-scaling for compositional-like interpretation in percentages.
#' Useful for simplified Wilcoxon rank comparisons and such
#'
#' @examples
#' mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
#' out <- setaPercent(mat)
#' out$counts
#'
#' @export
setaPercent <- function(counts) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  # Example: compute % by columns (assuming columns = samples)
  col_sums <- colSums(counts)
  pct_counts <- t(t(counts) / col_sums * 100)
  list(method = "percent", counts = pct_counts)
}

#' log2(CPM) Transform
#'
#' Computes log2-based counts-per-million (CPM) for each column. Optionally uses provided size factors.
#'
#' @param counts A matrix of counts.
#' @param pseudocount Numeric. Added to counts to avoid \code{log2(0)}. Default is 1.
#' @param size_factors Optional numeric vector. If \code{NULL}, uses the column sums.
#' @param scale_factor Numeric. The scaling factor for "per million" style. Default is \code{1e6}.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{\code{"logCPM"}.}
#'   \item{counts}{A matrix of the same dimension, containing log2(CPM + pseudocount).}
#' }
#'
#' @details
#' A common RNA-seq transform is log2(CPM + 1). This variant allows adjusting the
#' pseudocount, size factors, and overall scale factor.
#'
#' @examples
#' mat <- matrix(1:4, nrow=2)
#' out <- setaLogCPM(mat)
#' out$counts
#'
#' @export
setaLogCPM <- function(counts, pseudocount = 1, size_factors = NULL,
                       scale_factor = 1e6) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (is.null(size_factors)) size_factors <- colSums(counts)
  # Add pseudocount, divide by size factor,
  # then multiply by size factor
  cpm <- sweep(counts + pseudocount, 2, 
               size_factors + pseudocount, FUN = "/") * scale_factor
  log_cpm <- log2(cpm)
  list(method = "logCPM", counts = log_cpm)
}

#' Wrapper for Common Compositional Transforms
#'
#' A convenience function that dispatches to one of the transforms:
#' CLR, ALR, ILR, percent, or logCPM.
#'
#' @param counts A matrix of counts.
#' @param method Which transform to apply. One of \code{"CLR"}, \code{"ALR"},
#'   \code{"ILR"}, \code{"percent"}, or \code{"logCPM"}.
#' @param ref Reference column (only used if \code{method="ALR"}).
#' @param taxTree Optional tree for ILR (not yet implemented).
#' @param pseudocount Numeric, used by CLR, ALR, ILR, logCPM. Default 1.
#' @param size_factors For logCPM scaling. If \code{NULL}, uses column sums.
#'
#' @return A list as returned by the respective transform function.
#'
#' @examples
#' mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
#' setaTransform(mat, method="CLR")
#' setaTransform(mat, method="percent")
#'
#' @export
setaTransform <- function(
  counts,
  method = c("CLR", "ALR", "ILR", "percent", "logCPM"),
  ref = NULL,
  taxTree = NULL,
  pseudocount = 1,
  size_factors = NULL
) {
  method <- match.arg(method)
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  
  switch(method,
    "CLR" = {
      return(setaCLR(counts, pseudocount = pseudocount))
    },
    "ALR" = {
      return(setaALR(counts, ref = ref, pseudocount = pseudocount))
    },
    "ILR" = {
      return(setaILR(counts, taxTree = taxTree, pseudocount = pseudocount))
    },
    "percent" = {
      return(setaPercent(counts))
    },
    "logCPM" = {
      return(setaLogCPM(counts, pseudocount = pseudocount, 
                        size_factors = size_factors))
    }
  )
}