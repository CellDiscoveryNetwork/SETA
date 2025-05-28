#' Compute Distance Matrix between Samples
#'
#' Calculates a pairwise distance matrix between samples based on user-specified
#' or default (\code{"euclidean"}) distance metrics. If used on CLR-transformed
#' data, the default Euclidean distance is the \emph{Aitchison distance},
#' which is commonly used in compositional data analysis (CoDA).
#'
#' @param transformed_counts Numeric matrix: rows as samples and columns as taxa
#'        (e.g., output from \code{setaCLR}, \code{setaTransform}, etc.).
#' @param method Character. Distance metric for \code{\link[stats]{dist}}.
#'        Default: \code{"euclidean"} See \code{\link[stats]{dist}} for options.
#'
#' @return A long-form \code{data.frame} with three columns:
#' \describe{
#'   \item{from}{Sample ID of the first sample in the pairwise comparison.}
#'   \item{to}{Sample ID of the second sample in the pairwise comparison.}
#'   \item{distance}{Numeric distance between the two samples.}
#' }
#'
#' @details
#' This function calculates distances \strong{between samples}.
#'
#' Output is a long-form structure convenient to merge with sample-level
#' metadata using \code{\link[base]{merge}} or \code{\link[dplyr]{left_join}}.
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)},
#' 44(2), 139-177.
#'
#' @examples
#' # Example CLR transformed data (2 samples, 3 taxa)
#' mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("Taxon1", "Taxon2", "Taxon3")
#' rownames(mat) <- c("SampleA", "SampleB")
#'
#' clr_mat <- setaCLR(mat, pseudocount = 0)
#'
#' # Calculate Euclidean (Aitchison) distance
#' dist_df <- setaDistances(clr_mat)
#' print(head(dist_df))
#'
#' @importFrom stats dist function
#' @export
setaDistances <- function(transformed_counts, method = "euclidean") {
  if (!is.matrix(transformed_counts$counts)) {
    stop("'transformed_counts' must be a numeric matrix 
            from setaTransform or other seta methods
            with samples in rows and taxa in columns.")
  }

  dist_mat <- dist(transformed_counts$counts, method = method)

  # Convert distance matrix to a long-form dataframe
  dist_df <- as.data.frame(as.table(as.matrix(dist_mat)))
  colnames(dist_df) <- c("from", "to", "distance")

  # Remove self-distances and duplicated pairs
  dist_df <- subset(dist_df, as.character(from) < as.character(to))

  rownames(dist_df) <- NULL
  dist_df
}
