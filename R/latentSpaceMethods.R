#' Compute a Latent Space from Transformed Counts
#'
#' Given an object produced by one of the
#' \code{seta*} transform functions (e.g., \code{setaCLR}),
#' this function applies a dimension reduction method (PCA, PCoA, or NMDS) to
#' \code{transform_obj$counts}.
#'
#' @param transform_obj A list returned by, e.g., \code{setaCLR}, \code{setaILR}
#'   containing a \code{counts} matrix, where rows are samples
#'   and columns are features (taxa or cell types).
#' @param method A string specifying the dimension reduction method.
#'   One of \code{"PCA"}, \code{"PCoA"}, or \code{"NMDS"}.
#' @param dims Integer number of dimensions to return. Default is 2.
#'
#' @return A list containing:
#' \describe{
#'   \item{method}{The chosen latent space method.}
#'   \item{latentSpace}{A data frame of coordinates in the chosen latent space
#'                      with \code{dims} columns.}
#'   \item{loadings}{For PCA, the loadings matrix. Otherwise \code{NA}.}
#'   \item{varExplained}{Variance explained (for PCA or PCoA)
#'                       or stress (for NMDS).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{PCA}: Uses \code{stats::prcomp}
#'                       on the rows of \code{transform_obj$counts}.
#'   \item \strong{PCoA}: Computes a distance matrix
#'                        via \code{stats::dist}, then
#'         applies classical multidimensional scaling (\code{stats::cmdscale}).
#'   \item \strong{NMDS}: Uses \code{MASS::isoMDS} to compute
#'                        non-metric MDS from the distance matrix.
#' }
#' Each method returns a data frame of coordinates in \code{latentSpace},
#' plus additional information specific to that method.
#'
#' @examples
#' set.seed(687)
#' mat <- matrix(rpois(20, lambda=5), nrow=4)  # small 4x5 matrix
#' colnames(mat) <- paste0("C", 1:5)
#' clr_out <- setaCLR(mat)
#' latent_pca <- setaLatent(clr_out, method="PCA", dims=2)
#' latent_pca$latentSpace
#'
#' @importFrom stats prcomp dist cmdscale
#' @importFrom MASS isoMDS
#' @export
setaLatent <- function(transform_obj,
                       method = c("PCA", "PCoA", "NMDS"),
                       dims = 2) {
    method <- match.arg(method)
    x <- transform_obj$counts
    if (!is.matrix(x)) stop("'transform_obj$counts' must be a matrix.")
    
    out <- list(method = method,
                latentSpace = NULL,
                loadings = NULL,
                varExplained = NULL)
    
    if (method == "PCA") {
        pca_res <- prcomp(x, center = TRUE, scale. = FALSE)
        out$latentSpace <- as.data.frame(pca_res$x[, seq_len(dims), drop = FALSE])
        out$loadings <- as.data.frame(
            pca_res$rotation[, seq_len(dims), drop = FALSE]
        )
        ve <- pca_res$sdev^2 / sum(pca_res$sdev^2)
        out$varExplained <- ve[seq_len(dims)]
    }
    
    if (method == "PCoA") {
        d <- dist(x)
        cmd <- cmdscale(d, k = dims, eig = TRUE)
        out$latentSpace <- as.data.frame(cmd$points)
        out$loadings <- NA
        eig_vals <- cmd$eig[cmd$eig > 0]
        out$varExplained <- eig_vals / sum(eig_vals)
    }
    
    if (method == "NMDS") {
        d <- dist(x)
        nmds <- isoMDS(d, k = dims)
        out$latentSpace <- as.data.frame(nmds$points)
        out$loadings <- NA
        out$varExplained <- data.frame(Stress = nmds$stress)
    }
    
    out
}
