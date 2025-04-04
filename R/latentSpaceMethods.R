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

  if(method == "PCA") {
    pca_res <- stats::prcomp(x, center = TRUE, scale. = FALSE)
    out$latentSpace <- as.data.frame(pca_res$x[, seq_len(dims), drop = FALSE])
    out$loadings <- as.data.frame(
        pca_res$rotation[, seq_len(dims), drop = FALSE]
        )
    ve <- pca_res$sdev^2 / sum(pca_res$sdev^2)
    out$varExplained <- ve[seq_len(dims)]
  }

  if(method == "PCoA") {
    d <- stats::dist(x)
    cmd <- stats::cmdscale(d, k = dims, eig = TRUE)
    out$latentSpace <- as.data.frame(cmd$points)
    out$loadings <- NA
    eig_vals <- cmd$eig[cmd$eig > 0]
    out$varExplained <- eig_vals / sum(eig_vals)
  }

  if(method == "NMDS") {
    requireNamespace("MASS", quietly = TRUE)
    d <- stats::dist(x)
    nmds <- MASS::isoMDS(d, k = dims)
    out$latentSpace <- as.data.frame(nmds$points)
    out$loadings <- NA
    out$varExplained <- data.frame(Stress = nmds$stress)
  }

  out
}