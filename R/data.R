#' @rdname data
#' @name data
#' @aliases mockSC mockSP getMGS
#' @title Synthetic single-cell, mixture and marker data
#'
#' @description
#' \code{mockSC} is designed to generate synthetic single-cell
#'   These data are not meant to represent biologically
#'   meaningful use-cases, but are solely intended for use in examples, for
#'   unit-testing, and to demonstrate \code{SETA}'s general functionality.
#'
#' @param ng,nc,nt,ns,nb integer scalar specifying the number
#'   of genes, cells, types (groups), number of samples and number of batches
#'   to simulate.
#'
#' @return
#' \itemize{
#' \item{\code{mockSC} returns a \code{SeuratObject}
#'   with rows = genes, columns = single cells, and cell metadata
#'   (\code{colData}) column \code{type} containing group identifiers.}
#'
#' @examples
#' sce <- mockSC()
#' df_long <- mockLong()
#' df_count <- mockCount(df_long)
#' df_count2 <- mockCount(sce@meta.data)
NULL

#' @rdname data
#' @importFrom SeuratObject CreateSeuratObject DefaultAssay
#' @importFrom Seurat NormalizeData ScaleData RunPCA
#' @importFrom stats rnbinom runif aggregate
#' @importFrom Matrix Matrix
#' @export

mockSC <- function(ng = 200, nc = 50, nt = 3, ns = 4, nb = 2) {
    z <- lapply(seq_len(nt), \(t) {
        # mu parameter size for rnbinom
        ms <- 2^runif(ng, 2, 10)
        # dispersion parameter for rnbinom
        ds <- 0.5 + 100 / ms
        y <- rnbinom(ng * nc, mu = ms, size = 1 / ds)
        y <- matrix(y, nrow = ng, ncol = nc)
        y <- Matrix(y, sparse = TRUE)
        dimnames(y) <- list(
            paste0("gene", seq_len(ng)),
            paste0("cell", seq_len(nc))
        )
        x <- CreateSeuratObject(counts = y)
        x$type <- factor(
            paste0("type", t),
            paste0("type", seq_len(nt))
        )
        x$batch <- factor(
            sample(paste0("Batch", seq_len(nb)),
                   nc,
                   replace = TRUE),
            paste0("Batch", seq_len(nb))
        )
        x$sample <- factor(
            sample(paste0("Sample", seq_len(ns)),
                   nc,
                   replace = TRUE),
            paste0("Sample", seq_len(ns))
        )
        return(x)
    })
    # Merge the individual cell types
    se <- merge(z[[1]], z[2:nt])

    # Preprocess data
    se <- NormalizeData(se, verbose = FALSE)
    se <- FindVariableFeatures(se, verbose = FALSE)
    se <- ScaleData(se, verbose = FALSE)
    se <- RunPCA(se, npcs = 5, verbose = FALSE)
    # Placeholder for @misc slot modifications if needed for tests
    se@misc$pvclust <- list()

    se
}

mockLong <- function(nc = 500, nt = 3, ns = 4, nb = 2) {

    data.frame(
        bc = paste0("cell", seq_len(nc)),
        type = sample(paste0("type", seq_len(nt)), nc, replace = TRUE),
        batch = sample(paste0("batch", seq_len(nb)), nc, replace = TRUE),
        sample = sample(paste0("sample", seq_len(ns)), nc, replace = TRUE)
    )

}

mockCount <- function(df) {
    aggregate(bc ~ type + sample + batch, data = mockLong(), FUN = length)
}


