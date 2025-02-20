#' @rdname data
#' @name data
#' @aliases mockSC mockSP getMGS
#' @title Synthetic single-cell, mixture and marker data
#'
#' @description
#' \code{mockSC} is designed to generate synthetic single-cell
#' These data are not meant to represent biologically
#' meaningful use-cases, but are solely intended for use in examples, for
#' unit-testing, and to demonstrate \code{SETA}'s general functionality.
#'
#' @param ng,nc,nt integer scalar specifying the number
#'   of genes, cells, types (groups) to simulate.
#'
#' @return
#' \itemize{
#' \item{\code{mockSC} returns a \code{SingleCellExperiment}
#'   with rows = genes, columns = single cells, and cell metadata
#'   (\code{colData}) column \code{type} containing group identifiers.}
#'
#' @examples
#' sce <- mockSC()
NULL

#' @rdname data
#' @importFrom SeuratObject CreateSeuratObject DefaultAssay
#' @importFrom Seurat NormalizeData ScaleData RunPCA
#' @importFrom stats rnbinom runif
#' @importFrom Matrix Matrix
#' @export

mockSC <- function(ng = 200, nc = 50, nt = 3) {
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
            sample(c("Batch1", "Batch2"),
                   nc,
                   replace = TRUE),
            c("Batch1", "Batch2")
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
