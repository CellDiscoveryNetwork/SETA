mockSC <- function(ng = 200, nc = 50, nt = 3, ns = 4, nb = 2) {
    if(!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use mockSC.")
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Matrix package must be installed.")
    }
    z <- lapply(seq_len(nt), \(t) {
        # mu parameter size for rnbinom
        ms <- 2^runif(ng, 2, 10)
        # dispersion parameter for rnbinom
        ds <- 0.5 + 100 / ms
        y <- stats::rnbinom(ng * nc, mu = ms, size = 1 / ds)
        y <- matrix(y, nrow = ng, ncol = nc)
        y <- Matrix::Matrix(y, sparse = TRUE)
        dimnames(y) <- list(
            paste0("gene", seq_len(ng)),
            paste0("cell", seq_len(nc))
        )
        x <- Seurat::CreateSeuratObject(counts = y)
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
    se <- Seurat::NormalizeData(se, verbose = FALSE)
    se <- Seurat::FindVariableFeatures(se, verbose = FALSE)
    se <- Seurat::ScaleData(se, verbose = FALSE)
    se <- Seurat::RunPCA(se, npcs = 5, verbose = FALSE)
    # Placeholder for @misc slot modifications if needed for tests
    se@misc$pvclust <- list()

    se
}

mockLong <- function(
    nc = 500,  # number of cells
    nt = 3,    # number of types
    ns = 4,    # number of samples
    nb = 2,    # number of batches
    useBatch = TRUE
    ) {
    df <- data.frame(
        bc    = paste0("cell", seq_len(nc)),
        type  = sample(paste0("type",   seq_len(nt)), nc, replace = TRUE),
        sample= sample(paste0("sample", seq_len(ns)), nc, replace = TRUE)
    )
    if (useBatch) {
        df$batch <- sample(paste0("batch", seq_len(nb)), nc, replace = TRUE)
    }
    df
}


mockCount <- function(df = mockLong()) {
    # If 'batch' exists, aggregate by (type, sample, batch);
    # otherwise, just (type, sample).
    groupVars <- c("type", "sample")
    if ("batch" %in% colnames(df)) {
        groupVars <- c(groupVars, "batch")
    }
    formulaStr <- paste0("bc ~ ", paste(groupVars, collapse = " + "))
    aggregate(
        stats::as.formula(formulaStr),
        data = df,
        FUN  = length
    )
}


mockSCE <- function(
    nc = 500,
    nt = 3,
    ns = 4,
    nb = 2,
    useBatch = TRUE
) {
    if(!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'SingleCellExperiment' must be installed to use mockSCE.")
    }
    df <- mockLong(nc = nc, nt = nt, ns = ns, nb = nb, useBatch = useBatch)
    bc <- df$bc
    mat <- matrix(stats::rpois(nrow(df) * 20, lambda = 5), nrow = 20)
    rownames(mat) <- paste0("gene", seq_len(20))
    colnames(mat) <- bc
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = mat),
        colData = df
    )
}