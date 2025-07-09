mockSC <- function(ng = 200,   # genes
                   nc = 50,    # cells per fine‑type
                   nt = 3,     # # fine‑types ("type1", …)
                   ns = 4,     # # samples
                   nb = 2) {   # # batches
  stopifnot(requireNamespace("Seurat",  quietly = TRUE),
            requireNamespace("SeuratObject",  quietly = TRUE),
            requireNamespace("Matrix",  quietly = TRUE))
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    skip("SeuratObject not available for mock object")
  }
  ## 1) create global IDs --------------------------------------------------
  type_levels <- paste0("type",  seq_len(nt))
  maps        <- makeTypeHierarchy(type_levels)

  cell_ids <- unlist(lapply(seq_len(nt), function(t)
    paste0("cell", seq_len(nc), "_t", t)))         # unique over all types
  gene_ids <- paste0("gene", seq_len(ng))

  ## 2) counts matrix ------------------------------------------------------
  counts_vec <- stats::rpois(length(gene_ids) * length(cell_ids), lambda = 10)
  counts     <- Matrix::Matrix(
                 matrix(counts_vec, nrow = ng, dimnames = list(gene_ids, cell_ids)),
                 sparse = TRUE)

  ## 3) per‑cell metadata --------------------------------------------------
  fine_type <- rep(type_levels, each = nc)               # length == length(cell_ids)
  mid_type  <- maps$mid  [fine_type]
  broad_type<- maps$broad[fine_type]

  meta <- data.frame(
    type       = fine_type,
    fine_type  = fine_type,
    mid_type   = mid_type,
    broad_type = broad_type,
    batch      = sample(paste0("Batch",   seq_len(nb)), length(cell_ids), TRUE),
    sample     = sample(paste0("Sample",  seq_len(ns)), length(cell_ids), TRUE),
    row.names  = cell_ids,
    check.names = FALSE
  )

  ## 4) build Seurat object (single assay, no extra layers) ---------------
  se <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)

  # light preprocessing (optional)
  se <- Seurat::NormalizeData(se, verbose = FALSE)
  se <- Seurat::FindVariableFeatures(se, verbose = FALSE)
  se <- Seurat::ScaleData(se, verbose = FALSE)
  se <- Seurat::RunPCA(se, npcs = 5, verbose = FALSE)

  se@misc$pvclust <- list()  # placeholder slot for downstream tests
  se
}

mockLong <- function(nc = 50, nt = 3, ns = 4, nb = 2, useBatch = TRUE) {
  type_levels <- paste0("type", seq_len(nt))
  maps        <- makeTypeHierarchy(type_levels)

  df <- data.frame(
    bc       = paste0("cell", seq_len(nc)),
    type     = sample(type_levels, nc, TRUE),
    sample   = sample(paste0("sample", seq_len(ns)), nc, TRUE),
    stringsAsFactors = FALSE
  )
  if (useBatch)
    df$batch <- sample(paste0("batch", seq_len(nb)), nc, TRUE)

  df$fine_type  <- df$type
  df$mid_type   <- maps$mid  [df$type]
  df$broad_type <- maps$broad[df$type]
  df
}


mockCount <- function(df = mockLong()) {
  groupVars <- c("type", "sample")
  if ("batch" %in% names(df)) groupVars <- c(groupVars, "batch")
  aggregate(bc ~ ., data = df[, c("bc", groupVars)], FUN = length)
}

mockSCE <- function(nc = 500, nt = 3, ns = 4, nb = 2, useBatch = TRUE) {
  stopifnot(requireNamespace("SingleCellExperiment", quietly = TRUE))
  df  <- mockLong(nc, nt, ns, nb, useBatch)
  mat <- matrix(stats::rpois(nc * 20, lambda = 5), 20,
                dimnames = list(paste0("gene", seq_len(20)), df$bc))
  SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = mat),
    colData = df
  )
}

makeTypeHierarchy <- function(type_levels) {
    n  <- length(type_levels)
    i  <- seq_along(type_levels)
    list(
        mid   = setNames(paste0("mid",   ceiling(i / 2)),           type_levels),
        broad = setNames(paste0("broad", ifelse(i <= n / 2, 1, 2)), type_levels)
    )
}