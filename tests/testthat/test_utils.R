test_that("setaCounts returns a correctly‑shaped matrix from metadata frames", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")

  ## build three metadata frames -------------------------------------------
  sce <- mockSCE()
  sce_meta <- as.data.frame(SingleCellExperiment::colData(sce))
  sce_meta$bc <- rownames(sce_meta)
  seu_meta <- mockSC()@meta.data
  seu_meta$bc <- rownames(seu_meta)
  long_df  <- mockLong()

  ## helper to compute expected dims ---------------------------------------
  expect_dims <- function(df, sample_col = "sample", type_col = "type") {
    c(length(unique(df[[sample_col]])),
      length(unique(df[[type_col]])))
  }

  for (meta in list(sce_meta, seu_meta, long_df)) {
    mat <- setaCounts(meta)
    expect_true(is.matrix(mat))
    expect_equal(dim(mat), expect_dims(meta))
  }

  ## rownames barcode shortcut --------------------------------------------
  df_rn            <- long_df
  rownames(df_rn)  <- df_rn$bc
  df_rn$bc         <- NULL
  mat_rn           <- setaCounts(df_rn, bc_col = "rownames")
  mat_regular      <- setaCounts(long_df)
  expect_identical(mat_rn, mat_regular)
})

test_that("setaTaxonomyDF builds one‑to‑one taxonomy frames", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")

  ## resolution hierarchy (coarse to fine) ----------------------------------
  res_cols <- c("broad_type", "mid_type", "fine_type")

  sce <- mockSCE()
  sce_meta <- as.data.frame(SingleCellExperiment::colData(sce))
  sce_meta$bc <- rownames(sce_meta)
  seu_meta <- mockSC()@meta.data
  seu_meta$bc <- rownames(seu_meta)
  long_df  <- mockLong()

  ## helper to test a metadata frame ---------------------------------------
  check_tax <- function(meta, bc = "bc") {
    tax <- setaTaxonomyDF(meta, resolution_cols = res_cols, bc_col = bc)
    expect_true(is.data.frame(tax))
    expect_equal(ncol(tax), length(res_cols))
    expect_equal(colnames(tax), res_cols)
    expect_equal(nrow(tax), length(unique(meta[[tail(res_cols, 1)]])))
    expect_identical(rownames(tax), tax[[tail(res_cols, 1)]])
  }

  lapply(list(sce_meta, seu_meta, long_df), check_tax)

  ## rownames barcode shortcut --------------------------------------------
  df_rn            <- long_df
  rownames(df_rn)  <- df_rn$bc
  df_rn$bc         <- NULL
  check_tax(df_rn, bc = "rownames")
})