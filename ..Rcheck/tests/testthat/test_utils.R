test_that("setaCounts extracts matrix from multiple object types", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")

  sce <- mockSCE()
  seu <- mockSC()
  df  <- mockLong()

  mat_sce     <- setaCounts(sce)
  mat_seurat  <- setaCounts(seu)
  mat_df      <- setaCounts(df)

  expect_true(is.matrix(mat_sce))
  expect_true(is.matrix(mat_seurat))
  expect_true(is.matrix(mat_df))
  # Additional checks can be made if we know expected dimensions
})