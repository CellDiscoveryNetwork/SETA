test_that("setaCounts extracts matrix from multiple object types", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")

  sce <- mockSCE()
  seu <- mockSC()
  df  <- mockLong()
  mtx <- matrix(1:4, nrow = 2)

  mat_sce     <- setaCounts(sce)
  mat_seurat  <- setaCounts(seu)
  mat_df      <- setaCounts(df)
  mat_matrix  <- setaCounts(mtx)

  expect_true(is.matrix(mat_sce))
  expect_true(is.matrix(mat_seurat))
  expect_true(is.matrix(mat_df))
  expect_true(is.matrix(mat_matrix))
  # Additional checks can be made if we know expected dimensions
})