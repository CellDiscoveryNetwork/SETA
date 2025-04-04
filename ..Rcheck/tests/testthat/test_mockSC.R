test_that("mockSC returns a Seurat object", {
    skip_if_not_installed("Seurat")
    x <- mockSC()
    expect_true("Seurat" %in% class(x))
    expect_true(ncol(x@meta.data) > 0)
})