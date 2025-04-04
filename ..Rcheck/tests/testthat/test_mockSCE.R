test_that("mockSCE returns a SingleCellExperiment object", {
    skip_if_not_installed("SingleCellExperiment")
    x <- mockSCE()
    expect_true("SingleCellExperiment" %in% class(x))
    expect_true(ncol(colData(x)) > 0)
})