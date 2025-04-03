test_that("mockLong returns correct structure", {
    x <- mockLong()
    expect_true(is.data.frame(x))
    expect_true(all(c("bc", "type", "sample", "batch") %in% colnames(x)))
})
