# Generate mockLong data
set.seed(123)
sce <- mockSCE()
df <- setaMetadata(
    x = data.frame(colData(sce)),
    sample_col = "sample",
    meta_cols = c("batch", "broad_type"))
# Test data
test_that("mockLong returns correct structure", {
    expect_true(is.data.frame(x))
    expect_true(all(c("bc", "type", "sample", "batch") %in% colnames(x)))
    expect_equal(dim(x), c(50, 7))
    expect_true(all(x$type %in% paste0("type", 1:3)))
    expect_true(all(x$sample %in% paste0("sample", 1:4)))
    expect_true(all(x$batch %in% paste0("batch", 1:2)))
})


