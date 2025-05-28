# Generate mockLong data
set.seed(123)
x <- mockLong()
xc <- mockCount(x)
test_that("mockLong returns correct structure", {
    expect_true(is.data.frame(x))
    expect_true(all(c("bc", "type", "sample", "batch") %in% colnames(x)))
    expect_equal(dim(x), c(50, 7))
    expect_true(all(x$type %in% paste0("type", 1:3)))
    expect_true(all(x$sample %in% paste0("sample", 1:4)))
    expect_true(all(x$batch %in% paste0("batch", 1:2)))
})

test_that("mockCount returns correct structure", {
    expect_true(is.data.frame(xc))
    expect_true(all(c("bc", "type", "sample", "batch") %in% colnames(xc)))
    expect_equal(dim(xc), c(21, 4))
    expect_true(all(xc$type %in% paste0("type", 1:3)))
    expect_true(all(xc$sample %in% paste0("sample", 1:4)))
    expect_true(all(xc$batch %in% paste0("batch", 1:2)))
})

test_that("mockCount returns correct data", {
    # Testing the 1st row
    expect_identical(xc[1, 1], "type1")
    expect_identical(xc[1, 2], "sample1")
    expect_identical(xc[1, 3], "batch1")
    expect_identical(xc[1, 4], as.integer(1))
    # testing a 12th row
    expect_identical(xc[12, 1], "type2")
    expect_identical(xc[12, 2], "sample1")
    expect_identical(xc[12, 3], "batch2")
    expect_identical(xc[12, 4], as.integer(2))
    
})

# test
makeTypeHierarchy