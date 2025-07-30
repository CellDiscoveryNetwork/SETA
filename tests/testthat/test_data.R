# Generate mockLong data
set.seed(123)
x <- mockLong()
xc <- mockCount(x)

# Test data
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
test_that("makeTypeHierarchy returns correct structure", {
  types <- c("A", "B", "C", "D")
  result <- makeTypeHierarchy(types)

  expect_type(result, "list")
  expect_named(result, c("mid", "broad"))
  expect_length(result$mid, length(types))
  expect_length(result$broad, length(types))
  expect_named(result$mid, types)
  expect_named(result$broad, types)
})

test_that("makeTypeHierarchy assigns correct mid and broad levels", {
  types <- c("X", "Y", "Z", "W")
  result <- makeTypeHierarchy(types)

  expect_equal(result$mid,   c(X = "mid1", Y = "mid1", Z = "mid2", W = "mid2"))
  expect_equal(result$broad, c(X = "broad1", Y = "broad1", Z = "broad2", W = "broad2"))
})

# Test mock single cell data
set.seed(687)
se <- mockSC() # create Seurat object
sc <- mockSCE() # Create SingleCellExperiment object

test_that("mockSC returns a Seurat object", {
    skip_if_not_installed("Seurat")
    expect_true("Seurat" %in% class(se))
    expect_identical(ncol(se), 150)
    expect_identical(nrow(se), 200)
    expect_identical(
        colnames(se@meta.data),
        c("orig.ident", "nCount_RNA", "nFeature_RNA", "type",
          "fine_type", "mid_type", "broad_type", "batch", "sample"))
})


test_that("mockSCE returns a SingleCellExperiment object", {
    skip_if_not_installed("SingleCellExperiment")
    expect_true("SingleCellExperiment" %in% class(sc))
    expect_true(ncol(SingleCellExperiment::colData(sc)) > 0)
    expect_equal(ncol(sc), 500)
    expect_equal(nrow(sc), 20)
    expect_identical(
        colnames(SingleCellExperiment::colData(sc)),
        c("bc", "type", "sample", "batch",
          "fine_type", "mid_type", "broad_type"))
})