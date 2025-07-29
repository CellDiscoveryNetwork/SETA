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
    expect_equal(ncol(sc), 500)
    expect_equal(nrow(sc), 20)
    expect_identical(
        colnames(SingleCellExperiment::colData(sc)),
        c("bc", "type", "sample", "batch",
          "fine_type", "mid_type", "broad_type"))
})
