# Test mock single cell data
set.seed(123)
se <- mockSC() # create Seurat object
sc <- mockSCE() # Create SingleCellExperiment object

test_that("mockSC returns a Seurat object", {
    expect_true("Seurat" %in% class(se))
    expect_identical(ncol(se), 150)
    expect_identical(nrow(se), 200)
    expect_identical(
        colnames(se@meta.data),
        c("orig.ident", "nCount_RNA", "nFeature_RNA", "type",
          "fine_type", "mid_type", "broad_type", "batch", "sample"))
    expect_identical(names(se@misc), "pvclust")
    expect_type(se@misc$pvclust, "list")
})


test_that("mockSCE returns a SingleCellExperiment object", {
    expect_true("SingleCellExperiment" %in% class(sc))
    expect_equal(ncol(sc), 150)
    expect_equal(nrow(sc), 200)
    expect_identical(
        colnames(colData(sc)),
        c("orig.ident", "nCount_RNA", "nFeature_RNA", "type",
          "fine_type", "mid_type", "broad_type", "batch", "sample"))
})
