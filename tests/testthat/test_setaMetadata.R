# Generate mockLong data
set.seed(123)
sce <- mockSCE()

test_that("setaMetadata errors on multiplexed samples", {
  df <- as.data.frame(SingleCellExperiment::colData(sce))
  expect_error(
    setaMetadata(df,
                 sample_col = "sample",
                 meta_cols = c("batch")),
    "Are your samples multiplexed"
  )
})

meta_df <- as.data.frame(SingleCellExperiment::colData(sce))
meta_df$sample_batch <- paste(meta_df$sample, meta_df$batch, sep = "_")

df2 <- setaMetadata(
  x          = meta_df,
  sample_col = "sample_batch",
  meta_cols  = c("batch")
)

# Test data
test_that("setaMetadata returns correct structure", {
  expect_true(is.data.frame(df2))
  expect_true(all(c("sample_id", "batch") %in% colnames(df2)))
  expect_equal(dim(df2), c(8, 2))
  expect_true(all(df2$batch  %in% paste0("batch", 1:2)))
})