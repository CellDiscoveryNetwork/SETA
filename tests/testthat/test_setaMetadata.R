# Generate mockLong data
set.seed(687)
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

test_that("setaMetadata removes sample_col from meta_cols automatically", {
  expect_true(all(c("sample_id", "batch") %in% colnames(df2)))
  expect_false("sample_batch" %in% colnames(df2))
})

test_that("setaMetadata throws if any meta_cols not in data", {
  df <- mockLong()
  expect_error(
    setaMetadata(df, sample_col = "sample", meta_cols = c("type", "nonexistent")),
    "not in your data"
  )
})

test_that("setaMetadata handles different column types correctly", {
  df <- mockLong()
  df$bool_col <- df$sample %in% unique(df$sample)[1:2]  # logical
  df$int_col <- as.integer(as.numeric(as.factor(df$sample)))  # integer
  df$num_col <- as.numeric(as.factor(df$sample))  # numeric
  df$char_col <- as.character(df$sample)  # character
  df$fact_col <- factor(df$sample, levels = unique(df$sample))  # factor

  df_unique <- df[!duplicated(df$sample), ]

  meta <- setaMetadata(
    x = df_unique,
    sample_col = "sample",
    meta_cols = c("bool_col", "int_col", "num_col", "char_col", "fact_col")
  )

  expect_type(meta$bool_col, "logical")
  expect_type(meta$int_col, "integer")
  expect_type(meta$num_col, "double")
  expect_type(meta$char_col, "character")
  expect_s3_class(meta$fact_col, "factor")
  expect_equal(levels(meta$fact_col), levels(df$fact_col))
})
