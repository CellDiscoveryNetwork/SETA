set.seed(687)

test_that("setaCounts returns a correctly‑shaped matrix from metadata frames", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")

  sce <- mockSCE()
  sce_meta <- as.data.frame(SingleCellExperiment::colData(sce))
  sce_meta$bc <- rownames(sce_meta)
  seu_meta <- mockSC()@meta.data
  seu_meta$bc <- rownames(seu_meta)
  long_df  <- mockLong()

  expect_dims <- function(df, sample_col = "sample", type_col = "type") {
    c(length(unique(df[[sample_col]])),
      length(unique(df[[type_col]])))
  }

  for (meta in list(sce_meta, seu_meta, long_df)) {
    mat <- setaCounts(meta)
    expect_true(is.matrix(mat))
    expect_equal(dim(mat), expect_dims(meta))
  }

  df_rn            <- long_df
  rownames(df_rn)  <- df_rn$bc
  df_rn$bc         <- NULL
  mat_rn           <- setaCounts(df_rn, bc_col = "rownames")
  mat_regular      <- setaCounts(long_df)
  expect_identical(mat_rn, mat_regular)
})

test_that("setaCounts handles factor columns", {
  df <- data.frame(
    bc = paste0("cell", 1:6),
    type = factor(c("A", "A", "B", "B", "C", "C")),
    sample = factor(c("S1", "S1", "S2", "S2", "S3", "S3"))
  )
  expect_message(mat <- setaCounts(df), "Converting factor class columns to character")
  expect_equal(dim(mat), c(3, 3))  # 3 samples x 3 types
})

test_that("setaCounts warns on invalid sample IDs", {
  df <- data.frame(
    bc = paste0("cell", 1:6),
    type = rep("A", 6),
    sample = c("S(1)", "S*2", "S/3", "S@4", "S-4", "S_5")  # "*", "@", "()", "/" are invalid
  )
  expect_warning(setaCounts(df), "special characters")
})

test_that("setaTaxonomyDF builds one‑to‑one taxonomy frames", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")

  res_cols <- c("broad_type", "mid_type", "fine_type")

  sce <- mockSCE()
  sce_meta <- as.data.frame(SingleCellExperiment::colData(sce))
  sce_meta$bc <- rownames(sce_meta)
  seu_meta <- mockSC()@meta.data
  seu_meta$bc <- rownames(seu_meta)
  long_df  <- mockLong()

  check_tax <- function(meta, bc = "bc") {
    tax <- setaTaxonomyDF(meta, resolution_cols = res_cols, bc_col = bc)
    expect_true(is.data.frame(tax))
    expect_equal(ncol(tax), length(res_cols))
    expect_equal(colnames(tax), res_cols)
    expect_equal(nrow(tax), length(unique(meta[[tail(res_cols, 1)]])))
    expect_identical(rownames(tax), tax[[tail(res_cols, 1)]])
  }

  lapply(list(sce_meta, seu_meta, long_df), check_tax)

  df_rn            <- long_df
  rownames(df_rn)  <- df_rn$bc
  df_rn$bc         <- NULL
  check_tax(df_rn, bc = "rownames")
})

test_that("setaTaxonomyDF errors on non-1:1 finest mapping", {
  df <- data.frame(
    bc = paste0("cell", 1:4),
    fine_type = c("FT1", "FT1", "FT2", "FT2"),
    mid_type = c("MT1", "MT2", "MT2", "MT2"),
    broad_type = c("BT1", "BT2", "BT2", "BT2")
  )
  expect_error(
    setaTaxonomyDF(df, resolution_cols = c("broad_type", "mid_type", "fine_type")),
    "Finest labels mapping to >1 coarser combo"
  )
})

test_that("setaTaxonomyDF errors on NA in finest label", {
  df <- data.frame(
    bc = paste0("cell", 1:3),
    fine_type = c("FT1", NA, "FT3"),
    mid_type = c("MT1", "MT2", "MT3"),
    broad_type = c("BT1", "BT2", "BT3")
  )
  expect_error(
    setaTaxonomyDF(df, resolution_cols = c("broad_type", "mid_type", "fine_type")),
    "NA in the finest label"
  )
})

test_that("resolveGroup resolves character leaf names correctly", {
  counts <- matrix(1, nrow = 3, ncol = 4)
  colnames(counts) <- c("AT1", "AT2", "Fib1", "Fib2")

  result <- resolveGroup(c("AT1", "Fib2"), counts)
  expect_equal(result, c(1, 4))
})

test_that("resolveGroup resolves higher-level labels via taxonomyDF", {
  counts <- matrix(1, nrow = 3, ncol = 4,
                   dimnames = list(NULL, c("AT1", "AT2", "Fib1", "Fib2")))
  taxDF <- data.frame(broad_type = c("Epi", "Epi", "Stroma", "Stroma"),
                      row.names  = colnames(counts))

  result <- resolveGroup("Stroma", counts, taxDF, taxonomy_col = "broad_type")
  expect_equal(result, c(3, 4))  # Fib1 and Fib2
})

test_that("resolveGroup handles numeric indices with bounds checking", {
  counts <- matrix(1, nrow = 3, ncol = 4)
  result <- resolveGroup(c(2, 4), counts)
  expect_equal(result, c(2, 4))

  expect_error(resolveGroup(5, counts), "Numeric group spec out of range")
})

test_that("resolveGroup errors on invalid input type", {
  counts <- matrix(1, nrow = 2, ncol = 2)
  expect_error(
    resolveGroup(list("bad_input"), counts),
    "Group spec must be character or numeric"
  )
})

test_that("resolveGroup errors when matching labels not in counts colnames", {
  counts <- matrix(1, nrow = 2, ncol = 2,
                   dimnames = list(NULL, c("G1", "G2")))
  taxonomyDF <- data.frame(group = c("A", "B"),
                           row.names = c("X", "Y"))  # X, Y not in counts
  expect_error(
    resolveGroup("A", counts, taxonomyDF, taxonomy_col = "group"),
    "Unresolved taxa/labels"
  )
})

test_that("taxonomy_to_tbl_graph builds tree with metadata", {
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")

  df <- data.frame(
    broad = c("Epi", "Epi", "Stroma"),
    mid   = c("Alv", "Alv", "Fib"),
    fine  = c("AT1", "AT2", "Fib1")
  )

  tg <- taxonomy_to_tbl_graph(df, columns = c("broad", "mid", "fine"), root_name = "Root")

  expect_s3_class(tg, "tbl_graph")
  nodes <- as.data.frame(tg, "nodes")

  expect_true("name" %in% names(nodes))
  expect_true("broad" %in% names(nodes))
  expect_equal(nodes$name[1], "Root")
})

test_that("taxonomy_to_tbl_graph errors with no columns", {
  df <- data.frame()
  expect_error(taxonomy_to_tbl_graph(df, columns = character(0)),
               "Need at least one column")
})