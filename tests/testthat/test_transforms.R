set.seed(687)

test_that("setaCLR returns correct known results", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
  # Hand-calculated CLR for this 2x2:
  # Row1: log(1/1.414214) = -0.3465739, log(2/1.414214) = 0.346573
  # Row2: log(4/5.656854) = -0.3465735, log(8/5.656854) = 0.346574
  expected <- matrix(c(-0.3465739, 0.346573, -0.346574, 0.346574),
                     nrow = 2, byrow = TRUE)
  out <- setaCLR(mat, pseudocount = 0)
  expect_equal(out$counts, expected, tolerance = 1e-3)
  expect_equal(out$method, "CLR")
})

test_that("setaALR returns correct known results", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
  colnames(mat) <- c("A", "B")
  # ALR with "A" as ref => log(B/A)
  # Row1: log(2/1)=0.69315, Row2: log(8/4)=0.69315
  expected <- matrix(c(0.69315, 0.69315), nrow = 2, byrow = TRUE)
  colnames(expected) <- "B"
  out <- setaALR(mat, ref = "A", pseudocount = 0)
  expect_equal(out$counts, expected, tolerance = 1e-3)
  expect_match(out$method, "ALR_ref=A")
})

test_that("setaILR with Helmert basis and boxcox_p != 0 works", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
  out <- setaILR(mat, boxcox_p=0.5)
  expect_match(out$method, "ILR_Helmert (boxcox_p=0.5)", fixed = TRUE)
  expect_equal(dim(out$counts), c(2, 1))
})

test_that("setaPercent yields correct percentages", {
  # Example 2x2
  mat <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  # column sums = (4, 6), so percentages for row1 col1 => (1/4)*100=25
  # row1 col2 => (2/6)*100=33.333...
  # row2 col1 => (3/4)*100=75, row2 col2 => (4/6)*100=66.666...
  expected <- matrix(c(33.3333, 66.6667, 42.8571, 57.1429), nrow = 2, byrow = TRUE)
  out <- setaPercent(mat)
  expect_equal(out$counts, expected, tolerance = 1e-3)
  expect_equal(out$method, "percent")
})

test_that("setaLogCPM normalizes as expected with default size factors (n cells per sample)", {
  # 2x2 matrix
  mat <- matrix(c(10, 100, 10, 100), nrow = 2, byrow = TRUE)
  # colSums => c(110, 110), so cpm for row1 col1 => (10 / 110) * 1e6 = ~90909...
  # log2(90909 + pseudocount) => ~16.47 if pseudocount=1
  out <- setaLogCPM(mat)
  expect_equal(out$method, "logCPM")
  expect_equal(dim(out$counts), dim(mat))
  # We won't hand-calc exact log2 values but we can do a quick sanity check:
  expect_true(all(out$counts > 0))  # Because big counts => logs are positive
})


test_that("Transforms work on mock SCE, Seurat, and long data", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")
  sce <- mockSCE()
  seu <- mockSC()
  df  <- mockLong()
  matSCE <- setaCounts(as.data.frame(SummarizedExperiment::colData(sce)),
                       bc = "bc")
  matSeurat <- setaCounts(seu@meta.data, bc = "rownames")
  matDF <- setaCounts(df)
  outSCE <- setaCLR(matSCE)
  outSeurat <- setaCLR(matSeurat)
  outDF <- setaCLR(matDF)
  expect_equal(outSCE$method, "CLR")
  expect_equal(outSeurat$method, "CLR")
  expect_equal(outDF$method, "CLR")
})

test_that("setaTransform works with all methods", {
  mat <- matrix(c(1, 3, 2, 4), nrow = 2, byrow = TRUE)
  # CLR
  res_clr <- setaTransform(mat, method = "CLR", pseudocount = 0.5)
  expect_equal(res_clr$method, "CLR")

  # ALR
  colnames(mat) <- c("TypeA", "TypeB")
  res_alr <- setaTransform(mat, method = "ALR", ref = "TypeA", pseudocount = 1)
  expect_match(res_alr$method, "ALR_ref=TypeA")

  # ILR
  res_ilr <- setaTransform(mat, method = "ILR", pseudocount = 1)
  expect_match(res_ilr$method, "ILR_Helmert")

  # percent
  res_pct <- setaTransform(mat, method = "percent")
  expect_equal(res_pct$method, "percent")

  ## balance
  # Mock 2 samples × 4 taxa  (AT1 AT2 Fib1 Fib2)
  cnt   <- matrix(c(1,2,3,4, 5,6,7,8), nrow = 2, byrow = TRUE)
  colnames(cnt) <- c("AT1","AT2","Fib1","Fib2")

  meta  <- data.frame(
    broad_type = c("Epi","Epi","Stroma","Stroma"),
    row.names  = colnames(cnt)
  )

  res_bal <- setaTransform(
    cnt,
    method            = "balance",
    balances          = list(num = "Epi", denom = "Stroma"),
    taxonomyDF        = meta,
    taxonomy_col      = "broad_type"
  )

  expect_equal(res_bal$method, "balance")
  expect_equal(ncol(res_bal$counts), 1)
})

test_that("setaBalance returns expected log‑ratio", {

  ## ------------------------------------------------------------------ ##
  ## 1.  Toy metadata  (matches help page) -----------------------------
  meta <- data.frame(
    bc          = paste0("cell", 1:6),
    fine_type   = c("AT1","AT2","AT1","Fib1","Fib1","AT2"),
    mid_type    = c("Alv","Alv","Alv","Fib","Fib","Alv"),
    broad_type  = c("Epi","Epi","Epi","Stroma","Stroma","Epi")
  )

  taxDF <- setaTaxonomyDF(meta,
            resolution_cols = c("broad_type","mid_type","fine_type"))

  ## leaves = AT1, AT2, Fib1  (3 leaves) -------------------------------
  leaves <- rownames(taxDF)

  ## ------------------------------------------------------------------ ##
  ## 2. Fake counts  (2 samples × 3 leaves) ----------------------------
  set.seed(123)
  cnt <- matrix(rpois(2 * length(leaves), 10), nrow = 2,
                dimnames = list(paste0("S", 1:2), leaves))


  ## ------------------------------------------------------------------ ##
  ## 3. Expected manual calculation ------------------------------------
  gm <- function(x) exp(mean(log(x + 1)))   # pseudocount = 1
  expected <- vapply(1:nrow(cnt), function(i) {
    log(gm(cnt[i, c("AT1","AT2")]) /
        gm(cnt[i, "Fib1"]))
  }, numeric(1))

  ## ------------------------------------------------------------------ ##
  ## 4. setaBalance call -----------------------------------------------
  out <- setaBalance(
           cnt,
           balances = list(num = "Epi", denom = "Stroma"),
           taxonomyDF   = taxDF,
           taxonomy_col = "broad_type",
           pseudocount  = 1
         )

  ## run setaBalance ------------------------------------------------------
  out <- setaBalance(cnt,
           balances   = list(num = "Epi",
                             denom = "Stroma"),
           taxonomyDF = taxDF,
           taxonomy_col = "broad_type",
           pseudocount = 1)

  expect_equal(out$method, "balance")
  expect_equal(as.numeric(out$counts[,1]), expected, tolerance = 1e-12)
})