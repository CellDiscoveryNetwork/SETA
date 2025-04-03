test_that("setaCLR returns correct known results", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
  # Hand-calculated CLR for this 2x2:
  # Row1: log(1/1.4142) = -0.3466, log(2/1.4142) = 0.3466
  # Row2: log(4/5.6569) = -0.3483, log(8/5.6569) = 0.3448
  expected <- matrix(c(-0.3466, 0.3466, -0.3483, 0.3448),
                     nrow = 2, byrow = TRUE)
  out <- setaCLR(mat)
  expect_equal(out$counts, expected, tolerance = 1e-3)
  expect_equal(out$method, "CLR")
})

test_that("setaALR returns correct known results", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
  colnames(mat) <- c("A", "B")
  # ALR with "A" as ref => log(B/A)
  # Row1: log(2/1)=0.6931, Row2: log(8/4)=0.6931
  expected <- matrix(c(0.6931, 0.6931), nrow = 2, byrow = TRUE)
  out <- setaALR(mat, ref = "A")
  expect_equal(out$counts, expected, tolerance = 1e-3)
  expect_match(out$method, "ALR_ref=A")
})

test_that("setaILR with Helmert basis works on small matrix", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
  out <- setaILR(mat)
  expect_equal(dim(out$counts), c(2, 1))
  expect_equal(out$method, "ILR_Helmert")
})

test_that("Transforms work on mock SCE, Seurat, and long data", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")
  sce <- mockSCE()
  seu <- mockSC()
  df  <- mockLong()
  matSCE <- setaCounts(sce)
  matSeurat <- setaCounts(seu)
  matDF <- setaCounts(df)
  outSCE <- setaCLR(matSCE)
  outSeurat <- setaCLR(matSeurat)
  outDF <- setaCLR(matDF)
  expect_equal(outSCE$method, "CLR")
  expect_equal(outSeurat$method, "CLR")
  expect_equal(outDF$method, "CLR")
})