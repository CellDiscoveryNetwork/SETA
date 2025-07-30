test_that("setaDistances works with valid CLR-transformed matrix", {
    # Simulated CLR output
    clr_mat <- list(
        counts = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
    )
    rownames(clr_mat$counts) <- c("SampleA", "SampleB")
    colnames(clr_mat$counts) <- c("Taxon1", "Taxon2", "Taxon3")
    result <- setaDistances(clr_mat)
    expect_s3_class(result, "data.frame")
    expect_named(result, c("from", "to", "distance"))
    expect_equal(nrow(result), 1)  # only 1 pair for 2 samples
    expect_true(result$from == "SampleA" && result$to == "SampleB")
    expect_equal(result$distance, sqrt(sum((clr_mat$counts[1, ] - clr_mat$counts[2, ])^2)))
})

test_that("setaDistances removes self and symmetric duplicates", {
    mat <- matrix(runif(12), nrow = 4)
    rownames(mat) <- paste0("Sample", 1:4)
    clr_mat <- list(counts = mat)
    result <- setaDistances(clr_mat)
    # 4 samples = choose(4, 2) = 6 pairs
    expect_equal(nrow(result), 6)
    expect_false(any(result$from == result$to))
    expect_true(all(as.character(result$from) < as.character(result$to)))
})

test_that("setaDistances respects distance method", {
    mat <- matrix(1:9, nrow = 3)
    rownames(mat) <- c("A", "B", "C")
    clr_mat <- list(counts = mat)

    result1 <- setaDistances(clr_mat, method = "euclidean")
    result2 <- setaDistances(clr_mat, method = "manhattan")

    expect_true(all(result1$distance != result2$distance))
})

test_that("setaDistances fails gracefully on bad input", {
    bad_input <- list(counts = as.data.frame(matrix(1:4, nrow = 2)))

    expect_error(setaDistances(bad_input), 
                 "'transformed_counts' must be a numeric matrix")
})
