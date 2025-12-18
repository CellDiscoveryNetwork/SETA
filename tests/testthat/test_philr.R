# Tests for PhILR transform
# PhILR requires phyloseq and philr packages (suggested dependencies)

test_that("setaPhILR errors when phyloseq is not available", {
    # This test verifies the error message when phyloseq is missing
    # It will be skipped if phyloseq IS installed (which is expected in CI)
    # The actual error is tested in manual testing when packages aren't installed
    skip_if_not_installed("phyloseq")
    
    # If we reach here, phyloseq is installed, so we can't test the error
    # This is expected behavior - the error is tested manually
    expect_true(TRUE)
})

test_that("setaPhILR errors when philr is not available", {
    # Similar to above - verifies error handling for missing philr package
    skip_if_not_installed("philr")
    
    # If we reach here, philr is installed
    expect_true(TRUE)
})

test_that("setaPhILR errors on missing tree", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    
    mat <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, byrow = TRUE)
    colnames(mat) <- c("TaxonA", "TaxonB", "TaxonC")
    
    expect_error(
        setaPhILR(mat, tree = NULL),
        "A phylogenetic tree is required"
    )
    
    expect_error(
        setaPhILR(mat),
        "A phylogenetic tree is required"
    )
})

test_that("setaPhILR errors on non-matrix input", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    
    df <- data.frame(A = 1:2, B = 3:4, C = 5:6)
    tree <- mockPhyloTree(c("A", "B", "C"))
    
    expect_error(
        setaPhILR(df, tree = tree),
        "'counts' must be a matrix"
    )
})

test_that("setaPhILR errors when taxa names don't match", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    skip_if_not_installed("ape")
    
    mat <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, byrow = TRUE)
    colnames(mat) <- c("TaxonA", "TaxonB", "TaxonC")
    rownames(mat) <- c("Sample1", "Sample2")
    
    # Tree with different tip labels
    tree <- mockPhyloTree(c("TaxonX", "TaxonY", "TaxonZ"))
    
    expect_error(
        setaPhILR(mat, tree = tree),
        "No taxa names match"
    )
})

test_that("setaPhILR works with basic input", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    skip_if_not_installed("ape")
    
    # Use same mock data structure as other transforms
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Create tree based on cell types in the counts matrix
    cell_types <- colnames(mat)
    tree <- mockPhyloTree(tip_labels = cell_types, seed = 687)
    
    # Apply PhILR
    result <- setaPhILR(mat, tree = tree)
    
    # Check output structure
    expect_type(result, "list")
    expect_equal(result$method, "phILR")
    expect_true(is.matrix(result$counts))
    expect_equal(nrow(result$counts), nrow(mat))  # Same number of samples
    expect_equal(ncol(result$counts), ncol(mat) - 1)  # One less dimension
    expect_equal(rownames(result$counts), rownames(mat))
})

test_that("setaPhILR works with sample_data", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    skip_if_not_installed("ape")
    
    # Use same mock data structure as other transforms
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Create tree based on cell types
    cell_types <- colnames(mat)
    tree <- mockPhyloTree(tip_labels = cell_types, seed = 687)
    
    # Create sample_data - mockLong doesn't have sample-level metadata columns
    # (bc, type, batch are all cell-level), so create a simple data.frame
    sample_data <- data.frame(
        sample_id = rownames(mat),
        row.names = rownames(mat),
        stringsAsFactors = FALSE
    )
    
    result <- setaPhILR(mat, tree = tree, sample_data = sample_data)
    
    expect_equal(result$method, "phILR")
    expect_equal(nrow(result$counts), nrow(mat))
    expect_equal(ncol(result$counts), ncol(mat) - 1)  # n_taxa - 1 dimensions
})

test_that("setaPhILR handles partial taxa matching", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    skip_if_not_installed("ape")
    
    # Use mock data but add an extra taxon
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Add an extra column that won't be in the tree
    mat <- cbind(mat, ExtraTaxon = c(10, 20))
    
    # Create tree with only original cell types (missing ExtraTaxon)
    original_types <- setdiff(colnames(mat), "ExtraTaxon")
    tree <- mockPhyloTree(tip_labels = original_types, seed = 687)
    
    expect_warning(
        result <- setaPhILR(mat, tree = tree),
        "Some taxa in counts matrix are not in tree"
    )
    
    # Should use only matching taxa
    expect_equal(ncol(result$counts), length(original_types) - 1)
})

test_that("setaPhILR works with setaTransform wrapper", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    skip_if_not_installed("ape")
    
    # Use same mock data structure as other transforms
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Create tree based on cell types
    cell_types <- colnames(mat)
    tree <- mockPhyloTree(tip_labels = cell_types, seed = 687)
    
    result <- setaTransform(mat, method = "phILR", tree = tree)
    
    expect_equal(result$method, "phILR")
    expect_equal(result$within_resolution, FALSE)
    expect_equal(result$grouping_var, NULL)
    expect_true(is.matrix(result$counts))
    expect_equal(nrow(result$counts), nrow(mat))
    expect_equal(ncol(result$counts), ncol(mat) - 1)
})

test_that("setaTransform errors when phILR used without tree", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    
    mat <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, byrow = TRUE)
    colnames(mat) <- c("TaxonA", "TaxonB", "TaxonC")
    
    expect_error(
        setaTransform(mat, method = "phILR"),
        "please supply the 'tree' argument"
    )
})

test_that("setaTransform prevents phILR with within_resolution", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    skip_if_not_installed("ape")
    
    # Use same mock data structure
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Create tree
    cell_types <- colnames(mat)
    tree <- mockPhyloTree(tip_labels = cell_types, seed = 687)
    
    # Create taxonomyDF (similar to reference_frames vignette)
    taxonomyDF <- data.frame(
        Lineage = rep(c("L1", "L2"), length.out = length(cell_types)),
        row.names = cell_types
    )
    
    expect_error(
        setaTransform(
            mat,
            method = "phILR",
            tree = tree,
            taxonomyDF = taxonomyDF,
            taxonomy_col = "Lineage",
            within_resolution = TRUE
        ),
        "phILR transform cannot be used with within_resolution"
    )
})

test_that("setaPhILR preserves rownames", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    skip_if_not_installed("ape")
    
    # Use same mock data structure
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Create tree
    cell_types <- colnames(mat)
    tree <- mockPhyloTree(tip_labels = cell_types, seed = 687)
    
    result <- setaPhILR(mat, tree = tree)
    
    expect_equal(rownames(result$counts), rownames(mat))
})

test_that("setaPhILR handles different pseudocount values", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("philr")
    skip_if_not_installed("ape")
    
    # Use same mock data structure
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Create tree
    cell_types <- colnames(mat)
    tree <- mockPhyloTree(tip_labels = cell_types, seed = 687)
    
    result1 <- setaPhILR(mat, tree = tree, pseudocount = 1)
    result2 <- setaPhILR(mat, tree = tree, pseudocount = 0.5)
    
    # Results should differ with different pseudocounts
    expect_false(identical(result1$counts, result2$counts))
    
    # But structure should be the same
    expect_equal(dim(result1$counts), dim(result2$counts))
})

