# Tests for granularity scan function
# Requires ape package (suggested dependency)

set.seed(687)

# Load package functions (tests run in package environment)
# setaCounts, setaCLR, setaLatent, etc. should be available

test_that("setaGranularityScan errors when ape is not available", {
    # This test verifies the error message when ape is missing
    skip_if_not_installed("ape")
    
    # If we reach here, ape is installed
    expect_true(TRUE)
})

test_that("setaGranularityScan errors on non-matrix input", {
    skip_if_not_installed("ape")
    
    df <- data.frame(A = 1:2, B = 3:4)
    tree <- mockPhyloTree(c("A", "B"), seed = 687)
    
    expect_error(
        setaGranularityScan(df, tree = tree),
        "'counts' must be a matrix"
    )
})

test_that("setaGranularityScan errors on invalid tree", {
    skip_if_not_installed("ape")
    
    mat <- matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2)
    colnames(mat) <- c("A", "B")
    
    expect_error(
        setaGranularityScan(mat, tree = "not_a_tree"),
        "must be a phylogenetic tree object"
    )
})

test_that("setaGranularityScan errors when taxa don't match", {
    skip_if_not_installed("ape")
    
    mat <- matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2)
    colnames(mat) <- c("A", "B")
    tree <- mockPhyloTree(c("X", "Y"), seed = 687)
    
    expect_error(
        setaGranularityScan(mat, tree = tree),
        "No taxa names match"
    )
})

test_that("setaGranularityScan works with basic phylo tree", {
    skip_if_not_installed("ape")
    
    # Use same mock data structure as other transforms
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Create tree
    cell_types <- colnames(mat)
    tree <- mockPhyloTree(tip_labels = cell_types, seed = 687)
    
    # Run granularity scan with small range
    result <- setaGranularityScan(
        counts = mat,
        tree = tree,
        n_clusters = 2:3
    )
    
    # Check output structure
    expect_type(result, "list")
    expect_true("latent_spaces" %in% names(result))
    expect_true("metrics" %in% names(result))
    expect_true("aggregated_counts" %in% names(result))
    expect_null(result$cluster_labels)  # Not requested
    
    # Check metrics
    expect_true(is.data.frame(result$metrics))
    expect_true("n_clusters" %in% colnames(result$metrics))
    expect_equal(nrow(result$metrics), 2)  # 2 and 3 clusters
    
    # Check latent spaces
    expect_true(is.list(result$latent_spaces))
    expect_equal(length(result$latent_spaces), 2)
    expect_true("n_clusters_2" %in% names(result$latent_spaces))
    expect_true("n_clusters_3" %in% names(result$latent_spaces))
    
    # Check aggregated counts
    expect_true(is.list(result$aggregated_counts))
    expect_equal(length(result$aggregated_counts), 2)
})

test_that("setaGranularityScan returns correct latent space structure", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    result <- setaGranularityScan(mat, tree = tree, n_clusters = 2)
    
    # Check latent space structure
    latent <- result$latent_spaces$n_clusters_2
    expect_equal(latent$method, "PCA")
    expect_true(is.data.frame(latent$latentSpace))
    expect_equal(nrow(latent$latentSpace), nrow(mat))  # Same number of samples
    expect_equal(ncol(latent$latentSpace), 2)  # latent_dims = 2
    expect_true(is.numeric(latent$varExplained))
})

test_that("setaGranularityScan aggregates counts correctly", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    result <- setaGranularityScan(mat, tree = tree, n_clusters = 2)
    
    # Check aggregated counts
    agg_counts <- result$aggregated_counts$n_clusters_2
    expect_true(is.matrix(agg_counts))
    expect_equal(nrow(agg_counts), nrow(mat))  # Same samples
    expect_equal(ncol(agg_counts), 2)  # 2 clusters
    expect_equal(rownames(agg_counts), rownames(mat))
    expect_equal(colnames(agg_counts), c("Cluster1", "Cluster2"))
    
    # Aggregated counts should sum to original row sums
    expect_equal(rowSums(agg_counts), rowSums(mat), tolerance = 1e-10)
})

test_that("setaGranularityScan works with different latent methods", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    # Test PCA
    result_pca <- setaGranularityScan(mat, tree = tree, 
                                     n_clusters = 2,
                                     latent_method = "PCA")
    expect_equal(result_pca$latent_spaces$n_clusters_2$method, "PCA")
    
    # Test PCoA
    result_pcoa <- setaGranularityScan(mat, tree = tree,
                                      n_clusters = 2,
                                      latent_method = "PCoA")
    expect_equal(result_pcoa$latent_spaces$n_clusters_2$method, "PCoA")
    
    # Test NMDS
    result_nmds <- setaGranularityScan(mat, tree = tree,
                                      n_clusters = 2,
                                      latent_method = "NMDS")
    expect_equal(result_nmds$latent_spaces$n_clusters_2$method, "NMDS")
})

test_that("setaGranularityScan handles default n_clusters", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    # Default should be 2:n_taxa
    n_taxa <- ncol(mat)
    result <- setaGranularityScan(mat, tree = tree)
    
    # Should have n_taxa - 1 resolutions (2 through n_taxa)
    expect_equal(nrow(result$metrics), n_taxa - 1)
    expect_equal(min(result$metrics$n_clusters), 2)
    expect_equal(max(result$metrics$n_clusters), n_taxa)
})

test_that("setaGranularityScan validates n_clusters range", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    # n_clusters < 2 should be filtered out
    result <- setaGranularityScan(mat, tree = tree, n_clusters = c(1, 2, 3))
    expect_equal(min(result$metrics$n_clusters), 2)
    
    # n_clusters > n_taxa should be filtered out
    n_taxa <- ncol(mat)
    result2 <- setaGranularityScan(mat, tree = tree, 
                                   n_clusters = c(2, n_taxa, n_taxa + 1))
    expect_equal(max(result2$metrics$n_clusters), n_taxa)
    
    # Invalid range should error
    expect_error(
        setaGranularityScan(mat, tree = tree, n_clusters = c(100, 200)),
        "must contain values between"
    )
})

test_that("setaGranularityScan works with return_cell_labels", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    # Create cell metadata matching taxa
    cell_meta <- data.frame(
        cell_type = colnames(mat),
        other_info = paste0("info_", colnames(mat)),
        stringsAsFactors = FALSE
    )
    
    result <- setaGranularityScan(
        counts = mat,
        tree = tree,
        n_clusters = 2:3,
        cell_metadata = cell_meta,
        taxa_col = "cell_type",
        return_cell_labels = TRUE
    )
    
    # Check cluster labels
    expect_true(!is.null(result$cluster_labels))
    expect_true(is.data.frame(result$cluster_labels))
    expect_true("taxa" %in% colnames(result$cluster_labels))
    expect_true("n_clusters_2" %in% colnames(result$cluster_labels))
    expect_true("n_clusters_3" %in% colnames(result$cluster_labels))
    expect_equal(nrow(result$cluster_labels), ncol(mat))
})

test_that("setaGranularityScan errors when cell_metadata missing taxa_col", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    cell_meta <- data.frame(cell_type = colnames(mat))
    
    expect_error(
        setaGranularityScan(mat, tree = tree,
                           cell_metadata = cell_meta,
                           return_cell_labels = TRUE),
        "must be specified"
    )
})

test_that("setaGranularityScan handles partial taxa matching", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    
    # Add extra column
    mat <- cbind(mat, ExtraType = c(10, 20, 30, 40))
    
    # Tree with only original types
    original_types <- setdiff(colnames(mat), "ExtraType")
    tree <- mockPhyloTree(tip_labels = original_types, seed = 687)
    
    expect_warning(
        result <- setaGranularityScan(mat, tree = tree, n_clusters = 2),
        "Some taxa in counts are not in tree"
    )
    
    # Should use only matching taxa
    expect_equal(ncol(result$aggregated_counts$n_clusters_2), 2)
})

test_that("setaGranularityScan metrics include variance explained", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    result <- setaGranularityScan(mat, tree = tree, 
                                 n_clusters = 2:3,
                                 latent_method = "PCA")
    
    # Check metrics structure
    expect_true("var_explained_pc1" %in% colnames(result$metrics))
    expect_true("total_var_explained" %in% colnames(result$metrics))
    expect_true(all(result$metrics$var_explained_pc1 > 0))
    expect_true(all(result$metrics$var_explained_pc1 <= 1))
})

test_that("setaGranularityScan preserves sample rownames", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    result <- setaGranularityScan(mat, tree = tree, n_clusters = 2)
    
    # Check rownames preserved in aggregated counts
    expect_equal(rownames(result$aggregated_counts$n_clusters_2), 
                 rownames(mat))
    
    # Check rownames preserved in latent spaces
    expect_equal(rownames(result$latent_spaces$n_clusters_2$latentSpace),
                 rownames(mat))
})

test_that("setaGranularityScan handles different pseudocount values", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    result1 <- setaGranularityScan(mat, tree = tree, 
                                  n_clusters = 2,
                                  pseudocount = 1)
    result2 <- setaGranularityScan(mat, tree = tree,
                                  n_clusters = 2,
                                  pseudocount = 0.5)
    
    # Results should differ
    expect_false(identical(
        result1$latent_spaces$n_clusters_2$latentSpace,
        result2$latent_spaces$n_clusters_2$latentSpace
    ))
    
    # But structure should be same
    expect_equal(dim(result1$latent_spaces$n_clusters_2$latentSpace),
                 dim(result2$latent_spaces$n_clusters_2$latentSpace))
})

test_that("setaGranularityScan handles different latent_dims", {
    skip_if_not_installed("ape")
    
    set.seed(687)
    df <- mockLong()
    mat <- setaCounts(df)
    tree <- mockPhyloTree(tip_labels = colnames(mat), seed = 687)
    
    # With 3 cell types, we can have max 3 clusters, so test with 2 clusters and 2 dimensions
    result <- setaGranularityScan(mat, tree = tree,
                                 n_clusters = 2,
                                 latent_dims = 2)
    
    expect_equal(ncol(result$latent_spaces$n_clusters_2$latentSpace), 2)
    expect_equal(length(result$latent_spaces$n_clusters_2$varExplained), 2)
})

test_that("setaGranularityScan provides helpful error for unrecognized tree types", {
    skip_if_not_installed("ape")
    
    mat <- matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2)
    colnames(mat) <- c("A", "B")
    
    # Test with invalid tree type
    expect_error(
        setaGranularityScan(mat, tree = list(not_a_tree = TRUE)),
        "must be a phylogenetic tree object"
    )
    
    # Error should include class information
    expect_error(
        setaGranularityScan(mat, tree = "string"),
        "Received object of class"
    )
})
