#' Granularity Scan Across Tree Resolutions
#'
#' Performs a granularity scan by cutting a phylogenetic tree at different
#' numbers of clusters and computing latent spaces for each resolution.
#' This allows exploration of how compositional structure changes across
#' different levels of taxonomic aggregation.
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#'   Column names (taxa names) must match the tip labels of the phylogenetic tree.
#' @param tree A phylogenetic tree object. Can be:
#'   \itemize{
#'     \item class \code{phylo} (from \code{ape} package)
#'     \item class \code{Node} (from \code{data.tree} package, e.g., from ARBOL)
#'   }
#'   The tip labels must match \code{colnames(counts)}.
#' @param n_clusters Integer vector. Range of cluster numbers to test.
#'   Default is \code{2:length(colnames(counts))} (all possible resolutions).
#' @param cell_metadata Optional. A data.frame with cell-level metadata.
#'   If provided, must contain a column matching taxa names (for output a).
#'   Should have one row per cell/taxa.
#' @param taxa_col Character. Column name in \code{cell_metadata} that contains
#'   taxa names matching \code{colnames(counts)}. Only used if \code{cell_metadata}
#'   is provided.
#' @param latent_method Character. Method for latent space analysis.
#'   One of \code{"PCA"}, \code{"PCoA"}, or \code{"NMDS"}. Default is \code{"PCA"}.
#' @param latent_dims Integer. Number of dimensions for latent space.
#'   Default is 2.
#' @param pseudocount Numeric. Pseudocount for CLR transform. Default is 1.
#' @param return_cell_labels Logical. If \code{TRUE} and \code{cell_metadata}
#'   is provided, returns cluster labels per cell for each resolution.
#'   Default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{cluster_labels}{If \code{return_cell_labels = TRUE} and
#'     \code{cell_metadata} provided: A data.frame with cluster assignments
#'     for each cell/taxa at each resolution. Otherwise \code{NULL}.}
#'   \item{latent_spaces}{A list of \code{setaLatent} output objects, one
#'     for each resolution. Names are \code{"n_clusters_<N>"}.}
#'   \item{metrics}{A data.frame with metrics across resolutions:
#'     \itemize{
#'       \item \code{n_clusters}: Number of clusters at this resolution
#'       \item \code{var_explained_pc1}: Variance explained by PC1 (if PCA)
#'       \item \code{var_explained_pc2}: Variance explained by PC2 (if PCA)
#'       \item \code{total_var_explained}: Total variance explained (if PCA/PCoA)
#'       \item \code{stress}: Stress value (if NMDS)
#'     }}
#'   \item{aggregated_counts}{A list of aggregated count matrices, one per
#'     resolution. Names are \code{"n_clusters_<N>"}.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Converts the phylo tree to hclust format
#'   \item For each number of clusters in \code{n_clusters}:
#'     \itemize{
#'       \item Cuts the tree using \code{cutree()}
#'       \item Aggregates counts by cluster (sums within each cluster)
#'       \item Applies CLR transform to aggregated counts
#'       \item Computes latent space using specified method
#'       \item Stores results and metrics
#'     }
#' }
#'
#' This allows users to explore how compositional structure changes as taxa
#' are aggregated at different resolutions of the tree.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   # Create example data
#'   set.seed(687)
#'   mat <- matrix(rpois(20, lambda = 10), nrow = 4, ncol = 5)
#'   colnames(mat) <- paste0("Type", 1:5)
#'   rownames(mat) <- paste0("Sample", 1:4)
#'   
#'   # Create a tree
#'   tree <- ape::rtree(n = 5, tip.label = colnames(mat))
#'   
#'   # Run granularity scan
#'   scan_result <- setaGranularityScan(mat, tree, n_clusters = 2:4)
#'   
#'   # View metrics
#'   scan_result$metrics
#'   
#'   # Access latent space for a specific resolution
#'   scan_result$latent_spaces$n_clusters_3
#' }
#' }
#'
#' @importFrom stats cutree hclust dist
#' @export
setaGranularityScan <- function(counts,
                                tree,
                                n_clusters = NULL,
                                cell_metadata = NULL,
                                taxa_col = NULL,
                                latent_method = c("PCA", "PCoA", "NMDS"),
                                latent_dims = 2,
                                pseudocount = 1,
                                return_cell_labels = FALSE) {
    
    # Check for required packages
    if (!requireNamespace("ape", quietly = TRUE)) {
        stop("The 'ape' package is required for granularity scan.\n",
             "Please install it with: install.packages('ape')")
    }
    
    # Validate inputs
    if (!is.matrix(counts))
        stop("'counts' must be a matrix with samples in rows and taxa in columns.")
    
    # Handle different tree formats
    tree_phylo <- NULL
    
    # Check for phylo first
    if (inherits(tree, "phylo")) {
        tree_phylo <- tree
    } else {
        # Check for Node object (data.tree R6 class)
        # Node objects have "Node" in their class and are R6 objects
        is_node <- FALSE
        tree_classes <- class(tree)
        
        # Check for Node class (multiple ways to be safe)
        if ("Node" %in% tree_classes || 
            any(grepl("Node", tree_classes, fixed = TRUE)) ||
            (is.object(tree) && "R6" %in% tree_classes && 
             exists("name", envir = tree, inherits = FALSE))) {
            is_node <- TRUE
        }
        
        # Also check if conversion functions exist - if they do, likely a Node
        if (!is_node && (exists("as.phylo.NodeI", mode = "function") ||
                         exists("as.phylo.Node", mode = "function"))) {
            # Try conversion - if it works, treat as Node
            as_phylo_test <- NULL
            if (exists("as.phylo.NodeI", mode = "function")) {
                as_phylo_test <- get("as.phylo.NodeI", mode = "function")
            } else {
                as_phylo_test <- get("as.phylo.Node", mode = "function")
            }
            tryCatch({
                test_result <- as_phylo_test(tree)
                if (inherits(test_result, "phylo")) {
                    is_node <- TRUE
                }
            }, error = function(e) {
                # Not a Node or conversion failed
            })
        }
        
        if (is_node) {
            # Convert Node (data.tree) to phylo
            # ARBOL and other packages may provide conversion methods
            if (!requireNamespace("ape", quietly = TRUE)) {
                stop("The 'ape' package is required to convert Node objects to phylo.\n",
                     "Please install it with: install.packages('ape')")
            }
            
            # Try multiple conversion methods
            conversion_success <- FALSE
            
            # Method 1: Try as.phylo.Node or as.phylo.NodeI (ARBOL or similar packages)
            # Check in global environment and loaded packages
            as_phylo_node <- NULL
            
            # Try as.phylo.NodeI first (from user's old code)
            if (exists("as.phylo.NodeI", mode = "function")) {
                as_phylo_node <- get("as.phylo.NodeI", mode = "function")
            } else if (exists("as.phylo.Node", mode = "function")) {
                as_phylo_node <- get("as.phylo.Node", mode = "function")
            } else {
                # Try to find in loaded namespaces (e.g., ARBOL package)
                loaded_ns <- loadedNamespaces()
                for (ns in loaded_ns) {
                    if (exists("as.phylo.NodeI", where = paste0("package:", ns), mode = "function")) {
                        as_phylo_node <- get("as.phylo.NodeI", envir = asNamespace(ns))
                        break
                    } else if (exists("as.phylo.Node", where = paste0("package:", ns), mode = "function")) {
                        as_phylo_node <- get("as.phylo.Node", envir = asNamespace(ns))
                        break
                    }
                }
            }
            
            if (!is.null(as_phylo_node)) {
                tryCatch({
                    tree_phylo <- as_phylo_node(tree)
                    if (inherits(tree_phylo, "phylo")) {
                        conversion_success <- TRUE
                    }
                }, error = function(e) {
                    # Continue to next method
                })
            }
            
            # Method 2: Try ape::as.phylo if Node has standard structure
            if (!conversion_success) {
                tryCatch({
                    tree_phylo <- ape::as.phylo(tree)
                    if (inherits(tree_phylo, "phylo")) {
                        conversion_success <- TRUE
                    }
                }, error = function(e) {
                    # Continue to next method
                })
            }
            
            # Method 3: Try data.tree::as.phylo if available
            if (!conversion_success && requireNamespace("data.tree", quietly = TRUE)) {
                if (exists("as.phylo", where = "package:data.tree")) {
                    tryCatch({
                        tree_phylo <- data.tree::as.phylo(tree)
                        if (inherits(tree_phylo, "phylo")) {
                            conversion_success <- TRUE
                        }
                    }, error = function(e) {
                        # Continue to error message
                    })
                }
            }
            
            if (!conversion_success) {
                stop("Cannot automatically convert Node object to phylo.\n",
                     "Please convert your Node object to phylo first. Options:\n",
                     "  1. If using ARBOL: tree <- as.phylo.NodeI(your_node_object)\n",
                     "  2. If using data.tree: tree <- data.tree::as.phylo(your_node_object)\n",
                     "  3. Or provide a phylo object directly.\n",
                     "If you have a custom conversion function, please convert before calling setaGranularityScan().")
            }
            
            if (!inherits(tree_phylo, "phylo")) {
                stop("Tree conversion did not produce a phylo object. ",
                     "Please convert your Node object to phylo before calling this function.")
            }
        } else {
            # Provide helpful error message with class information
            tree_class <- paste(class(tree), collapse = ", ")
            stop("'tree' must be a phylogenetic tree object (class 'phylo') or Node object (class 'Node').\n",
                 "Received object of class: ", tree_class, "\n",
                 "If this is a Node object from ARBOL, please convert it first:\n",
                 "  tree <- as.phylo.NodeI(your_node_object)\n",
                 "Then pass the phylo object to setaGranularityScan().")
        }
    }
    
    latent_method <- match.arg(latent_method)
    
    # Ensure counts is numeric matrix (handle table class from setaCounts)
    if (inherits(counts, "table")) {
        counts <- unclass(counts)
        class(counts) <- "matrix"
        mode(counts) <- "numeric"
    } else if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }
    if (!is.numeric(counts)) {
        counts <- matrix(as.numeric(counts), nrow = nrow(counts), ncol = ncol(counts),
                        dimnames = dimnames(counts))
    }
    
    # Check tree tip labels match column names
    taxa_in_counts <- colnames(counts)
    taxa_in_tree <- tree_phylo$tip.label
    
    if (is.null(taxa_in_counts))
        stop("'counts' must have column names (taxa names) that match tree tip labels.")
    
    if (is.null(taxa_in_tree))
        stop("Tree must have tip labels.")
    
    matched <- intersect(taxa_in_counts, taxa_in_tree)
    if (length(matched) == 0)
        stop("No taxa names match between counts matrix and tree.")
    
    if (length(matched) < length(taxa_in_counts)) {
        warning("Some taxa in counts are not in tree. Using only matching taxa.")
        counts <- counts[, matched, drop = FALSE]
    }
    
    if (length(matched) < length(taxa_in_tree)) {
        tree_phylo <- ape::keep.tip(tree_phylo, matched)
    }
    
    # Set default n_clusters if not provided
    n_taxa <- ncol(counts)
    if (is.null(n_clusters)) {
        n_clusters <- 2:n_taxa
    }
    
    # Validate n_clusters
    n_clusters <- sort(unique(as.integer(n_clusters)))
    n_clusters <- n_clusters[n_clusters >= 2 & n_clusters <= n_taxa]
    
    if (length(n_clusters) == 0)
        stop("'n_clusters' must contain values between 2 and ncol(counts).")
    
    # Convert tree to hclust for cutree
    # as.hclust.phylo requires ultrametric tree, so use distance-based approach
    # Compute cophenetic distances and convert to hclust
    cophenetic_dist <- ape::cophenetic.phylo(tree_phylo)
    hclust_tree <- hclust(as.dist(cophenetic_dist), method = "average")
    
    # Initialize output structures
    latent_spaces <- list()
    aggregated_counts_list <- list()
    metrics_rows <- list()
    
    # Handle cell labels if requested
    cell_labels_df <- NULL
    if (return_cell_labels && !is.null(cell_metadata)) {
        if (is.null(taxa_col))
            stop("If 'cell_metadata' is provided, 'taxa_col' must be specified.")
        
        if (!taxa_col %in% colnames(cell_metadata))
            stop("'taxa_col' not found in 'cell_metadata'.")
        
        # Initialize cell labels data.frame
        # Use row names from cell_metadata if available, otherwise create unique IDs
        if (is.null(rownames(cell_metadata)) || any(duplicated(rownames(cell_metadata)))) {
            # Create unique row names if not present or duplicated
            cell_ids <- paste0("cell_", seq_len(nrow(cell_metadata)))
        } else {
            cell_ids <- rownames(cell_metadata)
        }
        
        cell_labels_df <- data.frame(
            taxa = cell_metadata[[taxa_col]],
            stringsAsFactors = FALSE,
            row.names = cell_ids
        )
    }
    
    # Loop over each number of clusters
    for (n_clust in n_clusters) {
        # Cut tree to n_clust clusters
        cluster_assignments <- cutree(hclust_tree, k = n_clust)
        
        # Store cluster assignments for cell labels if requested
        if (!is.null(cell_labels_df)) {
            # Map cluster assignments to cell metadata
            cluster_col <- paste0("n_clusters_", n_clust)
            cell_labels_df[[cluster_col]] <- cluster_assignments[cell_labels_df$taxa]
        }
        
        # Aggregate counts by cluster
        # Create mapping: each taxon -> its cluster
        cluster_map <- cluster_assignments[colnames(counts)]
        
        # Aggregate: sum counts within each cluster
        aggregated <- matrix(0, nrow = nrow(counts), ncol = n_clust)
        colnames(aggregated) <- paste0("Cluster", seq_len(n_clust))
        rownames(aggregated) <- rownames(counts)
        
        for (clust in seq_len(n_clust)) {
            taxa_in_clust <- names(cluster_map)[cluster_map == clust]
            if (length(taxa_in_clust) > 0) {
                aggregated[, clust] <- rowSums(counts[, taxa_in_clust, drop = FALSE])
            }
        }
        
        # Store aggregated counts
        aggregated_counts_list[[paste0("n_clusters_", n_clust)]] <- aggregated
        
        # Apply CLR transform
        clr_result <- setaCLR(aggregated, pseudocount = pseudocount)
        
        # Compute latent space
        latent_result <- setaLatent(clr_result, method = latent_method, dims = latent_dims)
        
        # Store latent space result
        latent_spaces[[paste0("n_clusters_", n_clust)]] <- latent_result
        
        # Extract metrics
        metrics_row <- data.frame(
            n_clusters = n_clust,
            stringsAsFactors = FALSE
        )
        
        if (latent_method == "PCA") {
            if (length(latent_result$varExplained) >= 1) {
                metrics_row$var_explained_pc1 <- latent_result$varExplained[1]
            }
            if (length(latent_result$varExplained) >= 2) {
                metrics_row$var_explained_pc2 <- latent_result$varExplained[2]
            }
            metrics_row$total_var_explained <- sum(latent_result$varExplained[seq_len(latent_dims)])
        } else if (latent_method == "PCoA") {
            if (length(latent_result$varExplained) >= 1) {
                metrics_row$var_explained_pc1 <- latent_result$varExplained[1]
            }
            if (length(latent_result$varExplained) >= 2) {
                metrics_row$var_explained_pc2 <- latent_result$varExplained[2]
            }
            metrics_row$total_var_explained <- sum(latent_result$varExplained[seq_len(latent_dims)])
        } else if (latent_method == "NMDS") {
            if (is.data.frame(latent_result$varExplained) && "Stress" %in% colnames(latent_result$varExplained)) {
                metrics_row$stress <- latent_result$varExplained$Stress
            }
        }
        
        metrics_rows[[length(metrics_rows) + 1]] <- metrics_row
    }
    
    # Combine metrics
    metrics <- do.call(rbind, metrics_rows)
    
    # Prepare output
    result <- list(
        cluster_labels = cell_labels_df,
        latent_spaces = latent_spaces,
        metrics = metrics,
        aggregated_counts = aggregated_counts_list
    )
    
    result
}
