#' Extract Taxonomic Counts from Various Single Cell Objects
#'
#' Given a \code{SingleCellExperiment}, \code{Seurat}, or
#' long-form \code{data.frame}, creates a type-by-sample matrix of cell counts.
#' Users can specify the column names for cell types, samples, and barcodes.
#'
#' @param obj Either a \code{SingleCellExperiment}, a \code{Seurat} object, or
#'   a \code{data.frame}. This function requires specific metadata or columns
#'   to function correctly.
#' @param cell_type_col The name of the column representing cell types.
#'   Default is "type".
#' @param sample_col The name of the column representing sample identifiers.
#'   Default is "sample".
#' @param bc_col The name of the column representing barcodes (only needed for
#'   \code{data.frame} input). Default is "bc".
#'
#' @return A matrix whose rows are samples and whose columns are cell types,
#'   with entries of the count of unique barcodes per type-sample combination.
#'
#' @details
#' \itemize{
#'   \item \strong{SingleCellExperiment} Reads \code{colData(obj)}
#'   for the specified \code{cell_type_col} and \code{sample_col}.
#'   \item \strong{Seurat} Uses \code{obj@meta.data} for the specified
#'   \code{cell_type_col} and \code{sample_col}.
#'   \item \strong{data.frame} Counts cells per specified \code{cell_type_col}
#'   and \code{sample_col}.
#' }
#' If the specified columns are missing, an error is thrown.
#'
#' @examples
#' # For a data.frame with custom column names:
#' set.seed(687)
#' df <- data.frame(
#'   barcode = paste0("cell", 1:10),
#'   cellType = sample(c("Tcell", "Bcell"), 10, TRUE),
#'   sampleID = sample(c("sample1","sample2"), 10, TRUE)
#' )
#' cmat <- setaCounts(df,
#'                    cell_type_col = "cellType",
#'                    sample_col = "sampleID",
#'                    bc_col = "barcode")
#' print(cmat)
#'
#'
#' @export
setaCounts <- function(obj,
                       cell_type_col = "type",
                       sample_col = "sample",
                       bc_col = "bc") {
  # SingleCellExperiment import
  if ("SingleCellExperiment" %in% class(obj)) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' must be installed.", call. = FALSE)
    }
    coldata <- SingleCellExperiment::colData(obj)
    if (!all(c(cell_type_col, sample_col) %in% colnames(coldata))) {
      stop(sprintf("colData(obj) must contain '%s' and '%s' columns.",
                   cell_type_col,
                   sample_col),
           call. = FALSE)
    }
    return(as.matrix(table(coldata[[sample_col]], coldata[[cell_type_col]])))
  }

  # Seurat import
  if ("Seurat" %in% class(obj)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' must be installed.", call. = FALSE)
    }
    meta <- obj@meta.data
    if (!all(c(cell_type_col, sample_col) %in% colnames(meta))) {
      stop(sprintf("obj@meta.data must contain '%s' and '%s' columns.",
                   cell_type_col,
                   sample_col),
           call. = FALSE)
    }
    return(as.matrix(table(meta[[sample_col]], meta[[cell_type_col]])))
  }

  # Long-form data.frame import
  if (is.data.frame(obj)) {
    requiredCols <- c(bc_col, cell_type_col, sample_col)
    if (!all(requiredCols %in% colnames(obj))) {
      stop(sprintf("data.frame must have columns: '%s', '%s', '%s' at minimum.",
                   bc_col,
                   cell_type_col,
                   sample_col),
           call. = FALSE)
    }
    # Deduplicate by barcode so each cell is counted once
    df <- unique(obj[, c(bc_col, cell_type_col, sample_col)])
    return(as.matrix(table(df[[sample_col]], df[[cell_type_col]])))
  }
  stop("Unsupported object type.")
}

#' Build a Taxonomy Data Frame at Multiple Resolutions (Base R version)
#'
#' Given a Seurat object, SingleCellExperiment, or long-form \code{data.frame},
#' this function constructs a "taxonomy" data frame. Each row corresponds to
#' a unique "lowest-level" label, and each column is one resolution level
#' specified in \code{resolution_cols}. The last element of
#' \code{resolution_cols} is considered the "lowest-level" (finest) label,
#' which will become the row name of the returned data frame.
#'
#' @param obj One of:
#'   \itemize{
#'     \item A \strong{Seurat} object, with \code{@meta.data} containing columns
#'       named in \code{resolution_cols}.
#'     \item A \strong{SingleCellExperiment} object, with \code{colData()}
#'       containing columns named in \code{resolution_cols}.
#'     \item A long-form \code{data.frame}, with one row per cell, containing
#'       at least a \code{bc} column plus all columns in \code{resolution_cols}.
#'   }
#' @param resolution_cols A character vector of colnames indicating different
#'   resolutions of cell-type (or lineage) labels. The \strong{last} element in
#'   this vector is treated as the "lowest-level" (finest) label, which will
#'   become the row name of the returned data frame.
#'
#' @return A \code{data.frame} where each row is one unique value of the
#'   "lowest-level" label. The columns are all the entries in
#'   \code{resolution_cols}. Row names are set to the lowest-level label.
#'   If a single lowest-level label maps to multiple coarser-level labels,
#'   the function throws an error.
#'
#' @examples
#' # Long-form example
#' df_long <- data.frame(
#'     bc = paste0("cell", 1:6),
#'     fine_type = c("AlveolarType1","AlveolarType2","AlveolarType1",
#'                   "Fibroblast1","Fibroblast1","AlveolarType2"),
#'     mid_type  = c("Alveolar","Alveolar","Alveolar","Fibroblast",
#'                   "Fibroblast","Alveolar"),
#'     broad_type= c("Epithelial","Epithelial","Epithelial","Stromal",
#'                   "Stromal","Epithelial")
#' )
#'
#' # Build a taxonomy data frame specifying fine->mid->broad hierarchy
#' taxDF <- setaTaxonomyDF(
#'     df_long,
#'     resolution_cols = c("fine_type","mid_type","broad_type")
#' )
#' taxDF
#'
#' @export
setaTaxonomyDF <- function(
    obj,
    resolution_cols = c("fine_type", "mid_type", "broad_type")
) {
    ########################################################################
    ## 1) Extract or build a data.frame with columns: bc + resolution_cols ##
    ########################################################################
    if (!is.vector(resolution_cols) || length(resolution_cols) < 1) {
        stop("'resolution_cols' must be a character vector
              with at least 1 entry.")
    }

    meta <- NULL

    # CASE A: Long-form data.frame
    if (is.data.frame(obj)) {
        # Must have at least 'bc' plus all resolution_cols
        req_cols <- c("bc", resolution_cols)
        missing_cols <- setdiff(req_cols, colnames(obj))
        if (length(missing_cols) > 0) {
            stop(
                "Long-form data.frame is missing columns: ",
                paste(missing_cols, collapse = ", ")
            )
        }
        meta <- obj[, req_cols, drop = FALSE]

    # CASE B: Seurat
    } else if (inherits(obj, "Seurat")) {
        meta_seu <- obj@meta.data
        missing_cols <- setdiff(resolution_cols, colnames(meta_seu))
        if (length(missing_cols) > 0) {
            stop(
                "Seurat object is missing these columns in @meta.data: ",
                paste(missing_cols, collapse = ", ")
            )
        }
        # We'll add a 'bc' column from rownames
        meta_seu$bc <- rownames(meta_seu)
        meta <- meta_seu[, c("bc", resolution_cols), drop = FALSE]

    # CASE C: SingleCellExperiment
    } else if (inherits(obj, "SingleCellExperiment")) {
        cd <- SummarizedExperiment::colData(obj)
        missing_cols <- setdiff(resolution_cols, colnames(cd))
        if (length(missing_cols) > 0) {
            stop(
                "SingleCellExperiment colData is missing columns: ",
                paste(missing_cols, collapse = ", ")
            )
        }
        # Build a base df: bc is colnames(obj), then resolution_cols
        bc_vec <- colnames(obj)
        meta_sce <- data.frame(
            bc = bc_vec,
            as.data.frame(cd[, resolution_cols, drop = FALSE]),
            stringsAsFactors = FALSE
        )
        meta <- meta_sce

    } else {
        stop(
            "Unsupported object type. Must be data.frame, ",
            "Seurat, or SingleCellExperiment."
        )
    }

    ############################################################################
    ## 2) Identify the fine (lowest-level) label last in resolution_cols      ##
    ############################################################################
    fine_label <- tail(resolution_cols, 1)

    if (any(is.na(meta[[fine_label]]))) {
        stop(
            "Some cells have an NA value for the '",
            fine_label,
            "' column."
        )
    }

    ########################################################################
    ## 3) Construct unique combos of all resolution_cols -> one row per   ##
    ##    distinct (fine_type, mid_type, broad_type, etc.) combo.         ##
    ########################################################################
    needed <- meta[, resolution_cols, drop = FALSE]

    unique_combos <- unique(needed)

    ########################################################################
    ## 4) Check that each fine_label is associated with exactly ONE combo ##
    ########################################################################
    fine_vals <- unique_combos[[fine_label]]
    freq_tab <- table(fine_vals)
    conflicts <- names(freq_tab)[freq_tab > 1]
    if (length(conflicts) > 0) {
        stop(
            "Some fine-level labels map to multiple coarser combos: ",
            paste(conflicts, collapse = ", ")
        )
    }

    ########################################################################
    ## 5) Sort the table by coarsest label,                               ##
    ##    Set row names to the fine_label and return the unique_combos df ##
    ########################################################################
    # unique_combos <- unique_combos[order(unique_combos[[1]]), ]
    rownames(unique_combos) <- unique_combos[[fine_label]]
    unique_combos
}

#' Convert Multi-Column Taxonomy to a Single-Root tbl_graph (with node metadata)
#'
#' This function takes a data frame describing a hierarchical taxonomy across
#' multiple columns (e.g., broad -> mid -> fine). Each row represents a unique
#' path through the hierarchy. The function introduces a single root node
#' (named \code{root_name}) above the first hierarchy column, then constructs
#' a directed tree in which each level connects to the next.
#' After building the graph, it appends node-level metadata by looking up
#' which rows (and columns) in \code{tax_df} contain each node. This allows
#' you to color or facet by different levels of the taxonomy when using \pkg{ggraph}.
#'
#' @param tax_df A data frame with one row per unique path in the hierarchy.
#'   For example, if your columns are \code{c(\"broad\",\"mid\",\"fine\")}, each
#'   row is a single path from \code{broad -> mid -> fine}.
#' @param columns A character vector of column names in \code{tax_df} to use.
#'   They should be ordered from the broadest level (first) to the finest level
#'   (last). If \code{NULL}, the function will use all columns of
#'   \code{tax_df} in their given order.
#' @param root_name A character string naming the artificial root node,
#'   inserted above the first hierarchy column. Default is \code{\"AllCells\"}.
#'
#' @return A \code{tbl_graph} object (directed) with a single root node. The node
#'   data includes extra columns corresponding to each level in \code{columns}.
#'   If a node corresponds to multiple categories at a given level, these are
#'   combined with \code{\"|\"}.
#'
#' @details
#' 1. The function first builds an edge list
#'    \enumerate{
#'      \item \code{Root -> level1} for each row
#'      \item \code{level1 -> level2}
#'      \item \ldots
#'      \item \code{level_{N-1} -> levelN}
#'    }
#'    and removes duplicates, creating a single connected tree.
#'
#' 2. It then \emph{annotates each node} with the best-known taxonomy data. For
#'    a node named \code{x}, we look up all rows of \code{tax_df} where
#'    \code{x} appears in \code{columns}, gather the distinct values from each
#'    \code{col}, and store them joined with \code{\"|\"} if more than one
#'    distinct value is found.
#'
#' This means if a node is shared among multiple broad categories (uncommon, but
#' possible), that node's \code{broad} column will contain something like
#' \code{\"Epithelial|Stromal\"}.
#'
#' @examples
#' # Minimal example with a 3-level hierarchy (broad -> mid -> fine)
#' tax_df_example <- data.frame(
#'     broad = c(\"Epithelial\", \"Epithelial\", \"Stromal\"),
#'     mid   = c(\"Alveolar\", \"Alveolar\", \"Fibroblast\"),
#'     fine  = c(\"AlveolarType1\", \"AlveolarType2\", \"Fibroblast1\"),
#'     stringsAsFactors = FALSE
#' )
#'
#' library(tidygraph)
#' library(ggraph)
#' library(ggplot2)
#'
#' # Build a single-root tree and incorporate node metadata
#' tbl_g <- taxonomy_to_tbl_graph(
#'     tax_df_example,
#'     columns   = c(\"broad\", \"mid\", \"fine\"),
#'     root_name = \"AllCells\"
#' )
#'
#' # Inspect node data (metadata for each node)
#' as.data.frame(tbl_g, \"nodes\")
#'
#' # Visualize with ggraph, coloring by 'broad' level
#' ggraph(tbl_g, layout = \"tree\") +
#'     geom_edge_diagonal() +
#'     geom_node_point(aes(color = broad), size = 3) +
#'     geom_node_text(aes(label = name), vjust = 1, hjust = 0.5) +
#'     theme_minimal() +
#'     labs(title = \"Single-Root Taxonomy Tree\")
#'
#' @export
taxonomy_to_tbl_graph <- function(tax_df,
                                  columns   = NULL,
                                  root_name = "AllCells") {
    if (!requireNamespace("tidygraph", quietly = TRUE)) {
        stop("The 'tidygraph' package is required for taxonomy_to_tbl_graph().")
    }
    if (is.null(columns)) {
        columns <- colnames(tax_df)
    }
    if (length(columns) < 1) {
        stop("Need at least one column in 'columns'.")
    }

    ############################
    ## 1) Build the edge list ##
    ############################
    edge_list <- data.frame(
        from = character(0),
        to   = character(0),
        stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(tax_df))) {
        path <- c(root_name, as.character(tax_df[i, columns, drop = TRUE]))
        # e.g. c("AllCells", "Epithelial", "Alveolar", "AlveolarType1")

        # Build edges between consecutive elements in path
        for (j in seq_along(path)[-length(path)]) {
            edge_list <- rbind(
                edge_list,
                data.frame(from = path[j], to = path[j+1],
                           stringsAsFactors = FALSE)
            )
        }
    }
    # Remove duplicates
    edge_list <- unique(edge_list)

    # Create tbl_graph
    tg <- tidygraph::as_tbl_graph(edge_list, directed = TRUE)

    ###############################
    ## 2) Add node-level metadata
    ###############################
    # We'll define columns for each taxonomy level. For a node 'x'
    #  - if x == root_name, we store NA or "Root"
    #  - else find all rows of tax_df where x appears in row's columns,
    #    gather distinct values for each col, combine with '|'.

    node_names <- tg %>% tidygraph::activate("nodes") %>% dplyr::pull(name)
    # We'll build a data frame with the same number of rows as node_names
    node_info <- data.frame(name = node_names, stringsAsFactors = FALSE)

    for (col in columns) {
        node_info[[col]] <- vapply(
            node_info$name,
            FUN.VALUE = character(1),
            FUN = function(x) {
                if (identical(x, root_name)) {
                    # Root node => NA
                    return(NA_character_)
                }
                # Which rows contain 'x' in column set?
                # E.g., if x == "AlveolarType1" or "Epithelial"
                idx <- which(apply(tax_df[, columns, drop=FALSE], 1,
                                  function(rowvals) x %in% rowvals))
                if (length(idx) == 0) {
                    # Possibly an error or leftover node that wasn't in tax_df
                    return(NA_character_)
                }
                # Gather distinct values of tax_df[idx, col]
                vals <- unique(tax_df[idx, col])
                if (length(vals) > 1) {
                    # If node is used in multiple places with different col entries
                    return(paste(vals, collapse="|"))
                }
                vals
            }
        )
    }

    # Bind these columns into the node data of tbl_graph
    # We'll do: tg <- tg %>% tidygraph::left_join(node_info, by="name")
    # But 'left_join' is from dplyr, so we can do it with 'bind_nodes'
    # from tidygraph 1.2 or we can do a direct manual approach:
    tg <- tg %>%
        tidygraph::activate("nodes") %>%
        dplyr::left_join(node_info, by="name")

    tg
}
