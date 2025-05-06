#' Extract Taxonomic Counts from Various Single Cell Objects
#'
#' Given a long-form \code{data.frame},
#' creates a type-by-sample matrix of cell counts.
#' Users can specify the column names for cell types, samples, and barcodes.
#'
#' @param obj A long-form data.frame. Typically colData from a
#'   SingleCellExperiment or @meta.data from a Seurat object.
#' @param cell_type_col Column name for cell types (default "type")
#' @param sample_col Column name for sample IDs (default "sample")
#' @param bc_col Column name for barcodes (default "bc")
#'        Use `"rownames"` to extract barcodes from row names.
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
#' print(head(cmat))
#' @return A sample-by-celltype matrix of counts.
#' @export
setaCounts <- function(obj,
                       cell_type_col = "type",
                       sample_col   = "sample",
                       bc_col       = "bc") {
  if (!is.data.frame(obj)) {
    stop("Input must be a data.frame.
          If you want to use Seurat metadata or SCE colData,
          convert these to a dataframe and input them directly.")
  }

  # Handle special case for rownames
  if (identical(bc_col, "rownames")) {
    obj$.__barcodes__ <- rownames(obj)
    bc_col <- ".__barcodes__"
  }

  required <- c(bc_col, cell_type_col, sample_col)
  missing  <- setdiff(required, colnames(obj))
  if (length(missing) > 0) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "))
  }

  df  <- unique(obj[, required])
  mat <- as.matrix(table(df[[sample_col]], df[[cell_type_col]]))
  return(mat)
}

#' Build a taxonomy data frame at multiple resolutions
#'
#' setaTaxonomyDF() converts **one long-form metadata data.frame**-typically
#' colData(sce), seu@meta.data, or any frame you already have, into a tidy
#' taxonomy table.  Each row corresponds to a unique value of the *finest*
#' label (the **last** element of `resolution_cols`), and every coarser label
#' sits in its own column.
#'
#' ## What the input must contain
#' * exactly **one row per cell**
#' * at least one **barcode** column (default `"bc"`).
#'   Pass `bc_col = "rownames"` if barcodes live in `rownames(obj)`.
#' * **all** columns listed in `resolution_cols`
#'
#' No `Seurat`/`SingleCellExperiment` objects are accepted here: extract their
#' metadata/colData first, then hand it in as a `data.frame`
#'
#' ## Value
#' A `data.frame` whose **rownames** are the finest label. If any finest label
#' maps to more than one set of coarser labels the function should stop with an
#' informative error.
#'
#' @param obj A data.frame or similar object containing cell metadata.
#' @param resolution_cols A character vector of column names
#'        indicating hierarchical taxonomy (from broad to fine).
#' @param bc_col Optional. The name of the column containing barcodes,
#'                         or "rownames" if they are row names.
#'
#' @examples
#' meta <- data.frame(
#'   bc          = paste0("cell", 1:6),
#'   fine_type   = c("AT1","AT2","AT1","Fib1","Fib1","AT2"),
#'   mid_type    = c("Alv","Alv","Alv","Fib","Fib","Alv"),
#'   broad_type  = c("Epi","Epi","Epi","Stroma","Stroma","Epi")
#' )
#' setaTaxonomyDF(meta,
#'                resolution_cols = c("broad_type","mid_type","fine_type"))
#'
#' ## barcodes can be in rownames with bc_col = "rownames" (as in Seurat Object)
#' rownames(meta) <- meta$bc
#' meta$bc <- NULL
#' setaTaxonomyDF(meta,
#'                resolution_cols = c("broad_type","mid_type","fine_type"),
#'                bc_col = "rownames")
#' @export
setaTaxonomyDF <- function(obj,
                           resolution_cols = c("fine_type",
                                               "mid_type",
                                               "broad_type"),
                           bc_col = "bc") {

  ## --------------------------------------------------------------------- ##
  ## 0) Basic argument sanity -------------------------------------------- ##
  if (!is.data.frame(obj))
    stop("obj must be a data.frame. Extract metadata first.")

  if (!is.character(resolution_cols) || length(resolution_cols) < 1)
    stop("resolution_cols must be a non-empty character vector.")

  ## --------------------------------------------------------------------- ##
  ## 1) Handle barcode column (incl. 'rownames' shortcut) ---------------- ##
  if (identical(bc_col, "rownames")) {
    obj$.__bc__ <- rownames(obj)
    bc_col <- ".__bc__"
  }

  req_cols <- c(bc_col, resolution_cols)
  miss     <- setdiff(req_cols, names(obj))
  if (length(miss))
    stop("Missing required column(s): ", paste(miss, collapse = ", "))

  meta <- obj[, req_cols, drop = FALSE]

  ## --------------------------------------------------------------------- ##
  ## 2) Validate the finest label ---------------------------------------- ##
  fine_label <- tail(resolution_cols, 1)
  if (anyNA(meta[[fine_label]]))
    stop("Some cells have NA in the finest label column ('", fine_label, "').")

  ## --------------------------------------------------------------------- ##
  ## 3) Unique combos & one to one mapping check ------------------------- ##
  combos <- unique(meta[, resolution_cols, drop = FALSE])

  bad <- names(which(table(combos[[fine_label]]) > 1))
  if (length(bad))
    stop("Finest labels mapping to >1 coarser combo: ",
         paste(bad, collapse = ", "))

  ## --------------------------------------------------------------------- ##
  ## 4) Return taxonomy frame ------------------------------------------- ##
  rownames(combos) <- combos[[fine_label]]
  return(combos)
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
