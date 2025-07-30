#' Extract Sample-Level Metadata from Various Objects
#'
#' This function extracts sample-level metadata
#' from a dataframe. It ensures that each metadata column contains
#' unique values per sample. If a metadata column contains
#' non-unique values within any sample, that column is excluded from the output,
#' and the user is notified via a warning. Useful when preparing metadata for
#' visualizations or analyses where sample-level inspection is required.
#'
#' @param x An object of class dataframe which contains cell-level or
#' sample-level metadata.
#' @param sample_col Character. The sample identifier column name in source_obj.
#' Default is 'sample_id'. This column is used to group the metadata.
#' @param meta_cols Character vector. Names of metadata columns to retain.
#' If NULL, all columns present in the source object are considered.
#' However, only those columns where all entries are identical within
#' each sample are included in the final output.
#'
#' @return A dataframe where each row corresponds to a unique sample and each
#' column represents a metadata variable that has uniform values within samples.
#' Columns with non-unique values within any sample are
#' excluded, and a warning lists these columns.
#'
#' @examples
#' # Using a Seurat object
#' \donttest{
#' # if (requireNamespace("SeuratObject", quietly = TRUE)) {
#' # meta_df <- setaMetadata(seurat_obj@meta.data,
#' #                         sample_col="donor_id",
#' #                         meta_cols=c("disease", "Severity"))
#' # }
#' #
#' # Using a SingleCellExperiment object with default parameters
#' # if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#' # meta_df <- setaMetadata(data.frame(colData(sce_obj)))
#' #
#' # Using a dataframe and extracting all possible metadata columns
#' # meta_df <- setaMetadata(df)
#' # }
#' }
#' @export
setaMetadata <- function(x,
                         sample_col = "Sample ID",
                         meta_cols) {
    stopifnot(
        is.data.frame(x),
        length(sample_col) == 1L, is.character(sample_col),
        sample_col %in% names(x),
        is.character(meta_cols), length(meta_cols) > 0L
    )
    if (anyNA(x[[sample_col]])) stop("`", sample_col, "` contains NA.")
    
    meta_cols <- setdiff(unique(meta_cols), sample_col)
    missing <- setdiff(meta_cols, names(x))
    if (length(missing)) stop(
        "The following meta_cols are not in your data: ",
        paste(missing, collapse = ", "),
        call. = FALSE
    )
    
    samples <- unique(as.character(x[[sample_col]]))
    out <- data.frame(sample_id = samples, stringsAsFactors = FALSE)
    if (is.factor(x[[sample_col]])) {
        out$sample_id <- factor(out$sample_id,
                                levels = levels(x[[sample_col]]))
    }
    
    for (col in meta_cols) {
        col_data <- x[[col]]
        vals <- vapply(samples, function(sid) {
            v <- col_data[x[[sample_col]] == sid]
            u <- unique(v)
            if (length(u) != 1L) {
                stop(sprintf(
                    "Column '%s' has multiple values for sample '%s': %s.\n
                    Are your samples multiplexed? 
                    If so, please supply a sample identifier unique 
                    to each sample X pool.",
                    col, sid, paste(u, collapse = ", ")
                ), call. = FALSE)
            }
            as.character(u)
        }, FUN.VALUE = character(1))
        # preserve original type
        if (is.factor(col_data)) {
            out[[col]] <- factor(vals, levels = levels(col_data))
        } else if (is.integer(col_data)) {
            out[[col]] <- as.integer(vals)
        } else if (is.numeric(col_data)) {
            out[[col]] <- as.numeric(vals)
        } else if (is.logical(col_data)) {
            out[[col]] <- as.logical(vals)
        } else {
            out[[col]] <- vals
        }
    }
    
    out
}