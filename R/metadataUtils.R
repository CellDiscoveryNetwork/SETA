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
#' meta_df <- setaMetadata(seurat_obj@meta.data,
#'                         sample_col="donor_id",
#'                         meta_cols=c("disease", "Severity"))
#'
#' # Using a SingleCellExperiment object with default parameters
#' meta_df <- setaMetadata(data.frame(colData(sce_obj)))
#'
#' # Using a dataframe and extracting all possible metadata columns
#' meta_df <- setaMetadata(df)
#' }
#' @export
setaMetadata <- function(x,
                         sample_col = "Sample ID",
                         meta_cols = NULL) {
    stopifnot(
        is.data.frame(x), !is.null(meta_cols),
        sample_col %in% colnames(x),
        is.character(sample_col), is.character(meta_cols)
        )
    
    if (!all(meta_cols %in% colnames(x))) {
        missing_cols <- setdiff(meta_cols, colnames(x))
        stop("The following meta_cols are not in metadata: ",
             paste(missing_cols, collapse = ", "))
    }
    meta_df <- x[, c(sample_col, meta_cols), drop = FALSE]
    
    # Group by sample_col and compute summary keeping original types
    grouped <- split(meta_df, meta_df[[sample_col]])
    summarized_list <- lapply(grouped, function(df) {
        s <- lapply(df, function(x) {
            u <- unique(x)
            if (length(u) == 1) u else NA
        })
        as.data.frame(s, stringsAsFactors = FALSE)
    })
    summarized_meta <- do.call(rbind, summarized_list)
    # Drop columns with NA values (non-unique across the group)
    nu_cols <- colnames(summarized_meta)[colSums(is.na(summarized_meta)) > 0]
    if (length(nu_cols) > 1) {
        warning("Non-unique values in: ",
                paste(nu_cols[nu_cols != sample_col], collapse = ", "))
        keep <- !(colnames(summarized_meta) %in% nu_cols[nu_cols != sample_col])
        summarized_meta <- summarized_meta[, keep]
    }
    colnames(summarized_meta)[1] <- "sample_id"
    return(summarized_meta)
}
