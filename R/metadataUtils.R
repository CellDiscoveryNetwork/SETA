#' Extract Sample-Level Metadata from Various Objects
#'
#' This function extracts sample-level metadata from a given source object (Seurat, SingleCellExperiment, 
#' or dataframe). It ensures that each metadata column contains unique values per sample. If a metadata 
#' column contains non-unique values within any sample, that column is excluded from the output, and the 
#' user is notified via a warning. This function is particularly useful when preparing metadata for 
#' high-level visualizations or analyses where sample-level aggregation is required.
#'
#' @param source_obj An object of class Seurat, SingleCellExperiment, or a dataframe, which contains
#' cell-level or sample-level metadata.
#' @param sample_col Character. The name of the column in source_obj that indicates the sample identifier. 
#' Default is 'sample_id'. This column is used to group the metadata and check for uniqueness.
#' @param meta_cols Character vector. Names of metadata columns to retain. If NULL, all columns present in 
#' the source object are considered. However, only those columns where all entries are identical within 
#' each sample are included in the final output.
#'
#' @return A dataframe where each row corresponds to a unique sample and each column represents a metadata 
#' variable that has uniform values within samples. Columns with non-unique values within any sample are 
#' excluded, and a warning lists these columns.
#'
#' @examples
#' # Using a Seurat object
#' \donttest{
#' meta_df <- setaMetadata(seurat_obj,
#'                         sample_col="donor_id",
#'                         meta_cols=c("disease", "Severity"))
#'
#' # Using a SingleCellExperiment object with default parameters
#' meta_df <- setaMetadata(sce_obj)
#'
#' # Using a dataframe and extracting all possible metadata columns
#' meta_df <- setaMetadata(data_frame)
#' }
#' @export
setaMetadata <- function(source_obj, sample_col = "Sample ID", meta_cols = NULL) {
  if (inherits(source_obj, "Seurat")) {
    meta_df <- source_obj@meta.data
  } else if (inherits(source_obj, "SingleCellExperiment")) {
    meta_df <- as.data.frame(SummarizedExperiment::colData(source_obj))
  } else if (is.data.frame(source_obj)) {
    meta_df <- source_obj
  } else {
    stop("Unsupported object type. Must be Seurat, SingleCellExperiment, or dataframe.")
  }
  
  if (!sample_col %in% colnames(meta_df))
    stop("sample_col not found in metadata columns.")
  
  if (!is.null(meta_cols)) {
    if (!all(meta_cols %in% colnames(meta_df))) {
      missing_cols <- setdiff(meta_cols, colnames(meta_df))
      stop("The following meta_cols are not in metadata: ", paste(missing_cols, collapse = ", "))
    }
    meta_df <- meta_df[, c(sample_col, meta_cols), drop = FALSE]
  }
  
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
    warning("Non-unique values in: ", paste(nu_cols[nu_cols != sample_col], collapse = ", "))
    summarized_meta <- summarized_meta[, !(colnames(summarized_meta) %in% nu_cols[nu_cols != sample_col])]
  }
  
  colnames(summarized_meta)[1] <- "Sample ID"
  return(summarized_meta)
}
