#' Centered Log-Ratio (CLR) Transform
#' Applies a CLR transform to a matrix of counts.
#' Samples should be in rows and taxa (cell types) in columns.
#' For each sample, the transform computes
#' \eqn{\mathrm{CLR}(x)_i = \log\left(\frac{x_i + \text{pseudocount}}{g(x + \text{pseudocount})}\right)},
#' where \eqn{g(x + \text{pseudocount})} is the geometric mean of the row.
#'
#' @param counts An integer matrix of celltype counts with samples in rows.
#' @param pseudocount Numeric.
#'        Added to all entries to avoid \code{log(0)}. Default is 1.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{A string indicating the transform ("CLR").}
#'   \item{counts}{A matrix of the same dimensions as the input after CLR transform.}
#' }
#'
#' @details
#' The CLR transform is defined sample-wise as:
#' \deqn{\mathrm{CLR}(x)_{ij} = \log\left(\frac{x_{ij} + \text{pseudocount}}{g_i}\right)}{
#' \log\left(\frac{x_{ij} + \text{pseudocount}}{g_i}\right)}
#' where \eqn{g_i = \exp\left(\frac{1}{p}\sum_{j=1}^{p}\log(x_{ij} + \text{pseudocount})\right)} for sample \(i\)
#' and \(p\) is the number of taxa.
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society.
#' Series B (Methodological)}, 44(2), 139-177.
#'
#' @examples
#' # Example matrix with 2 samples and 2 taxa:
#' mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("Taxon1", "Taxon2")
#' out <- setaCLR(mat, pseudocount = 0)
#' out$counts
#' @name setaCLR
#' @export
setaCLR <- function(counts, pseudocount = 1) {
    if (!is.matrix(counts)) stop("'counts' must be a matrix.")
    counts <- counts + pseudocount
    # Compute geometric mean for each sample (row)
    gm <- exp(apply(log(counts), 1, mean))
    # Subtract log(geometric mean) from each log-transformed element in each row
    clr_mat <- sweep(log(counts), 1, log(gm), FUN = "-")
    # Restore column names
    colnames(clr_mat) <- colnames(counts)
    list(method = "CLR", counts = clr_mat)
}

#' Isometric Log-Ratio (ILR) Transform
#' Applies the ILR transform to an integer counts matrix.
#' For each sample (row), the data are log-transformed
#' (with an optional Box Cox like transformation)
#' then projected onto an orthonormal Helmert basis,
#' reducing dimensionality by one.
#'
#' @param counts An integer matrix of celltype counts with samples in rows.
#' @param boxcox_p Numeric. If nonzero, a Box Cox type transform
#'   is applied to the log-values. Default is 0 (no Box Cox transformation).
#' @param taxTree Unused. Reserved for future taxonomic-balance approaches.
#' @param pseudocount Numeric. Added to avoid \code{log(0)}. Default is 1.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{A string indicating the ILR transform.
#'                 If \code{boxcox_p} is nonzero,
#'                 the value is indicated in the method string.}
#'   \item{counts}{A matrix of ILR-transformed values with
#'                 \code{ncol(counts) - 1} columns
#'                 and the same number of rows (samples) as the input.}
#' }
#'
#' @details
#' The ILR transform is computed as follows:
#' \enumerate{
#'   \item Add a pseudocount and take the natural logarithm:
#'     \deqn{y = \log(x + \text{pseudocount})}
#'   \item If \code{boxcox_p != 0}, apply the Box Cox like transform:
#'     \deqn{y = \frac{\exp(p \, y) - 1}{p}}
#'   \item Project the log-transformed data onto an orthonormal
#'         Helmert basis computed via QR decomposition.
#' }
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society.
#'       Series B (Methodological)}, 44(2), 139-177.
#'
#' @examples
#' # Example matrix: rows are samples, columns are cell types.
#' mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("A", "B")
#' # ILR transformation reduces the dimension by 1.
#' out <- setaILR(mat, boxcox_p = 0, pseudocount = 1)
#' out$counts
#' @name setaILR
#' @export
setaILR <- function(counts, boxcox_p = 0, taxTree = NULL, pseudocount = 1) {
    if (!is.matrix(counts)) stop("'counts' must be a matrix.")
    if (!is.null(taxTree)) {
        message("A taxTree was provided but is not yet supported. 
             Defaulting to Helmert basis.")
    }
    counts <- counts + pseudocount
    log_x <- log(counts)
    if (boxcox_p != 0) {
        log_x <- (exp(boxcox_p * log_x) - 1) / boxcox_p
    }
    n_taxa <- ncol(counts)
    if (n_taxa < 2) stop("ILR requires at least 2 taxa (columns).")
    
    # Compute orthonormal Helmert basis
    H <- stats::contr.helmert(n_taxa)
    H <- qr.Q(qr(H))  # orthonormalize via QR decomposition
    # Directly project log-values onto the orthonormal Helmert basis
    ilr_mat <- log_x %*% H
    
    list(
        method = paste0("ILR_Helmert",
                        ifelse(boxcox_p != 0, 
                               paste0(" (boxcox_p=", boxcox_p, ")"),
                               "")
        ),
        counts = ilr_mat
    )
}

#' Additive Log-Ratio (ALR) Transform
#' Applies the ALR transform to an integer matrix of counts
#' using a specified reference taxon. Samples are in rows and taxa in columns.
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#' @param ref Either the reference taxon name (a character string,
#'            which must appear in \code{colnames(counts)})
#'            or the column index of the reference.
#' @param pseudocount Numeric.
#'        Added to every count to avoid \code{log(0)}. Default is 1.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{A string indicating the ALR transform
#'                 with the reference taxon.}
#'   \item{counts}{A matrix with one row per sample
#'                 and \eqn{(\text{n_taxa} - 1)} columns.}
#' }
#'
#' @details
#' Applies the ALR transform to an integer matrix of counts
#' using a specified reference taxon. Samples are in rows and taxa in columns.
#' For each sample, the transform computes:
#' \deqn{\mathrm{ALR}(x)_i = \log\left(\frac{x_i + \text{pseudocount}}{x_{ref} + \text{pseudocount}}\right)}{
#' \log\left(\frac{x_i + \text{pseudocount}}{x_{ref} + \text{pseudocount}}\right)}
#' for all taxa \(i\) except the reference.
#'
#' @examples
#' # Example with 2 samples and 2 taxa:
#' mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("TaxonA", "TaxonB")
#' # Using TaxonA as the reference.
#' out <- setaALR(mat, ref = "TaxonA", pseudocount = 0)
#' out$counts
#' @name setaALR
#' @export
setaALR <- function(counts, ref, pseudocount = 1) {
    if (!is.matrix(counts)) stop("'counts' must be a matrix.")
    if (missing(ref)) stop("Please specify a reference taxon (by name or index).")
    if (is.character(ref)) {
        if (!(ref %in% colnames(counts))) {
            stop("Reference taxon not found in colnames(counts).")
        }
        refCol <- which(colnames(counts) == ref)
    } else if (is.numeric(ref)) {
        if (ref < 1 || ref > ncol(counts)) {
            stop("Reference taxon index out of range.")
        }
        refCol <- ref
    } else {
        stop("'ref' must be a character or numeric.")
    }
    counts <- counts + pseudocount
    # For each sample, subtract log(value at reference taxon) from log(counts)
    alr_mat <- sweep(log(counts), 1, log(counts[, refCol]), FUN = "-")
    # Remove the reference taxon column from output
    alr_mat <- alr_mat[, -refCol, drop = FALSE]
    list(method = paste0("ALR_ref=",
                         ifelse(is.character(ref),
                                ref,
                                colnames(counts)[refCol])
    ),
    counts = alr_mat)
}

#' Percentage Transform
#' Converts each row (sample) of a counts matrix to percentages of its row sum.
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{The string \code{"percent"}.}
#'   \item{counts}{A matrix of the same dimensions as
#'                 \code{counts}, where each row sums to 100.}
#' }
#'
#' @details
#' Useful for simplified comparisons and as an input to non-parametric tests.
#'
#' @examples
#' mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
#' out <- setaPercent(mat)
#' out$counts
#' @name setaPercent
#' @export
setaPercent <- function(counts) {
    if (!is.matrix(counts)) stop("'counts' must be a matrix.")
    pct_counts <- sweep(counts, 1, rowSums(counts), FUN = "/") * 100
    list(method = "percent", counts = pct_counts)
}

#' log2(CPM) Transform
#' Computes the log2 counts-per-million (CPM) for each sample.
#' Samples are in rows and taxa in columns.
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#' @param pseudocount Numeric.
#'        Added to counts to avoid \code{log2(0)}. Default is 1.
#' @param size_factors Optional numeric vector of library sizes for each sample.
#'        If \code{NULL}, the row sums are used.
#' @param scale_factor Numeric.
#'        The scaling factor, typically 1e6 for CPM. Default is 1e6.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{The string \code{"logCPM"}.}
#'   \item{counts}{A matrix of the same dimensions
#'                 with log2-transformed CPM values.}
#' }
#'
#' @details
#' The transform is defined as:
#' \deqn{\log_2\left(\frac{(x + \text{pseudocount})}{\text{library size}} \times \text{scale_factor}\right),}{
#' \log_2\left(\frac{(x + \text{pseudocount})}{\text{library size}} \times \text{scale_factor}\right)}
#' where the library size is computed per sample.
#' @examples
#' mat <- matrix(c(10, 20, 100, 200), nrow = 2, byrow = TRUE)
#' out <- setaLogCPM(mat, pseudocount = 1)
#' out$counts
#' @name setaLogCPM
#' @export
setaLogCPM <- function(counts,
                       pseudocount = 1,
                       size_factors = NULL,
                       scale_factor = 1e6) {
    if (!is.matrix(counts)) stop("'counts' must be a matrix.")
    if (is.null(size_factors)) size_factors <- rowSums(counts)
    cpm <- sweep(counts + pseudocount,
                 1, 
                 size_factors + pseudocount,
                 FUN = "/") * scale_factor
    log_cpm <- log2(cpm)
    list(method = "logCPM", counts = log_cpm)
}

#' User‑defined balance transform (geometric‑mean log‑ratio)
#'
#' `setaBalance()` computes *one or more* biologically meaningful balances
#' (log‑ratios) from a count matrix.  Each balance is defined by two
#' groups of taxa: **numerator** (`num`) and **denominator** (`denom`).
#' Groups may be given as leaf names, higher‑level labels (resolved through a
#' `taxonomyDF`), or column indices. The resulting balance will be positive
#' if weighted in the numerator direction, and negative toward the denominator
#'
#' For every balance and every sample the function returns
#' \deqn{\log \frac{\mathrm{GM}(\mathrm{num})}{\mathrm{GM}(\mathrm{denom})}}
#' where GM() is the geometric mean of the (pseudocount‑adjusted) counts in
#' the respective group.
#'
#' @param counts Numeric matrix with **rows = samples, columns = taxa**.
#' @param balances A *single* balance (list with `num`, `denom`) **or**
#'        a *named list* of such lists for multiple balances.
#' @param taxonomyDF Optional.  A data frame from [setaTaxonomyDF()] used to
#'        expand higher‑level labels into their descendant leaves.
#' @param taxonomy_col Character.  Column in `taxonomyDF` whose values should
#'        match any higher‑level labels given in `balances`.
#' @param normalize_to_parent Logical (default `FALSE`). If `TRUE`, each sample
#'        is re‑closed to the sub‑composition formed by `num ∪ denom` before
#'        taking the log‑ratio – i.e. the balance is within the parent total.
#' @param pseudocount Numeric. Value added to every count to avoid
#'        `log(0)`.  Default `1`.
#'
#' @return A list with
#' \describe{
#'   \item{method}{`"balance"`}
#'   \item{counts}{Matrix **samples × balances**.  Column names are the
#'                 balance names (or `"Balance1"` if unnamed).}
#' }
#'
#' @examples
#' ## Toy metadata & taxonomy table (from setaTaxonomyDF documentation)
#' meta <- data.frame(
#'   bc          = paste0("cell", 1:6),
#'   fine_type   = c("AT1","AT2","AT1","Fib1","Fib1","AT2"),
#'   mid_type    = c("Alv","Alv","Alv","Fib","Fib","Alv"),
#'   broad_type  = c("Epi","Epi","Epi","Stroma","Stroma","Epi")
#' )
#' taxDF <- setaTaxonomyDF(meta,
#'              resolution_cols = c("broad_type","mid_type","fine_type"))
#'
#' ## Fake counts (2 samples × n_taxa leaves)
#' set.seed(687)
#' cnt <- matrix(rpois(2 * 3, 10), nrow = 2)
#' colnames(cnt) <- rownames(taxDF)
#'
#' ## (a) One balance: Epi vs Stroma (broad_type level)
#' bal1 <- list(num = "Epi", denom = "Stroma")
#' out1 <- setaBalance(cnt, bal1,
#'                     taxonomyDF = taxDF, taxonomy_col = "broad_type")
#' out1$counts
#'
#' ## (b) Two balances in one call
#' bals <- list(
#'   epi_vs_stroma = list(num = "Epi", denom = "Stroma"),
#'   AT1_vs_AT2    = list(num = "AT1", denom = "AT2")
#' )
#' out2 <- setaBalance(cnt, bals,
#'                     taxonomyDF = taxDF, taxonomy_col = "fine_type")
#' out2$counts
#' @export
setaBalance <- function(counts,
                        balances,
                        taxonomyDF          = NULL,
                        taxonomy_col        = NULL,
                        normalize_to_parent = FALSE,
                        pseudocount         = 1) {
    
    if (!is.matrix(counts))
        stop("'counts' must be a matrix (samples x taxa).")
    if (!is.list(balances) || length(balances) == 0)
        stop("'balances' must be a list.")
    
    ## Allow single unnamed balance ----------------------------------------
    single <- (!is.null(balances$num) && !is.null(balances$denom))
    if (single) balances <- list(balance = balances)
    
    bal_names <- names(balances)
    if (is.null(bal_names) || any(bal_names == ""))
        bal_names <- paste0("Balance", seq_along(balances))
    
    log_counts <- log(counts + pseudocount)
    out        <- matrix(NA_real_, nrow = nrow(counts), ncol = length(balances),
                         dimnames = list(rownames(counts), bal_names))
    
    for (i in seq_along(balances)) {
        
        b <- balances[[i]]
        if (is.null(b$num) || is.null(b$denom))
            stop("Each balance must have 'num' and 'denom' elements.")
        
        num_idx <- resolveGroup(b$num,   counts, taxonomyDF, taxonomy_col)
        den_idx <- resolveGroup(b$denom, counts, taxonomyDF, taxonomy_col)
        
        if (length(intersect(num_idx, den_idx)))
            stop("Numerator and denominator overlap in balance '", bal_names[i], "'")
        
        if (normalize_to_parent) {
            parent_idx <- c(num_idx, den_idx)
            parent_sum <- rowSums(counts[, parent_idx, drop = FALSE])
            log_counts[, parent_idx] <-
                log(counts[, parent_idx, drop = FALSE] / parent_sum + pseudocount)
        }
        
        out[, i] <- rowMeans(log_counts[, num_idx, drop = FALSE]) -
            rowMeans(log_counts[, den_idx, drop = FALSE])
    }
    
    list(method = "balance", counts = out)
}

#' Wrapper for Compositional Transforms with Optional Within-Lineage Resolutions
#' A convenience function that dispatches to one of the transforms:
#' CLR, ALR, ILR, percent, or logCPM. Note that the input \code{counts} matrix
#' should have rows as samples and columns as taxa. Optionally, you can supply
#' a taxonomy data frame to perform a within-lineage transform at a specified
#' resolution.
#'
#' @param counts A numeric matrix with rows as samples and columns as taxa.
#' @param method A character string specifying which transform to apply.
#'     One of \code{"CLR"}, \code{"ALR"}, \code{"ILR"}, \code{"percent"},
#'     \code{"logCPM"} or \code{"balance"}.
#' @param ref Reference taxon (only used if \code{method = "ALR"}). This can be
#'     a taxon name or a column index.
#' @param taxTree Optional tree for ILR (not yet implemented).
#' @param pseudocount Numeric, used by CLR, ALR, ILR, and logCPM. Default is 1.
#' @param size_factors For logCPM scaling. If \code{NULL}, uses row sums.
#' @param taxonomyDF Optional data frame specifying higher-level groupings
#'     for each taxon. Row names of \code{taxonomyDF} should match
#'     \code{colnames(counts)}.
#' @param balances For `"balance"`: a single balance list or a named list;
#' @param normalize_to_parent Logical, passed to [setaBalance()].
#' @param taxonomy_col The column of \code{taxonomyDF} indicating which lineage
#'     each taxon belongs to. Only used if \code{within_resolution = TRUE}.
#' @param within_resolution Logical. If \code{TRUE}, applies the transform
#'     within each lineage of taxa defined by \code{taxonomyDF[[taxonomy_col]]}
#'     separately, then merges them back into the original matrix structure.
#'     Default is \code{FALSE}. Ignored for `"balance"`.
#'
#' @return A list with the following elements:
#' \describe{
#'     \item{transform_method}{The core transform, e.g. \"CLR\", \"ALR\", etc.}
#'     \item{within_resolution}{Logical indicating if a within-lineage transform
#'         was used.}
#'     \item{grouping_var}{The name of the column in \code{taxonomyDF} used for
#'         grouping (lineages) if \code{within_resolution = TRUE}, otherwise
#'         \code{NULL}.}
#'     \item{counts}{The resulting matrix after transformation, with the same
#'         dimensions as the input \code{counts}.}
#' }
#'
#' @examples
#' mat <- matrix(c(1, 2, 4, 8, 3, 6, 9, 12),
#'               nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("TaxonA1", "TaxonA2", "TaxonB1", "TaxonB2")
#'
#' # Build a taxonomy data frame labeling lineages
#' df_lineage <- data.frame(
#'     Lineage = c("LineageA", "LineageA", "LineageB", "LineageB"),
#'     row.names = colnames(mat)
#' )
#'
#' # Apply CLR transform to all columns together
#' out1 <- setaTransform(mat, method = "CLR")
#'
#' # Apply CLR within each Lineage
#' out2 <- setaTransform(
#'     mat,
#'     method = "CLR",
#'     taxonomyDF = df_lineage,
#'     taxonomy_col = "Lineage",
#'     within_resolution = TRUE
#' )
#' @name setaTransform
#' @export
setaTransform <- function(
        counts,
        method         = c("CLR", "ALR", "ILR", "percent", "logCPM", "balance"),
        ref            = NULL,
        taxTree        = NULL,
        pseudocount    = 1,
        size_factors   = NULL,
        taxonomyDF     = NULL,
        taxonomy_col   = NULL,
        within_resolution   = FALSE,
        balances            = NULL,
        normalize_to_parent = FALSE
) {
    method <- match.arg(method)
    
    if (!is.matrix(counts))
        stop("'counts' must be a matrix with samples in rows and taxa in columns.")
    
    ##  Balances require their own block - balances can exist btwn clades, so
    ##  partitioning as below doesn't work
    if (method == "balance") {
        if (is.null(balances))
            stop("For method = 'balance' please supply the 'balances' argument.")
        
        res <- setaBalance(counts,
                           balances            = balances,
                           taxonomyDF          = taxonomyDF,
                           taxonomy_col        = taxonomy_col,
                           normalize_to_parent = normalize_to_parent,
                           pseudocount         = pseudocount)
        
        return(list(
            method            = res$method,
            within_resolution = FALSE,
            grouping_var      = NULL,
            counts            = res$counts
        ))
    }
    
    # All other methods
    if (!within_resolution || is.null(taxonomyDF) || is.null(taxonomy_col)) {
        result <- switch(
            method,
            "CLR"     = setaCLR(counts, pseudocount = pseudocount),
            "ALR"     = setaALR(counts, ref = ref, pseudocount = pseudocount),
            "ILR"     = setaILR(counts, taxTree = taxTree, pseudocount = pseudocount),
            "percent" = setaPercent(counts),
            "logCPM"  = setaLogCPM(counts,
                                   pseudocount  = pseudocount,
                                   size_factors = size_factors)
        )
        return(list(
            method            = result$method,
            within_resolution = FALSE,
            grouping_var      = NULL,
            counts            = result$counts
        ))
    }
    
    ## Reference‑frames
    if (!all(colnames(counts) %in% rownames(taxonomyDF)))
        stop("Some colnames(counts) are not in rownames(taxonomyDF).")
    
    taxonomyDF   <- taxonomyDF[colnames(counts), , drop = FALSE]
    group_vector <- taxonomyDF[[taxonomy_col]]
    unique_groups <- unique(group_vector)
    
    newCounts    <- counts
    final_method <- NULL
    
    for (grp in unique_groups) {
        idx        <- which(group_vector == grp)
        subCounts  <- counts[, idx, drop = FALSE]
        
        result <- switch(
            method,
            "CLR"     = setaCLR(subCounts, pseudocount = pseudocount),
            "ALR"     = setaALR(subCounts, ref = ref, pseudocount = pseudocount),
            "ILR"     = setaILR(subCounts,
                                taxTree = taxTree,
                                pseudocount = pseudocount),
            "percent" = setaPercent(subCounts),
            "logCPM"  = setaLogCPM(subCounts,
                                   pseudocount  = pseudocount,
                                   size_factors = size_factors)
        )
        newCounts[, idx] <- result$counts
        if (is.null(final_method)) final_method <- result$method
    }
    
    list(
        method            = final_method,
        within_resolution = TRUE,
        grouping_var      = taxonomy_col,
        counts            = newCounts
    )
}