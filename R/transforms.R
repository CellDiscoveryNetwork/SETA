setaCLR <- function(counts) {
  # Basic sanity check
  if(!is.matrix(counts)) stop("'counts' must be a matrix.")
  
  # Centered log-ratio transform
  # Pseudocount = 1 to avoid issues with zeros
  counts[counts == 0] <- 1
  gm <- exp(rowMeans(log(counts))) 
  clr_mat <- log(counts / gm)
  
  # Return a list with the transformed matrix and metadata
  list(
    method = "CLR",
    counts = clr_mat
  )
}

setaILR <- function(counts, taxTree = NULL) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (!is.null(taxTree)) {
    message("A taxonomic tree was provided but 
             we haven't written the code for this yet.
             Defaulting to Helmert basis. (see CaCoA for more info)")
  }
  counts[counts == 0] <- 1
  gm <- exp(rowMeans(log(counts)))
  clr_counts <- log(counts / gm)
  n <- ncol(clr_counts)
  if (n < 2) stop("ILR requires at least 2 columns (cell types).")
  H <- stats::contr.helmert(n)
  for (j in seq_len(ncol(H))) {
    H[, j] <- H[, j] / sqrt(sum(H[, j]^2))
  }
  ilr_counts <- t(apply(clr_counts, 1, function(x) t(H) %*% x))
  list(method = "ILR_Helmert", counts = ilr_counts)
}

setaALR <- function(counts, ref) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (missing(ref)) stop("Please specify a reference cell type or column.")
  if (is.character(ref)) {
    if (!(ref %in% colnames(counts))) {
      stop("Reference celltype not found in column names.")
    }
    refCol <- which(colnames(counts) == ref)
  } else if (is.numeric(ref)) {
    if (ref < 1 || ref > ncol(counts)) {
        stop("Reference celltype index out of range.")
    }
    refCol <- ref
  } else {
    stop("'ref' must be a character or numeric.")
  }
  counts[counts == 0] <- 1
  refVec <- counts[, refCol]
  alrCounts <- log(counts / refVec)
  alrCounts <- alrCounts[, -refCol, drop = FALSE]
  refName <- if (is.character(ref)) ref else colnames(counts)[ref]
  list(method = paste0("ALR_ref=", refName), counts = alrCounts)
}