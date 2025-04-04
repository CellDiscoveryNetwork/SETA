setaCLR <- function(counts, pseudocount = 1) {
  # Basic sanity check
  if(!is.matrix(counts)) stop("'counts' must be a matrix.")
  
  # Centered log-ratio transform
  # Pseudocount = 1 to avoid issues with zeros
  counts <- counts + pseudocount
  gm <- exp(rowMeans(log(counts)))
  clr_mat <- log(counts / gm)
  
  # Return a list with the transformed matrix and metadata
  list(
    method = "CLR",
    counts = clr_mat
  )
}

setaILR <- function(counts, boxcox_p = 0, taxTree = NULL, pseudocount = 1) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (!is.null(taxTree)) {
    message("A taxTree was provided but is not yet used. Defaulting to Helmert basis.")
  }
  # Add pseudocount to avoid log(0)
  counts <- counts + pseudocount
  # Take the log
  y <- log(counts)
  # Optional Box-Cox–like step on the log-values if p != 0
  # This is suggested by whuber:
  # https://stats.stackexchange.com/questions/259208/how-to-perform-isometric-log-ratio-transformation
  if (boxcox_p != 0) {
    y <- (exp(boxcox_p * y) - 1) / boxcox_p
  }
  
  # Subtract row means to recenter
  y <- y - rowMeans(y, na.rm = TRUE)
  
  # Build the Helmert basis
  k <- ncol(y)
  if (k < 2) stop("ILR requires at least 2 columns (cell types).")
  
  # contr.helmert(k) yields a k x (k-1) matrix; 
  # we transpose and normalize each row by sqrt((2:k)*(2:k-1)).
  # This ensures an orthonormal basis for ILR rotation.
  H <- stats::contr.helmert(k)          # k x (k-1)
  H <- t(H) / sqrt((2:k)*(2:k-1))       # (k-1) x k
  
  # For interpretability, we rename the resulting columns.
  if (!is.null(colnames(counts)) && k > 1) {
    colnames(y)[-1] <- paste0(colnames(y)[-1], ".ILR")
  }
  
  # Apply the Helmert rotation
  ilr_vals <- y %*% t(H)
  
  # Return a list with metadata
  list(
    method = paste0("ILR_Helmert (boxcox_p=", boxcox_p, ")"),
    counts = ilr_vals
  )
}

setaALR <- function(counts, ref, pseudocount = 1) {
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
  counts <- counts + pseudocount
  refVec <- counts[, refCol]
  alrCounts <- log(counts / refVec)
  alrCounts <- alrCounts[, -refCol, drop = FALSE]
  refName <- if (is.character(ref)) ref else colnames(counts)[ref]
  list(method = paste0("ALR_ref=", refName), counts = alrCounts)
}

setaPercent <- function(counts) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  # Example: compute % by columns (assuming columns = samples)
  col_sums <- colSums(counts)
  pct_counts <- t(t(counts) / col_sums * 100)
  list(method = "percent", counts = pct_counts)
}

setaLogCPM <- function(counts, pseudocount = 1, size_factors = NULL,
                       scale_factor = 1e6) {
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  if (is.null(size_factors)) size_factors <- colSums(counts)
  # Add pseudocount, divide by size factor,
  # then multiply by size factor
  cpm <- sweep(counts + pseudocount, 2, 
               size_factors + pseudocount, FUN = "/") * scale_factor
  log_cpm <- log2(cpm)
  list(method = "logCPM", counts = log_cpm)
}

setaTransform <- function(
  counts, 
  method = c("CLR", "ALR", "ILR", "percent", "logCPM"),
  ref = NULL,
  taxTree = NULL,
  pseudocount = 1,
  size_factors = NULL
) {
  method <- match.arg(method)
  if (!is.matrix(counts)) stop("'counts' must be a matrix.")
  
  switch(method,
    "CLR" = {
      return(setaCLR(counts, pseudocount = pseudocount))
    },
    "ALR" = {
      return(setaALR(counts, ref = ref, pseudocount = pseudocount))
    },
    "ILR" = {
      return(setaILR(counts, taxTree = taxTree, pseudocount = pseudocount))
    },
    "percent" = {
      return(setaPercent(counts))
    },
    "logCPM" = {
      return(setaLogCPM(counts, pseudocount = pseudocount, 
                        size_factors = size_factors))
    }
  )
}