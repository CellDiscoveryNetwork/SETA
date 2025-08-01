% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{setaTaxonomyDF}
\alias{setaTaxonomyDF}
\title{Build a taxonomy data frame at multiple resolutions}
\usage{
setaTaxonomyDF(
  obj,
  resolution_cols = c("fine_type", "mid_type", "broad_type"),
  bc_col = "bc"
)
}
\arguments{
\item{obj}{A data.frame or similar object containing cell metadata.}

\item{resolution_cols}{A character vector of column names
indicating hierarchical taxonomy (from broad to fine).}

\item{bc_col}{Optional. The name of the column containing barcodes,
or "rownames" if they are row names.}
}
\value{
A `data.frame` with one row per unique value of the finest label
        (the last entry in `resolution_cols`), and one column for each
        resolution level. The row names are set to the finest label values.
        If any finest label maps to more than one combination of coarser
        labels, the function stops with an informative error.
}
\description{
setaTaxonomyDF() converts **one long-form metadata data.frame**-typically
colData(sce), seu@meta.data, or any frame you already have, into a tidy
taxonomy table.  Each row corresponds to a unique value of the *finest*
label (the **last** element of `resolution_cols`), and every coarser label
sits in its own column.
}
\details{
## What the input must contain
* exactly **one row per cell**
* at least one **barcode** column (default `"bc"`).
  Pass `bc_col = "rownames"` if barcodes live in `rownames(obj)`.
* **all** columns listed in `resolution_cols`

No `Seurat`/`SingleCellExperiment` objects are accepted here: extract their
metadata/colData first, then hand it in as a `data.frame`

## Value
A `data.frame` whose **rownames** are the finest label. If any finest label
maps to more than one set of coarser labels the function should stop with an
informative error.
}
\examples{
meta <- data.frame(
  bc          = paste0("cell", 1:6),
  fine_type   = c("AT1","AT2","AT1","Fib1","Fib1","AT2"),
  mid_type    = c("Alv","Alv","Alv","Fib","Fib","Alv"),
  broad_type  = c("Epi","Epi","Epi","Stroma","Stroma","Epi")
)
setaTaxonomyDF(meta,
               resolution_cols = c("broad_type","mid_type","fine_type"))

## barcodes can be in rownames with bc_col = "rownames" (as in Seurat Object)
rownames(meta) <- meta$bc
meta$bc <- NULL
setaTaxonomyDF(meta,
               resolution_cols = c("broad_type","mid_type","fine_type"),
               bc_col = "rownames")
}
