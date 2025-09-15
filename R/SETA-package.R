#' Single Cell Ecological Taxonomic Analysis
#'
#' @description
#' SETA provides tools for compositional analysis of single-cell RNA-seq data, 
#' applying ecological principles and compositional data analysis (CoDA) methods 
#' to understand cell-type abundance patterns. The package offers a comprehensive 
#' workflow for extracting cell-type counts, applying various compositional 
#' transforms, performing latent space analysis, and visualizing results.
#'
#' @details
#' SETA treats cell-type abundance data as compositional data, similar to how 
#' ecologists analyze species abundance in environmental samples. This approach 
#' is appropriate because cell-type proportions sum to 1 (or 100%), and changes 
#' in one cell type affect all others.
#'
#' The package implements several compositional transforms:
#' \itemize{
#'   \item \strong{CLR (Centered Log-Ratio)}: Centers log-transformed data around the geometric mean
#'   \item \strong{ALR (Additive Log-Ratio)}: Uses a reference cell type as denominator
#'   \item \strong{ILR (Isometric Log-Ratio)}: Projects data onto orthonormal basis
#'   \item \strong{Balance transforms}: User-defined log-ratios between groups of cell types
#' }
#'
#' SETA also provides multi-resolution analysis capabilities, allowing users to 
#' analyze data at different taxonomic levels (e.g., broad cell types vs. specific subtypes).
#'
#' Key functions include:
#' \itemize{
#'   \item \code{\link{setaCounts}}: Extract cell-type count matrices from single-cell objects
#'   \item \code{\link{setaTransform}}: Apply compositional transforms (CLR, ALR, ILR, balance)
#'   \item \code{\link{setaLatent}}: Perform dimensionality reduction (PCA, PCoA, NMDS)
#'   \item \code{\link{setaDistances}}: Calculate compositional distances between samples
#'   \item \code{\link{setaTaxonomyDF}}: Create hierarchical taxonomies for multi-resolution analysis
#'   \item \code{\link{taxonomy_to_tbl_graph}}: Convert taxonomies to graph objects for visualization
#' }
#'
#' For detailed examples, see the package vignettes:
#' \itemize{
#'   \item \code{vignette("introductory_vignette", package = "SETA")}
#'   \item \code{vignette("comparing_samples", package = "SETA")}
#'   \item \code{vignette("reference_frames", package = "SETA")}
#' }
#'
#' The package is designed to be educational, teaching users about compositional 
#' data analysis principles while providing practical tools for single-cell research.
#'
#' @author
#' Kyle Kimler \email{kkimler@broadinstitute.org} (aut, cre)
#' 
#' Marc Elosua-Bayes (aut)
#'
#' @references
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data.
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 
#' 44(2), 139-177.
#'
#' Greenacre, M. (2018). Compositional Data Analysis in Practice. 
#' CRC Press.
#'
#' Pawlowsky-Glahn, V., Egozcue, J. J., & Tolosana-Delgado, R. (2015). 
#' \emph{Modeling and Analysis of Compositional Data}. John Wiley & Sons.
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/CellDiscoveryNetwork/SETA}
#'   \item Report bugs at \url{https://github.com/CellDiscoveryNetwork/SETA/issues}
#' }
#'
#' @keywords package compositional single-cell ecology
#' @name SETA-package
#' @aliases SETA
NULL
