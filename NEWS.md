# SETA NEWS

## Version 0.99.2 (2025-09-15)

### Bioconductor Preparation
- Updated package for Bioconductor review
- Added package-level documentation (`?SETA`)
- Bioconductor installation instructions
- Added proper error handling for optional dependencies
- Created NEWS file for change tracking

### Documentation Improvements
- Added detailed motivation and comparison sections to introductory vignette
- Vignette descriptions with step-by-step explanations
- Added backticks for code formatting throughout vignettes
- Added BiocStyle dependency for proper vignette formatting

### Package Structure
- Created `R/SETA-package.R` with roxygen2 documentation
- Added proper error handling for `SeuratObject` and `SingleCellExperiment` dependencies
- Updated DESCRIPTION with BiocStyle dependency
- package-level help with overview

## Version 0.99.1 (2025-08-26)

### New Features
- Initial submission of SETA package for compositional analysis of single-cell RNA-seq data
- Added compositional transform functions (CLR, ALR, ILR, balance)
- Implemented latent space analysis methods (PCA, PCoA, NMDS)
- Created multi-resolution analysis capabilities with hierarchical taxonomies
- Added distance calculation functions for compositional data
- Included vignettes with educational content

### Functions Added
- `setaCounts()`: Extract cell-type count matrices from single-cell objects
- `setaTransform()`: Apply various compositional transforms
- `setaCLR()`: Centered log-ratio transformation
- `setaALR()`: Additive log-ratio transformation
- `setaILR()`: Isometric log-ratio transformation
- `setaBalance()`: User-defined balance transformations
- `setaPercent()`: Percentage transformation
- `setaLogCPM()`: Log counts-per-million transformation
- `setaLatent()`: Dimensionality reduction methods
- `setaDistances()`: Calculate compositional distances
- `setaTaxonomyDF()`: Create hierarchical taxonomies
- `taxonomy_to_tbl_graph()`: Convert taxonomies to graph objects
- `setaMetadata()`: Extract sample-level metadata
- `resolveGroup()`: Resolve group specifications for balance transforms

### Documentation
- Added package documentation
- Created three detailed vignettes:
  - Introductory vignette with basic workflow
  - Comparing samples vignette with statistical analysis
  - Reference frames vignette with multi-resolution analysis
- Added function-level documentation with examples
- Included mathematical foundations and references

### Data
- Added mock data functions for testing and examples:
  - `mockSC()`: Mock Seurat object
  - `mockSCE()`: Mock SingleCellExperiment object
  - `mockCount()`: Mock count matrix
  - `mockLong()`: Mock long-form data frame

### Dependencies
- Core dependencies: dplyr, MASS, Matrix, stats, tidygraph, rlang
- Suggested dependencies: BiocStyle, caret, corrplot, ggplot2, ggraph, knitr, methods, patchwork, reshape2, rmarkdown, SeuratObject, Seurat, SummarizedExperiment, SingleCellExperiment, TabulaMurisSenisData, tidyr, tidytext, testthat

### Bioconductor Integration
- Prepared package for Bioconductor submission
- Added BiocStyle dependency for vignette formatting
- Updated README with Bioconductor installation instructions
- Created package-level documentation for `?SETA`
- Added NEWS file for change tracking

### Testing
- Test suite covering all major functions
- Tests for error handling and edge cases
- Validation of mathematical correctness for transforms
- Mock data testing for integration with single-cell objects
