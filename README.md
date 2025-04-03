# SETA: Ecological Compositional Analysis of scRNA-seq Data

SETA aims to make compositional analysis user friendly and easy to understand with visualization of single-cell RNA-seq data, allowing researchers to easily assess cell-type proportions and distribution changes across biological conditions.

## Project Status

**Under Construction**

- Based on a non-Bioconductor-compliant package, [SETA](https://github.com/jo-m-lab/SETA)

## Planned Features

- Import functions for `SingleCellExperiment` or `Seurat` objects
- Utilities for compositional analysis: proportionality, compositional analysis, sample-level trajectories
- Modular architecture to integrate readily with standard scRNA-seq toolchains

# To Do List

- Bioconductor Guidelines
  - Organize code into R/, man/, vignettes/
  - Include unit tests and pass R CMD check --as-cran
  - Provide comprehensive documentation and examples
- Versioning
  - Adopt Bioconductor versioning scheme (e.g., x.y.z) 
  - Increment version on each commit to reflect development milestones
- Continuous Integration
  - Configure GitHub Actions or similar for automated checks
  - Validate R CMD check, test coverage, and code style
- Composition Object Creation
  - Implement reference frame for normalization
  - Integrate latent space methods (PCA, PCoA, NMDS, RDA, PLS-DA)
  - Add trajectory capabilities
- Compositional Transforms
  - ALR
  - CLR
  - ILR (with or without balances)
- Methods for Cell Type Trees
  - Build hierarchical or phylogenetic structures from input counts
- Analysis Methods
  - Compare distances in latent or full space (Aitchison, Euclidean)
  - Compare loadings of different latent spaces
  - Compare results of rank-based tests (Wilcoxon/MWU)
  - Build proportionality or Pearson correlation networks of cell type compositions

## Installation

Until acceptance on Bioconductor, you can install from GitHub:

```{r}
install.packages("remotes")
remotes::install_github("CellDiscoveryNetwork/cdnSETA")
```

Contributions are welcome. Please open an issue or pull request with any suggestions or enhancements. 

## License

Released under an MIT open-source license
