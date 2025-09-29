<img src="man/figures/SETAsmall.jpg?raw=true" align="right" width=250px>  

![Downloads](https://img.shields.io/github/downloads/CellDiscoveryNetwork/SETA/total)
![GitHub stars](https://img.shields.io/github/stars/CellDiscoveryNetwork/SETA?style=social)
![R-CMD-check](https://github.com/CellDiscoveryNetwork/SETA/workflows/R-CMD-check/badge.svg)

## SETA: Ecological Compositional Analysis of scRNA-seq Data

SETA aims to make compositional analysis user friendly and easy to understand by breaking down a full analysis into educational pieces. SETA includes vignettes for visualization of compositional analysis of single-cell RNA-seq data, aiming to make it easier to perform sample-level unsupervised analysis of single-cell data via sample embeddings and distances. 

## Project Status

- Based on a non-Bioconductor-compliant package, [SETA](https://github.com/jo-m-lab/SETA)

## Planned Features

- proportionality networks
- sample-level trajectories
- vegan ecological latent space methods and metrics (like unifrac)

# To Do List

- Compositional Space Calculation
  - Latent space methods (RDA, PLS-DA, tensors!) - vegan and otherwise
  - Add trajectory capabilities
- Compositional Transforms
  - ILR with balances
  - ideas welcome
- Methods for Cell Type Trees
  - Addition of metadata to tree objects
- Analysis Methods
  - Build proportionality or Pearson correlation networks of cell type compositions
  - Tensors and complex modeling
- Vignettes
  - Proportionality Networks
  - Multi-view tensor sample-level analysis

## Installation

### From Bioconductor (recommended)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SETA")
```

### From GitHub (development version)

```r
install.packages("remotes")
remotes::install_github("CellDiscoveryNetwork/SETA")
```

Contributions are welcome. Please open an issue or pull request with any suggestions or enhancements. 

## License

Released under an MIT open-source license
