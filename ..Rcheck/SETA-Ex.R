pkgname <- "SETA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SETA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("setaALR")
### * setaALR

flush(stderr()); flush(stdout())

### Name: setaALR
### Title: Additive Log-Ratio (ALR) Transform
### Aliases: setaALR

### ** Examples

# Example with 2 samples and 2 taxa:
mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
colnames(mat) <- c("TaxonA", "TaxonB")
# Using TaxonA as the reference.
out <- setaALR(mat, ref = "TaxonA", pseudocount = 0)
out$counts




cleanEx()
nameEx("setaCLR")
### * setaCLR

flush(stderr()); flush(stdout())

### Name: setaCLR
### Title: Centered Log-Ratio (CLR) Transform
### Aliases: setaCLR

### ** Examples

# Example matrix with 2 samples and 2 taxa:
mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
colnames(mat) <- c("Taxon1", "Taxon2")
out <- setaCLR(mat, pseudocount = 0)
out$counts




cleanEx()
nameEx("setaCounts")
### * setaCounts

flush(stderr()); flush(stdout())

### Name: setaCounts
### Title: Extract Taxonomic Counts from Various Single Cell Objects
### Aliases: setaCounts

### ** Examples


# For a data.frame with custom column names:
df <- data.frame(
  barcode = paste0("cell", 1:10),
  cellType = sample(c("Tcell", "Bcell"), 10, TRUE),
  sampleID = sample(c("sample1","sample2"), 10, TRUE)
)
cmat <- setaCounts(df, cell_type_col = "cellType", sample_col = "sampleID", bc_col = "barcode")
print(cmat)





cleanEx()
nameEx("setaILR")
### * setaILR

flush(stderr()); flush(stdout())

### Name: setaILR
### Title: Isometric Log-Ratio (ILR) Transform
### Aliases: setaILR

### ** Examples

# Example matrix: rows are samples, columns are cell types.
mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
colnames(mat) <- c("A", "B")
# ILR transformation reduces the dimension by 1.
out <- setaILR(mat, boxcox_p = 0, pseudocount = 1)
out$counts




cleanEx()
nameEx("setaLatent")
### * setaLatent

flush(stderr()); flush(stdout())

### Name: setaLatent
### Title: Compute a Latent Space from Transformed Counts
### Aliases: setaLatent

### ** Examples

mat <- matrix(rpois(20, lambda=5), nrow=4)  # small 4x5 matrix
colnames(mat) <- paste0("C", 1:5)
clr_out <- setaCLR(mat)
latent_pca <- setaLatent(clr_out, method="PCA", dims=2)
latent_pca$latentSpace




cleanEx()
nameEx("setaLogCPM")
### * setaLogCPM

flush(stderr()); flush(stdout())

### Name: setaLogCPM
### Title: log2(CPM) Transform
### Aliases: setaLogCPM

### ** Examples

mat <- matrix(c(10, 20, 100, 200), nrow = 2, byrow = TRUE)
out <- setaLogCPM(mat, pseudocount = 1)
out$counts




cleanEx()
nameEx("setaMetadata")
### * setaMetadata

flush(stderr()); flush(stdout())

### Name: setaMetadata
### Title: Extract Sample-Level Metadata from Various Objects
### Aliases: setaMetadata

### ** Examples

# Using a Seurat object



cleanEx()
nameEx("setaPercent")
### * setaPercent

flush(stderr()); flush(stdout())

### Name: setaPercent
### Title: Percentage Transform
### Aliases: setaPercent

### ** Examples

mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
out <- setaPercent(mat)
out$counts




cleanEx()
nameEx("setaTransform")
### * setaTransform

flush(stderr()); flush(stdout())

### Name: setaTransform
### Title: Wrapper for Compositional Transforms
### Aliases: setaTransform

### ** Examples

mat <- matrix(c(1,2,4,8), nrow = 2, byrow = TRUE)
# Apply CLR transform:
setaTransform(mat, method = "CLR", pseudocount = 1)
# Apply percent transform:
setaTransform(mat, method = "percent")




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
