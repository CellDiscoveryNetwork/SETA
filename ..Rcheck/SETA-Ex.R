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
### Title: Additive Log-Ratio Transform
### Aliases: setaALR

### ** Examples

mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
colnames(mat) <- c("A", "B")
out <- setaALR(mat, ref="A", pseudocount=0)
out$counts




cleanEx()
nameEx("setaCLR")
### * setaCLR

flush(stderr()); flush(stdout())

### Name: setaCLR
### Title: Centered Log-Ratio Transform
### Aliases: setaCLR

### ** Examples

mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
out <- setaCLR(mat, pseudocount=0)
out$counts




cleanEx()
nameEx("setaCounts")
### * setaCounts

flush(stderr()); flush(stdout())

### Name: setaCounts
### Title: Extract Taxonomic Counts from Various Single Cell Objects
### Aliases: setaCounts

### ** Examples





cleanEx()
nameEx("setaILR")
### * setaILR

flush(stderr()); flush(stdout())

### Name: setaILR
### Title: Isometric Log-Ratio Transform
### Aliases: setaILR

### ** Examples

mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
out <- setaILR(mat, boxcox_p=0)
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

mat <- matrix(1:4, nrow=2)
out <- setaLogCPM(mat)
out$counts




cleanEx()
nameEx("setaPercent")
### * setaPercent

flush(stderr()); flush(stdout())

### Name: setaPercent
### Title: Percentage Transform
### Aliases: setaPercent

### ** Examples

mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
out <- setaPercent(mat)
out$counts




cleanEx()
nameEx("setaTransform")
### * setaTransform

flush(stderr()); flush(stdout())

### Name: setaTransform
### Title: Wrapper for Common Compositional Transforms
### Aliases: setaTransform

### ** Examples

mat <- matrix(c(1,2,4,8), nrow=2, byrow=TRUE)
setaTransform(mat, method="CLR")
setaTransform(mat, method="percent")




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
