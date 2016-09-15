RunATAC
-------

An R package for processing and plotting ATAC-seq data
------------------------------------------------------

This package is currently is it's infantcy. Don't expect any real functionality for a while...

Plotting ATAC-seq footprints using PWMs
---------------------------------------

In this example, we will plot the Tn5 insertion signal around CTCF motifs in ATAC-seq peaks

Get the CTCF motif

``` r
library(JASPAR2014)
library(TFBSTools)
library(seqLogo)
# Get the motif PWM
### Get the PWM data from the database and reformat for test functions
opts <- list()
opts[["all_versions"]] <- TRUE
opts[["collection"]] <- "CORE"
opts[["matrixtype"]] <- "PFM"
opts[["tax_group"]] <- "vertebrates"
pwm_data <- getMatrixSet(JASPAR2014, opts)
pwm_matrix_list <- Matrix(pwm_data)

# Get mofif PWM for CTCF
ctcf <- pwm_matrix_list$MA0139.1
rownames(ctcf) <- c("A", "C", "G", "T")
```
