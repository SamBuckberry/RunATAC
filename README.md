RunATAC
-------

### An R package for processing and plotting ATAC-seq data

This package is currently is it's infantcy. Don't expect any real functionality for a while...

Plotting ATAC-seq footprints using PWMs
---------------------------------------

In this example, we will plot the Tn5 insertion signal around CTCF motifs in ATAC-seq peaks

Get the CTCF motif

``` r
library(JASPAR2014)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, as.vector, cbind,
    ##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
    ##     grep, grepl, intersect, is.unsorted, lapply, lengths, Map,
    ##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    ##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    ##     setdiff, sort, table, tapply, union, unique, unlist, unsplit

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: IRanges

    ## Loading required package: XVector

``` r
library(TFBSTools)
library(seqLogo)
```

    ## Loading required package: grid

    ## 
    ## Attaching package: 'seqLogo'

    ## The following object is masked from 'package:TFBSTools':
    ## 
    ##     seqLogo

``` r
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
