library(RunATAC)
library(testthat)
context("GRanges operation are correct")

# Single position ranges
gr1 <- GRanges(seqnames = "chr1",
               ranges = IRanges(start = c(1:10), end = c(1:10)))

# >1 range
gr2 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 10))
gr3 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 6, end = 10))
gr4 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 11, end = 20))

test_that("Number of overlapping ranges calculation is correct", {
        expect_equal(calc_frip(gr1, gr2), 1)
        expect_equal(calc_frip(gr1, gr3), 0.5)
        expect_equal(calc_frip(gr1, gr4), 0)
})

