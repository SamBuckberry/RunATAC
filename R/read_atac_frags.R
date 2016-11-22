#' Read BAM file aligned fragments into GRanges object
#' 
#' @param bam_file Path to a BAM formatted file.
#' @param max_insert Maximum insert size allowed. Default is reccomended.
#' @param yieldSize The number of reads retreived from the BAM file in each chunk.
#' See ?BamFile in Rsamtools.
#' @param ... Arguments passed to ScanBamParam. 
#' @return A GRanges object of aligned fragments
#' @import Rsamtools
#' @import GenomicAlignments
#' @import GenomicRanges
#' @export
read_atac_frags <- function(bam_file, max_insert=2000, yieldSize=1e6, ...){
        
        # Check inputs
        if (class(bam_file) != "character" | length(bam_file) !=1) {
                stop("bam_file is not a character of length 1!")
        }
        
        if (class(max_insert) != "numeric" | length(max_insert) !=1) {
                stop("max_insert_size is not a numeric of length 1!")
        }
        
        if (class(yieldSize) != "numeric" | length(yieldSize) != 1) {
                stop("yieldSize is not a numeric of length 1!")
        }
        
        # Allow the setting of parameters for Bam file scan such as Which regions
        param <- Rsamtools::ScanBamParam(flag = scanBamFlag(isProperPair=TRUE), ...)
        
        message("Processing BAM file. This this may take a while for large files...")
        message("If you have lots of RAM, increase yieldSize to speed things up.")
        
        # Loop for reading the BAM file in chunks
        bf <- Rsamtools::BamFile(file = bam_file, yieldSize = yieldSize)
        open(bf)
        gr <- NULL
        repeat {
                chunk <- GenomicAlignments::readGAlignmentPairs(bf, param=param)
                if (length(chunk) == 0L)
                        break
                chunk_gr <- GenomicRanges::GRanges(chunk)
                if (is.null(gr)) {
                        gr <- chunk_gr
                } else {
                        gr <- c(gr, chunk_gr)
                }
        }
        close(bf)
        
        # Remove pairs with inserts larger than 2000 bases
        gr <- gr[width(gr) <= max_insert]
        
        return(gr)
}


