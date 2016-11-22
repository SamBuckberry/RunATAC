#' Read BAM file aligned fragments into GRanges object
#' 
#' @param bam_file Path to a BAM formatted file.
#' @param max_insert Maximum insert size allowed. Default is reccomended.
#' @return A GRanges object of aligned fragments
#' @export
read_atac_frags <- function(bam_file, max_insert=2000, mapq=20, yieldSize=1e6, ...){
        
        # Check inputs
        if (class(max_insert) != "numeric") {
                stop("max_insert_size is not numeric!")
        }
        
        
        # Allow the setting of parameters for Bam file scan such as Which regions
        param <- Rsamtools::ScanBamParam(mapqFilter = mapq, ...)
        
        message("Processing BAM file. This this may take a while for large files...")
        message("If you have lots of RAM, increase yieldSize to speed things up.")
        
        # Loop for reading the BAM file in chunks
        bf <- Rsamtools::BamFile(file = bam_file, yieldSize = yieldSize)
        open(bf)
        gr <- NULL
        repeat {
                chunk <- GenomicAlignments::readGAlignmentPairs(bf)
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


