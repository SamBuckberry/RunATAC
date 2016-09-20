#' Read BAM file Tn5 insertion positions into GRanges object
#' 
#' @param bam_file Path to a BAM formatted file.
#' @param yieldSize The number of records from the BAM file to read in each chunk. See ?BamFile.
#' @param ... Arguments passed to ScanBamParam.
#' @return A GRanges object of Tn5 insertion positions
#' @examples 
#' \dontrun{
#' tn <- import_atac_pos("test_100k.bam")
#' }
#' @export
read_atac_bam_tn <- function(bam_file, yieldSize=1e6, ...){
        
        # Allow the setting of parameters for Bam file scan such as Which regions
        param <- Rsamtools::ScanBamParam(...)
        
        message("Processing BAM file. This this may take a while for large files...")
        message("If you have lots of RAM, increase yieldSize to speed things up.")
        
        # Loop for reading the BAM file in chunks
        bf <- Rsamtools::BamFile(file = bam_file, yieldSize = yieldSize)
        open(bf)
        gr <- NULL
        repeat {
                chunk <- GenomicAlignments::readGAlignments(bf, param = param)
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
        
        gr <- GRangesList(gr) %>% unlist()
        
        message("Offsetting alignments to Tn5 integration site.")
        
        # Offset the reads to correspond to tn5 insertion site
        pos <- gr[strand(gr) == "+"] %>% 
                GenomicRanges::shift(shift=4) %>%
                GenomicRanges::resize(width = 2, fix = "start")
        
        neg <- gr[strand(gr) == "-"] %>%
                GenomicRanges::shift(shift = -5) %>%
                GenomicRanges::resize(width = 2, fix = "start")
        
        # Return the pos and neg strands together
        gr <- c(pos, neg)
        
        message("Done!")
        
        return(gr)
}
