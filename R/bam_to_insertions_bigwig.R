#' Write a bigwig file of Tn5 insertion depth
#' 
#' @param bam_file Path to an ATAC-seq paired-end BAM file.
#' @param file Path for output bigwig file. Extension should be .bw or .bigwig.
#' @param yieldSize Number of records to yield each time the file is read from with scanBam.
#' See ?BamFile
#' @param scale_cpm Logical. Should the output be scaled to number of 
#' million mapped reads?
#' @return No object is returned. A file is written. 
#' @export
bam_to_insertions_bw <- function(bam_file, file, yieldSize=1e6,
                                scale_cpm=FALSE, ...){
        
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
                GenomicRanges::resize(width = 1, fix = "start")
        
        neg <- gr[strand(gr) == "-"] %>%
                GenomicRanges::shift(shift = -5) %>%
                GenomicRanges::resize(width = 1, fix = "start")
        
        # Return the pos and neg strands together
        gr <- c(pos, neg)
        
        message("Calculating coverage")

        # Calculate Tn insertion coverage
        cov <- IRanges::coverage(gr)
        
        # Scale counts if specified
        if (scale_cpm == TRUE){
                cov <- (cov / length(gr)) * 1e+06
        }
        
        message("writing bigiwg file")
        
        # write Tn5 insertion bigwig
        rtracklayer::export.bw(object = cov, con = file)
        
        message("done!")
} 
