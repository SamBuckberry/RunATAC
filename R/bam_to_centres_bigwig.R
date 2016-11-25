#' Write a bigwig file of Tn5 insertion depth
#' 
#' @param bam_file Path to an ATAC-seq paired-end BAM file.
#' @param file Path for output bigwig file. Extension should be .bw or .bigwig.
#' @param max_insert The Maximum insert size for read fragments.
#' @param nuc_range The nucleosome spaning fragment size range.
#' @param scale_cpm Logical. Should the output be scaled to number of 
#' million mapped reads?
#' @param ... Arguments for readGAlignmentPairs.
#' @return No object is returned. A file is written. 
#' @export
bam_to_centres_bw <- function(bam_file, file, max_insert=2000,
                              nuc_range=c(189, 247), scale_cpm=FALSE, ...){

        # Check inputs
        if (class(max_insert) != "numeric") {
                stop("max_insert_size is not numeric!")
        }
        
        message("Processing BAM file. This this may take a while for large files...")

        # Read the bam file pairs to GRanges object
        gr <- GenomicAlignments::readGAlignmentPairs(file = bam_file, ...)
        gr <- GenomicRanges::GRanges(gr)
        
        
        message("Filtering for nucleosome spanning fragments")
        # Remove pairs with inserts larger than 2000 bases
        gr <- gr[width(gr) <= max_insert]
        
        # Remove inserts not in nucleosome range
        gr <- gr[width(gr) >= min(nuc_range) & width(gr) <= max(nuc_range)]
        
        
        # Reduce to fragment centre for nucleosome position
        gr <- GenomicRanges::resize(gr, width = 1, fix = 'center')
        
        
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