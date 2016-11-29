#' Calculate the PWM for Tn5 insertions.
#' 
#' @param ins_gr A GRanges object. 
#' Typically of ATAC-seq Tn5 insertions generated with the function read_atac_insertions. 
#' @param genome A BSgenome object.
#' @param ... Parameters for the PWM function in the Biostrings package.
#' @return A PWM matirix containing dinucleotide counts and precentages.
#' @export
calc_ins_pwm <- function(ins_gr, genome, ...){
        
        # Set width of 21bp 
        seq_width=21
        
        # Remove out-of-bounds ranges
        genome_ranges <- methods::as(GenomeInfoDb::seqinfo(genome), "GRanges")
        start(genome_ranges) <- start(genome_ranges) + seq_width
        end(genome_ranges) <- end(genome_ranges) - seq_width
        
        olap_within <- IRanges::overlapsAny(query = ins_gr,
                                            subject = genome_ranges,
                                            type = "within")
        
        # Make the 10 base sequence centred on insertion site
        ins_flank <- GenomicRanges::resize(x = ins_gr[olap_within],
                                           width = seq_width,
                                           fix = 'center')
        
        # Get the sequence for the ranges
        message("Retreiving DNA sequnence for flanking regions...")
        message("This may take a while for large datasets")
        flank_seq <- BSgenome::getSeq(x = genome, ins_flank)
        
        # Calculate the PWM
        message("Calculating the PWM...")
        #ins_pfm <- Biostrings::consensusMatrix(flank_seq)
        ins_pwm <- Biostrings::PWM(flank_seq, ...)
        
        return(ins_pwm)
}