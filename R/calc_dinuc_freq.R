#' Calculate the dinucleotide insertion frequency for Tn5 insertions.
#' 
#' @param ins_gr A GRanges object. 
#' Typically of ATAC-seq Tn5 insertions generated with the function read_atac_insertions. 
#' @param genome A BSgenome object.
#' @return A data.frame containing dinucleotide counts and precentages.
#' @export

calc_dinuc_freq <- function(ins_gr, genome){
        
        # Get the insertion base and preceeding base for dinucleotide frequency
        di_ins <- GenomicRanges::resize(x = ins, width = 2, fix = "center")
        
        # Get the dinucleotides
        ins_seq <- BSgenome::getSeq(x = genome, di_ins)
        
        # Calculate the dinucleotide frequencies
        ins_freq <- table(ins_seq) %>% data.frame()
        ins_freq$percent <- ins_freq$Freq / sum(ins_freq$Freq)
        
        colnames(ins_freq) <- c("dinucleotide", "count", "percentage")
        return(ins_freq)
}