#' Get the positions of the PWM in ATAC-seq peaks.
#' 
#' @param gr A GRanges object. Usually ranges of ATAC-seq peaks. 
#' @param pwm A positional weight matrix of class Matirx.
#' @param genome A BSgenome object. 
#' @param min.score Minimum alignment score for pwm.
#' @return A GRanges object with positions of PWM match in peaks.
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @import IRanges
#' @import BSgenome
#' @export
motif_gr <- function(gr, pwm, genome, min.score="80%"){
        
        if(class(gr) != "GRanges"){
                stop("Object gr is not of class GRanges.")
        }
        if(class(pwm) != "matrix"){
                stop("Object pwm is not of class matrix.")
        }
        if(class(genome) != "BSgenome"){
                stop("Object genome is not of class BSgenome.")
        }
        if(class(min.score) != "character"){
                stop("Object genome is not of class character.")
        }
        
        # Get the nucleotide sequences
        message("Retrieving sequences...")
        sequences <- BSgenome::getSeq(x = genome, names=gr)
        
        # Function to find motif in one sequence
        find_motif_start <- function(x)
        {
                motif_starts_pos <- Biostrings::matchPWM(pwm = pwm, subject = sequences[[x]],
                                                         min.score = min.score) %>% start()
                
                
                motif_starts_neg <- Biostrings::matchPWM(pwm = Biostrings::reverseComplement(pwm),
                                                         subject = sequences[[x]],
                                                         min.score = min.score) %>% start()
                
                starts_pos <- start(gr[x]) + motif_starts_pos
                starts_neg <- start(gr[x]) + motif_starts_neg
                starts <- c(starts_pos, starts_neg)
                
                if (length(starts) == 0){
                        out <- NULL
                } else {
                        ends <- starts + ncol(pwm)
                        out <- GenomicRanges::GRanges(seqnames = seqnames(gr[x]),
                                                      ranges = IRanges(start = starts,
                                                                       end = ends))
                }
                
                return(out)
        }
        
        # Apply the function across all sequences
        message("Searching for motif matches...")
        motif_ranges <- lapply(X = 1:length(sequences), FUN = find_motif_start)
        
        # Remove the NULLs from where no motif was detected
        is.NullOb <- function(x) is.null(motif_ranges[x]) | all(sapply(motif_ranges[x], is.null))
        no_keep <- lapply(1:length(motif_ranges), is.NullOb) %>% unlist() 
        
        # Convert to Granges object
        motif_ranges <- motif_ranges[!no_keep] %>% unlist() %>% GRangesList() %>% unlist()
        
        # Get the PWM alignment scores and add to GRanges
        match_seqs <- BSgenome::getSeq(x = genome, names=motif_ranges)
        
        get_scores <- function(x){
                score_pos <- PWMscoreStartingAt(pwm = pwm,
                                                subject = as.character(match_seqs[x]))
                score_neg <- PWMscoreStartingAt(pwm = pwm,
                                                subject = as.character(reverseComplement(
                                                        match_seqs[x])))
                score <- max(score_pos, score_neg)
                return(score)
        }
        
        message("Calculating PWM match scores...")
        match_scores <- lapply(1:length(match_seqs), get_scores) %>%
                unlist()
        

        motif_ranges$score <- match_scores
        
        
        # Decrease the width to motif centre
        motif_ranges <- resize(motif_ranges, width = 2, fix = 'center')
        
        # Return GRanges object
        return(motif_ranges)
}