#' Calculate the bias scores for the Tn5 PWM.
#' 
#' @param ins_pwm A positional weight matrix generated with the function calc_ins_pwm.
#' Must be a matrix of 21 columns and 4 rows.
#' @param gr A genomic ranges object of regions to obtain signal. All ranges should be the same width.
#' @param genome A BSgenome object.
#' @param range The distance from the centre of gr to flank. Default=100.
#' site for calculating the PWM. 
#' @return A matrix of per-base bias scores, with each row corresponding
#' to the centre of in 'gr' flanked up and downstream by the number of bases specified in 'range'.
#' @export
calc_profile_bias <- function(ins_pwm, gr, genome, range=100){
        
        # Check inputs
        if (ncol(ins_pwm) != 21) {
                stop("PWM does not have 21 columns")
        }
        
        if (nrow(ins_pwm) != 4) {
                stop("PWM does not have 4 rows")
        }
        
        region_gr <- GenomicRanges::resize(gr, width = range+(range+1), 
                                           fix = 'center')
        
        # Pad the end of the ranges with sequence for PWM scores 
        regions_pad <- region_gr
        strand(regions_pad) <- "*"
        end(regions_pad) <- end(regions_pad) + 10
        start(regions_pad) <- start(regions_pad) - 10
        
        # Get the sequences for the regions of interest
        message("Retreiving the sequences...")
        region_seq <-  BSgenome::getSeq(x = genome, regions_pad)
        
        
        # Calculate the PWM match along each sequence
        calc_seq_pwm_score <- function(x){
                
                dna <- region_seq[[x]]
                pwm_score <- Biostrings::PWMscoreStartingAt(pwm = ins_pwm,
                                                            subject = dna,
                                                            starting.at = 1:(length(dna)-20))
        }
        
        message("Calculating insertion bias scores...")
        scores <- lapply(1:length(region_seq), FUN = calc_seq_pwm_score)
        
        # Test if vectors of scores are of equal length to width of GRanges
        stopifnot(all(width(region_gr) == lengths(scores)))
        
        # Get the scores into a data.frame
        scores <- do.call(rbind, scores) %>% data.frame()
        colnames(scores) <- c(-range:range)
        
        return(scores)
}
