#' Calculate the PWM for Tn5 insertions.
#' 
#' @param ins_pwm A positional weight matrix generated with the function calc_ins_pwm.
#' Must be a matrix of 21 columns and 4 rows.
#' @param regions_gr A GRanges object of regions to calculate insertion bias scores. 
#' @param genome A BSgenome object.
#' site for calculating the PWM. 
#' @return A vector of per-base bias scores for each range in regions_gr.
#' @export
calc_ins_bias <- function(ins_pwm, regions_gr, genome){
        
        # Check inputs
        if (ncol(ins_pwm) != 21) {
                stop("PWM does not have 21 columns.")
        }
        
        if (nrow(ins_pwm) != 4) {
                stop("PWM does not have 4 rows.")
        }
        
        # Pad the end of the ranges with sequence for PWM scores 
        regions_pad <- regions_gr
        strand(regions_pad) <- "*"
        end(regions_pad) <- end(regions_pad) + 10
        start(regions_pad) <- start(regions_pad) - 10
        
        # Get the sequences for the regions of interest
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
        stopifnot(all(width(regions_gr) == lengths(scores)))
        
        return(scores)
}


# calc_ins_bias_cross_cor <- function(ins_pwm, regions_gr, genome){
#         
#         # Check inputs
#         if (ncol(ins_pwm) != 21) {
#                 stop("PWM does not have 21 columns.")
#         }
#         
#         if (nrow(ins_pwm) != 4) {
#                 stop("PWM does not have 4 rows.")
#         }
#         
#         # Pad the end of the ranges with sequence for PWM scores 
#         regions_pad <- regions_gr
#         strand(regions_pad) <- "*"
#         end(regions_pad) <- end(regions_pad) + 10
#         start(regions_pad) <- start(regions_pad) - 10
#         
#         # Get the sequences for the regions of interest
#         region_seq <-  BSgenome::getSeq(x = genome, regions_pad)
#         
#         
#         # Calculate the PWM match along each sequence
#         calc_seq_pwm_score <- function(x){
#                 
#                 dna <- region_seq[[x]]
#                 
#                 sub_seq <- dna[1:ncol(ins_pwm)]
#                 sub_seq <- DNAStringSet(sub_seq)
#                 sub_seq <- PWM(sub_seq)
#         
#                 
#                 x_ts <- as.ts(t(ins_pwm))
#                 y_ts <- as.ts(t(sub_seq))
#                 
#                 pp <- pacf(x = x_ts, y = y_ts, type = 'correlation')
#                         
#         }
#         
#         message("Calculating insertion bias scores...")
#         scores <- lapply(1:length(region_seq), FUN = calc_seq_pwm_score)
#         
#         # Test if vectors of scores are of equal length to width of GRanges
#         stopifnot(all(width(regions_gr) == lengths(scores)))
#         
#         return(scores)
# }



