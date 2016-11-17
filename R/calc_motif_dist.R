#' Calculate distances of Tn5 insertions from motif.
#' 
#' @param tn_signal_gr A GRanges object with tn5 insertions. See the read_atac_bam_tn function.
#' @param motif_pos_gr A GRanges object with the positions of the PWM in ATAC-seq peaks.
#' @param range The distance in bases flanking the the motif to calculate signal.
#' @return A numeric vector of distances of Tn5 insertion sites to motif centre.
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @import IRanges
#' @export
calc_motif_dist <- function(tn_signal_gr, motif_pos_gr, range=200){
        
        # Resize to range of signal collection
        motif_pos_gr <- GenomicRanges::resize(motif_pos_gr, width = range*2,
                                              fix = 'center')
        
        # Subset the Tn5 signal
        tn_signal_gr <- tn_signal_gr[IRanges::overlapsAny(query = tn_signal_gr,
                                                          subject = motif_pos_gr)]
        
        #For each motif, get the distances of nearby Tn5 insertions
        get_tn_dist <- function(x){
                
                # Subset to one region
                motif_pos <- motif_pos_gr[x]
                
                # Get the tn5 hits for range
                tn_sub <- tn_signal_gr[IRanges::overlapsAny(query = tn_signal_gr,
                                                            subject = motif_pos)]
                
                # Calculate the relative distance to the motif centre for insertions
                tn_dist <- (start(motif_pos) + range) - start(tn_sub)
                
                # Check if no Tn5 insertion events at position
                if(length(tn_dist) < 1){
                        tn_dist <- NA
                }
                
                return(tn_dist)
        }
        
        # Apply the function over all motif positions
        dists <- lapply(1:length(motif_pos_gr), get_tn_dist)
        dists <- dists[!is.na(dists)] %>% unlist()
        
        # Return distances
        return(dists)
}
