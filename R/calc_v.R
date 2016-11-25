#' Calculate distances of Tn5 insertions from motif.
#' 
#' @param frags_gr A GRanges object with ATAC-seq fragments.
#' See the read_atac_frags function.
#' @param motif_pos_gr A GRanges object with the positions of the PWM in ATAC-seq peaks.
#' See the motif_gr function.
#' @param flank Numeric of length 1. Number of bases flanking motif to plot.
#' @param max_frag Integer of length 1. The maximum fragment size to include in plot. 
#' @return A data.frame with 2 columns.
#' frag_size column contains the ATAC-seq fragment size.
#' The distance column contains the distance from the fragment centure to the centre of
#' the closest motif in motif_pos_gr. 
#' @import GenomicRanges
#' @export
calc_v <- function(frags_gr, motif_pos_gr, flank=250,
                   max_frag=300){
        
        # Resize fragments on centres
        frag_centre <- GenomicRanges::resize(x = frags_gr, width = 1,
                                             fix = 'center')
        strand(frag_centre) <- "*"
        
        # Reduce to fragments within range
        dat_range <- GenomicRanges::resize(x = motif_pos_gr, width = (flank+10)*2,
                                           fix = 'center')
        
        frag_hits <- IRanges::overlapsAny(query = frag_centre, subject = dat_range)
        frag_centre <- frag_centre[frag_hits]
        
        # Calculate distance to nearest motif
        nearest_ind <- GenomicRanges::nearest(x = frag_centre,
                                              subject = motif_pos_gr)
        
        nuc_dists <- GenomicRanges::start(frag_centre) - 
                GenomicRanges::start(motif_pos_gr[nearest_ind])         
        
        # Get the fragment widths
        frag_size <- GenomicRanges::width(frags_gr[frag_hits])
        
        stopifnot(length(frag_size) == length(nuc_dists))
        
        # Filter on size
        df <- data.frame(frag_size=frag_size, distance=nuc_dists)
        df <- df[df$distance %in% -flank:flank, ]
        df <- df[df$frag_size <= max_frag, ]
        df <- as.data.frame(df)
        
        return(df)
}