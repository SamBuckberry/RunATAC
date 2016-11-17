
plot_v <- function(frags_gr, motif_pos_gr, range=200){

        # Resize ranges to 1 for distance calculation
        frag_centres <- resize(frags_gr, width = 1, fix = 'center')        
        motif_centres <- resize(motif_pos_gr, width = 1, fix = 'center')
        
        frag_nearest <- nearest(frag_centres, subject = motif_centres)
        frag_dist <- frag_dist@elementMetadata$distance
        
                
        stopifnot(length(frag_dist) == length(frag_centres))
        
        frag_widths <- width(frags_gr)
        
        smoothScatter(x = frag_dist, y = frag_widths)                        
}