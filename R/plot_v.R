#' Calculate distances of Tn5 insertions from motif.
#' 
#' @param frags_gr A GRanges object with ATAC-seq fragments.
#' See the read_atac_frags function.
#' @param motif_pos_gr A GRanges object with the positions of the PWM in ATAC-seq peaks.
#' See the motif_gr function.
#' @param x_range Numeric of length 2. The range of data to include on the x axis.
#' @param max_frag Integer of length 1. The maximum fragment size to include in plot. 
#' @param ... arguments passed to ggplot::geom_point
#' @return An object of class 'ggplot'
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @import IRanges
#' @import ggplot2
#' @export
plot_v <- function(frags_gr, motif_pos_gr, x_range=c(-250, 250),
                   max_frag=300, ...){
        
        # Resize fragments on centres
        frag_centre <- GenomicRanges::resize(x = frags_gr, width = 1,
                                             fix = 'center')
        
        # Calculate distance to nearest motif
        nearest_ind <- GenomicRanges::nearest(x = frag_centre,
                                              subject = motif_pos_gr)
        
        nuc_dists <- GenomicRanges::start(frag_centre) - 
                GenomicRanges::start(motif_pos_gr[nearest_ind])         
        
        # Get the fragment widths
        frag_size <- width(frags_gr)
        
        # Filter on size
        df <- data.frame(frag_size=frag_size, distance=nuc_dists)
        df <- df[df$distance %in% x_range[1]:x_range[2], ]
        df <- df[df$frag_size <= max_frag, ]
        df <- as.data.frame(df)
        
        
        # plot the data
        gg <- ggplot2::ggplot(df, aes(x = distance, y = frag_size)) +        
                geom_point(...) +
                #stat_binhex() +
                scale_x_continuous(expand = c(0, 0)) + 
                scale_y_continuous(expand = c(0, 0)) +
                xlab("Relative distance (bp)") +
                ylab("Fragment size (bp)") +
                theme_bw() +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
        return(gg)
        
}