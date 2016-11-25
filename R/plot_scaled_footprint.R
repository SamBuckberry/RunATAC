#' Plot scaled ATAC-seq footprint signal.
#' @param tn_dist A vector of distances of tn5 insertions relative to motif centre.
#' See atac_motif_dist function.
#' @param window The plotting window in bases. 
#' @param ... Optional arguments passed to plot function.
#' @param scale_to_lib Should plot be scaled to library size?. Logical.
#' @param libsize The ATAC-seq library size for normalising signal to tn5 intergration events per million.
#' @param ylab Character. Y axis label.
#' @param xlab Character. X axis label.
#' @return An ATAC-seq footprint plot in the current graphics device. No values are returned.
#' @export
plot_scaled_footprint <- function(tn_dist, window=200, scale_to_lib=FALSE,
                                  libsize=1, xlab="", ylab="", ...){
        
        # Subset for window
        tn_dist <- tn_dist[abs(tn_dist) <= window/2]
        
        # Histogram
        breaks <- ((0 - (window / 2))-1) : ((0 + (window / 2))+1)
        h <- graphics::hist(tn_dist, plot = FALSE, breaks = breaks, include.lowest = FALSE)
        
        # Trim the first and last histogram bins
        h$counts <- h$counts[2:(length(h$counts)-1)]
        h$mids <- h$mids[2:(length(h$mids)-1)]
        
        # Scale
        if(scale_to_lib == TRUE){
                # Scale based on lib size to get insertions per million
                y_scaled <- (h$counts / (libsize/1e6))
        } else {
                #Scale the histogram from 0-1 to get fraction of signal
                maxs <- max(h$counts)
                mins <- min(h$counts)
                y_scaled <- scale(h$counts, center = mins, scale = maxs - mins)  
        }
        
        # Generate plot
        graphics::plot(h$mids, y_scaled, type="l", ...,
             ylab=ylab,
             xlab=xlab)
}