#' Plot insert size historgram
#' 
#' @param frag_gr A GRanges object of ATAC-seq fragments
#' @param ... Arguments for graphics::hist
#' @return A plot to the current graphics device
#' @export
plot_insert_size <- function(frag_gr, ...){
        
        # Check inputs
        if (class(frag_gr) != "GRanges") {
                stop("frag_gr is not a GRanges object!")
        }
        
        # Get the insert sizes
        widths <- width(frag_gr)
        
        graphics::hist(widths, xlab = "Insert size (bases)",
             breaks = 1000, ...)
}

