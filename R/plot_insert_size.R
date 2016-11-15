#' Plot insert size historgram
#' 
#' @param frag_gr A GRanges object of ATAC-seq fragments
#' @return A plot to the current graphics device
#' @export
#' @examples
#' \dontrun{
#' 
#' frags <- read_atac_frags("chr19.bam")
#' plot_insert_size(frags)
#' }
plot_insert_size <- function(frag_gr, ...){
        
        # Check inputs
        if (class(frag_gr) != "GRanges") {
                stop("frag_gr is not a GRanges object!")
        }
        
        # Get the insert sizes
        widths <- width(frag_gr)
        
        hist(widths, xlab = "Insert size (bases)",
             breaks = 1000, main="", ...)
}
