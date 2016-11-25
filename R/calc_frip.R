#' Calculate the fraction of reads in peaks.
#' 
#' @param signal A GRanges object. Typically of ATAC-seq Tn5 insertions. 
#' @param peaks A GRanges object. Typically of ATAC-seq peaks.
#' @return Fraction of ranges in signal overlapping peaks. Numeric.
#' @import IRanges
#' @export
calc_frip <- function(signal, peaks){
        
        if(class(signal) != "GRanges"){
                stop("Object gr is not of class GRanges.")
        }
        
        if(class(peaks) != "GRanges"){
                stop("Object gr is not of class GRanges.")
        }
        
        # Get the overlapping regions
        olaps <- IRanges::overlapsAny(query = signal, subject = peaks)
        
        # Calculate FRiP
        FrIP <- sum(olaps) / length(signal)
        return(FrIP)
}
