#' Get the signal from bigwig file for defined regions.
#' 
#' @param bigwig Path to a bigwig file. Also accpets http paths. 
#' @param gr A genomic ranges object of regions to obtain signal. All ranges should be the same width.
#' @param range The distance from the centre of gr to flank. Default=100. 
#' @param log Logical. Is the signal in the bigwig file in log space? Default=FALSE.
#' If log = TRUE, the scores in the bigwig file are transformed by e^scores.
#' @param aggregate Logical. Should data for each position be added to get an aggregate score?
#' @return A numeric vector of aggregate signal if aggregate = TRUE, which is the default. 
#' If aggregate = FALSE a data.frame of signal of all regions is returned.
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import tidyr
#' @import rtracklayer
range_summary <- function(bigwig, gr, range=100, log=FALSE, aggregate=FALSE){
        
        # Make all GRanges same width by fixing centre
        region_gr <- GenomicRanges::resize(gr, width = range+(range+1), 
                                           fix = 'center')
        
        # Get the signal for each base for the defined regions
        # Reurns a SimpleNumericList
        scores <- rtracklayer::import(con = bigwig, format = "Bigwig",
                                      which=region_gr, as="NumericList")
        
        # Transform scores from Log if log == TRUE
        if(log == TRUE){
                scores <- exp(scores)
        }
        
        # Reformat to data.frame
        scores <- as.data.frame(scores)
        
        # Set the relative postion
        range_intervals <- (0-range):(0+range)
        scores$position <- rep(range_intervals,
                               times=length(unique(scores$group)))
        
        # Drop the group_name
        scores$group_name <- NULL
        
        # Hack to get R CMD Check to pass
        position <- NULL
        value <- NULL
        
        # Reformat data into matrix
        scores <- tidyr::spread(scores, key = position, value = value)
        rownames(scores) <- scores$group
        scores$group <- NULL
        
        # Aggregate scores if aggregate=TRUE
        if(aggregate == TRUE){
                scores <- colSums(scores)
                rel_pos <- names(scores) %>% as.character %>% as.numeric()
                scores <- data.frame(position=rel_pos,
                                     signal=scores)
                rownames(scores) <- NULL
        }
        
        return(scores)
}
