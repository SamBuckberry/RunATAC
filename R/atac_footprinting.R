
#' Get the positions of the PWM in ATAC-seq peaks.
#' @param gr A GRanges object. Usually ranges of ATAC-seq peaks. 
#' @param pwm A positional weight matrix of class Matirx.
#' @param genome A BSgenome object. 
#' @param min.score Minimum alignment score for pwm. Default = "85%".
motif_gr <- function(gr, pwm, genome=Mmusculus, min.score="85%"){
        
        if(class(gr) != "GRanges"){
                stop("Object gr is not of class GRanges.")
        }
        if(class(pwm) != "matrix"){
                stop("Object pwm is not of class matrix.")
        }
        if(class(genome) != "BSgenome"){
                stop("Object genome is not of class BSgenome.")
        }
        if(class(min.score) != "character"){
                stop("Object genome is not of class character.")
        }
        
        # Get the nucleotide sequences
        sequences <- BSgenome::getSeq(x = genome, names=gr)
        
        # Function to find motif in one sequence
        find_motif_start <- function(x)
        {
                motif_starts <- matchPWM(pwm = pwm, subject = sequences[[x]],
                                         min.score = min.score) %>% start()
                starts <- start(gr[x]) + motif_starts
                if (length(starts) == 0){
                        out <- NULL
                } else {
                        ends <- starts + ncol(pwm)
                        out <- GRanges(seqnames = seqnames(gr[x]),
                                       ranges = IRanges(start = starts,
                                                        end = ends))
                }
                
                return(out)
        }
        
        # Apply the function across all sequences
        motif_ranges <- lapply(X = 1:length(sequences), FUN = find_motif_start)
        
        # Remove the NULLs from where no motif was detected
        is.NullOb <- function(x) is.null(motif_ranges[x]) | all(sapply(motif_ranges[x], is.null))
        no_keep <- lapply(1:length(motif_ranges), is.NullOb) %>% unlist() 
        
        # Convert to Granges object
        motif_ranges <- motif_ranges[!no_keep] %>% unlist() %>% GRangesList() %>% unlist()
        
        # Decrease the width to motif centre
        motif_ranges <- resize(motif_ranges, width = 2, fix = 'center')
        
        # Return GRanges object
        return(motif_ranges)
}

#' Calculate distances of Tn5 insertions from motif.
#' @param genome Object of class BSgenome.
#' @param tn_signal A GRanges object with tn5 insertions. See the read_atac_bam_tn function.
#' @param range The distance in bases flanking the the motif to calculate signal.
atac_motif_dist <- function(tn_signal_gr, motif_pos_gr, range=200){

        # Resize to range of signal collection
        motif_pos_gr <- GenomicRanges::resize(motif_pos_gr, width = range*2, fix = 'center')
        
        # Subset the Tn5 signal
        tn_signal_gr <- tn_signal_gr[overlapsAny(query = tn_signal_gr, subject = motif_pos_gr)]
        
        #For each motif, get the distances of nearby Tn5 insertions
        get_tn_dist <- function(x){
                
                # Subset to one region
                motif_pos <- motif_pos_gr[x]
                
                # Get the tn5 hits for range
                tn_sub <- tn_signal_gr[overlapsAny(query = tn_signal_gr, subject = motif_pos)]
                
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


#' Plot scaled ATAC-seq footprint signal.
#' @param tn_dist A vector of distances of tn5 insertions relative to motif centre.
#' See atac_motif_dist function.
#' @param window The plotting window in bases. 
#' @param ... Optional arguments passed to plot function.
plot_scaled_footprint <- function(tn_dist, window=200, ...){
        
        # Subset for window
        tn_dist <- tn_dist[abs(tn_dist) <= window/2]
        
        # Histogram
        h <- hist(tn_dist, plot = FALSE, breaks = window)
        
        #Scale the histogram from 0-1
        maxs <- max(h$counts)
        mins <- min(h$counts)
        y_scaled <- scale(h$counts, center = mins, scale = maxs - mins)
        
        # Generate plot
        plot(h$mids, y_scaled, type="l", ...,
             ylab="Fraction of signal",
             xlab="Distance to motif centre (bp)")
        }

