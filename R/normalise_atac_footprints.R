
atac_sig  <- "~/Desktop/browser_files/atac_d12_combined_replicates.ins.bigwig"
atac_bias <- "~/Desktop/mm10_bias_chr19.Scores.bigwig"
regions <- "~/polo_iPSC/ATACseq/processed_data/atac_peaks/atac_all_peaks_union.bed"
tf <- read.table("~/R_packages/RunATAC/inst/exdata/ctcf.pwm") %>% as.matrix()


#' Read BED formatted file to GRanges object
#' 
#' @param bed_file Path Path to BED formatted file.
#' @return A GRanges object
#' @examples
#' \dontrun{
#' import_bed("myBedFile.bed")
#' }
#' @note Only imports the first 3 columns of a bed file.
#' @export
read_bed <- function(bed_file){
        dat <- data.table::fread(bed_file, sep = "\t", header = FALSE)
        gr <- GenomicRanges::GRanges(seqnames = dat$V1,
                                     ranges = IRanges(start = dat$V2,
                                                      end = dat$V3))
        return(gr)
}

#' Get the positions of the PWM in ATAC-seq peaks.
#' 
#' @param gr A GRanges object. Usually ranges of ATAC-seq peaks. 
#' @param pwm A positional weight matrix of class Matirx.
#' @param genome A BSgenome object. 
#' @param min.score Minimum alignment score for pwm. Default = "85%".
#' @return A GRanges object with positions of PWM match in peaks.
#' @export
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @import IRanges
#' @examples
#' \dontrun{
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' ctcf <- read.table("/data/ctcf.pwm") %>% as.matrix()
#' peaks <- read_bed("/data/atac_peaks.bed")
#' mpos <- motif_gr(peaks, ctcf, Mmusculus)
#' }
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
                motif_starts <- Biostrings::matchPWM(pwm = pwm, subject = sequences[[x]],
                                                     min.score = min.score) %>% start()
                starts <- start(gr[x]) + motif_starts
                if (length(starts) == 0){
                        out <- NULL
                } else {
                        ends <- starts + ncol(pwm)
                        out <- GenomicRanges::GRanges(seqnames = seqnames(gr[x]),
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

regions <- read_bed(regions) 
regions <- regions[seqnames(regions) == "chr19"]
regions <-  motif_gr(regions, pwm = tf)


bigwig <- atac_sig
log=FALSE
window=100

atac_region_summary <- function(bigwig, regions, range=100, log=FALSE){
        
        # Resize to range of signal collection
        region_gr <- GenomicRanges::resize(regions, width = range*2, 
                                              fix = 'center')
        
        # Compute bias signal -----------------------------------------------------
        
        # Get the bias signal for regions
        gr <- rtracklayer::import(con = bigwig, format = "Bigwig",
                                    which=region_gr)
        
        # Transform scores from Log if log == TRUE
        if(log == TRUE){
                gr$score <- exp(gr$score)
        }
        
        get_scores <- function(x){
                
                # Subset to one region
                motif_pos <- region_gr[x]
                
                
                ########## Need to expand ranges to one base per range!!!!!!
                
                
                
                # Get the bias calculation for range
                sub <- gr[IRanges::overlapsAny(query = gr,
                                                      subject = motif_pos)]
                
                # Calculate the relative distance to the motif centre for insertions
                bias_dist <- (start(motif_pos) + range) - start(sub)
                
                # Get bias scores for relative position
                df <- data.frame(rel_pos=bias_dist, score=sub$score)
                
                return(df)
                
        }
        
        message("Calculating scores for each region...")
        scores <- lapply(1:length(region_gr), get_scores)
        scores <- do.call(rbind, scores)
        
        # Get aggregate bias signal
        score_sum <- scores %>%
                group_by(rel_pos) %>%
                summarise(yy=sum(score))
        
        colnames(score_sum) <- c("rel_pos", "score_aggregate")
        
        return(data.frame(score_sum))
        
}


ins <- atac_region_summary(bigwig = atac_sig, regions = regions, range = 100, log = FALSE)
bias <- atac_region_summary(bigwig = atac_bias, regions = regions, range = 100, log = TRUE)











bias_correct_transform_atac_motif <- function(atac_sig, atac_bias, regions, window=100){

        # Resize to range of signal collection
        motif_pos_gr <- GenomicRanges::resize(regions, width = window, 
                                              fix = 'center')
        
        # Compute bias signal -----------------------------------------------------
        
        # Get the bias signal for regions
        bias <- rtracklayer::import(con = atac_bias, format = "Bigwig",
                                    which=motif_pos_gr)
        
        # Transform scores from Log
        bias$score <- exp(bias$score)
        
        get_bias_scores <- function(x){
                
                # Subset to one region
                motif_pos <- motif_pos_gr[x]
                
                # Get the bias calculation for range
                bias_sub <- bias[IRanges::overlapsAny(query = bias,
                                                      subject = motif_pos)]
                
                # Calculate the relative distance to the motif centre for insertions
                bias_dist <- (start(motif_pos) + range) - start(bias_sub)
                
                # Get bias scores for relative position
                df <- data.frame(rel_pos=bias_dist, bias_score=bias_sub$score)
                
                return(df)
                
        }
        
        message("Calculating Tn5 insertion bias scores...")
        bias_scores <- lapply(1:length(motif_pos_gr), get_bias_scores)
        bias_scores <- do.call(rbind, bias_scores)
        
        # Get aggregate bias signal
        bias_sum <- bias_scores %>%
                group_by(rel_pos) %>%
                summarise(yy=sum(bias_score))
        
        colnames(bias_sum) <- c("rel_pos", "bias_aggregate")
        
        # Compute Tn5 signal -----------------------------------------------------
        
        # Subset the Tn5 signal
        tn_signal_gr <- tn_signal_gr[IRanges::overlapsAny(query = tn_signal_gr,
                                                          subject = motif_pos_gr)]
        
        #For each motif region, get the distances of nearby Tn5 insertions
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
        
        # Calculate the histogram
        hist_dists <- table(dists) %>% data.frame()
        
        colnames(hist_dists) <- c("rel_pos", "insertions")
        
        # Merge insertion frequency with bias 
        results <- merge.data.frame(bias_mean, hist_dists,
                                    by = "rel_pos", all = TRUE)
        
        results$rel_pos <- as.character(results$rel_pos) %>% as.numeric()
        results$bias_aggregate <- as.character(results$bias_aggregate) %>% as.numeric()
        results$insertions <- as.character(results$insertions) %>% as.numeric()
        
        results$norm <- results$insertions / results$bias_aggregate
        
        # sort the results
        results <- results[order(results$rel_pos, decreasing = FALSE), ]
        
        # scale the insertions from 0-1
        #         range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
        #         results$insertions_scaled <- range01(results$insertions, na.rm = TRUE)
        #         results$bias_scaled <- range01(results$bias_aggregate, na.rm = TRUE)
        #         
        #         results$normalised <- results$insertions / results$bias
        # Return results
        return(results)
        
}