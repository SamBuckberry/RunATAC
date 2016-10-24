### Get Tn5 bias track for regions

bias_bigwig <- "~/Desktop/mm10_bias_chr19.Scores.bigwig"
atac_regions <- "~/polo_iPSC/ATACseq/processed_data/atac_peaks/atac_all_peaks_union.bed"
atac_bam <- "~/Desktop/test.bam"
ctcf <- read.table("~/R_packages/RunATAC/inst/exdata/ctcf.pwm") %>% as.matrix()

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




#' Read BAM file Tn5 insertion positions into GRanges object
#' 
#' @param bam_file Path to a BAM formatted file.
#' @param yieldSize The number of records from the BAM file to read in each chunk. See ?BamFile.
#' @param ... Arguments passed to ScanBamParam.
#' @return A GRanges object of Tn5 insertion positions
#' @examples 
#' \dontrun{
#' tn <- import_atac_pos("test_100k.bam")
#' }
#' @export
read_atac_bam_tn <- function(bam_file, yieldSize=1e6, ...){
        
        # Allow the setting of parameters for Bam file scan such as Which regions
        param <- Rsamtools::ScanBamParam(...)
        
        message("Processing BAM file. This this may take a while for large files...")
        message("If you have lots of RAM, increase yieldSize to speed things up.")
        
        # Loop for reading the BAM file in chunks
        bf <- Rsamtools::BamFile(file = bam_file, yieldSize = yieldSize)
        open(bf)
        gr <- NULL
        repeat {
                chunk <- GenomicAlignments::readGAlignments(bf, param = param)
                if (length(chunk) == 0L)
                        break
                chunk_gr <- GenomicRanges::GRanges(chunk)
                if (is.null(gr)) {
                        gr <- chunk_gr
                } else {
                        gr <- c(gr, chunk_gr)
                }
        }
        close(bf)
        
        gr <- GRangesList(gr) %>% unlist()
        
        message("Offsetting alignments to Tn5 integration site.")
        
        # Offset the reads to correspond to tn5 insertion site
        pos <- gr[strand(gr) == "+"] %>% 
                GenomicRanges::shift(shift=4) %>%
                GenomicRanges::resize(width = 2, fix = "start")
        
        neg <- gr[strand(gr) == "-"] %>%
                GenomicRanges::shift(shift = -5) %>%
                GenomicRanges::resize(width = 2, fix = "start")
        
        # Return the pos and neg strands together
        gr <- c(pos, neg)
        
        message("Done!")
        
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


tn_signal_gr <- read_atac_bam_tn(bam_file = atac_bam)
atac_peaks <- read_bed(atac_regions)
atac_peaks <- atac_peaks[seqnames(atac_peaks) == "chr19"]

motif_pos_gr <- motif_gr(gr = atac_peaks, pwm = ctcf)

bias_correct_atac_motif <- function(tn_signal_gr, motif_pos_gr, bias_bigwig, range=200){
        
        # Resize to range of signal collection
        motif_pos_gr <- GenomicRanges::resize(motif_pos_gr, width = range*2,
                                              fix = 'center')
        
        # Compute bias signal -----------------------------------------------------
        
        # Get the bias signal for regions
        bias <- rtracklayer::import(con = bias_bigwig, format = "Bigwig",
                                    which=motif_pos_gr)
        
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
        
        bias_mean <- bias_scores %>%
                group_by(rel_pos) %>%
                summarise(yy=mean(bias_score))
        
        colnames(bias_mean) <- c("rel_pos", "bias_mean")

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
        results$bias_mean <- as.character(results$bias_mean) %>% as.numeric()
        results$insertions <- as.character(results$insertions) %>% as.numeric()
        
        # sort the results
        results <- results[order(results$rel_pos, decreasing = FALSE), ]
        
        # Return results
        return(results)

}



result <- bias_correct_atac_motif(tn_signal_gr, motif_pos_gr, bias_bigwig)

