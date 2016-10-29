
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
                motif_starts_pos <- Biostrings::matchPWM(pwm = pwm, subject = sequences[[x]],
                                                     min.score = min.score) %>% start()
                
                motif_starts_neg <- Biostrings::matchPWM(pwm = reverseComplement(pwm),
                                                         subject = sequences[[x]],
                                                         min.score = min.score) %>% start()
                
                starts_pos <- start(gr[x]) + motif_starts_pos
                starts_neg <- start(gr[x]) + motif_starts_neg
                starts <- c(starts_pos, starts_neg)
                
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

atac_region_summary <- function(bigwig, regions, range=100, log=FALSE){
        
        # Resize to range of signal collection
        region_gr <- GenomicRanges::resize(regions, width = range+(range+1), 
                                           fix = 'center')
        
        # Compute bias signal -----------------------------------------------------
        
        # Get the signal for each base for the defined regions
        # Reurns a SimpleNumericList
        scores <- rtracklayer::import(con = bigwig, format = "Bigwig",
                                      which=region_gr, as="NumericList")
        
        # Transform scores from Log if log == TRUE
        if(log == TRUE){
                scores <- exp(scores)
        }
        
        scores <- as.data.frame(scores)
        
        # Set the relative postion
        range_intervals <- (0-range):(0+range)
        scores$position <- rep(range_intervals, times=length(unique(scores$group)))
        
        scores <- scores[ ,c(4,3)]
        
        agg_scores <- aggregate(. ~ position, scores, sum)
        
        return(agg_scores)
}



atac_sig  <- "~/Desktop/atac_iPSC_combined_replicates.ins.bigwig"
atac_bias <- "~/Desktop/mm10_bias_chr19.Scores.bigwig"
regions <- "~/polo_iPSC/ATACseq/processed_data/atac_cluster_peaks/c_means_peaks/cluster_1.bed"
tf <- read.table("~/R_packages/RunATAC/inst/exdata/ctcf.pwm") %>% as.matrix()





regions <- read_bed(regions) 
regions <- regions[seqnames(regions) == "chr19"]
regions <-  motif_gr(regions, pwm = tf)

ins <- atac_region_summary(bigwig = atac_sig, regions = regions, range = 150, log = FALSE)
bias <- atac_region_summary(bigwig = atac_bias, regions = regions, range = 150, log = TRUE)

plot(ins$position, ins$value, type='l')
plot(bias$position, bias$value, type='l')

plot(bias$position, ins$value / bias$value, type='l')


load(file = "~/polo_iPSC/resources/pwm_matrix_list.Rda")
# Get mofif PWM for Oct4-Sox2
os_pwm <- pwm_matrix_list$MA0142.1

regions <- ("~/polo_iPSC/ATACseq/processed_data/atac_cluster_peaks/c_means_peaks/cluster_2.bed")
regions <- read_bed(regions)
reg_os <-  motif_gr(regions, pwm = os_pwm, min.score = "70%")

ins <- atac_region_summary(bigwig = atac_sig, regions = reg_os, range = 300, log = FALSE)
bias <- atac_region_summary(bigwig = atac_bias, regions = regions, range = 300, log = TRUE)

plot(ins$position, ins$value, type='l')
plot(bias$position, bias$value, type='l')

plot(bias$position, ins$value / bias$value, type='l')



