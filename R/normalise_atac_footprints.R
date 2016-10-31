
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


#' Get the signal from bigwig file for defined regions.
#' 
#' @param bigwig Path to a bigwig file. Also accpets http paths. 
#' @param gr A genomic ranges object of regions to obtain signal. All ranges should be the same width.
#' @param range The distance from the centre of gr to flank. Default=100. 
#' @param log Logical. Is the signal in the bigwig file in log space? Default=FALSE.
#' If log = TRUE, the scores in the bigwig file are transformed by e^scores.
#' @return A numeric vector of aggregate signal if aggregate = TRUE, which is the default. 
#' If aggregate = FALSE a data.frame of signal of all regions is returned.
#' @export
#' @import magrittr
#' @import GenomicRanges
#' @import IRanges
#' @import tidyr
#' @import rtracklayer
range_summary <- function(bigwig, gr, range=100, log=FALSE, aggregate=TRUE){
        
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

atac_sig  <- "~/Desktop/atac_iPSC_combined_replicates.ins.bigwig"
atac_bias <- "http://cpebrazor.ivec.org/public/listerlab/sam/polo_mm_iPSC/atac/mm10_bias.Scores.bigwig"
regions <- "~/polo_iPSC/ATACseq/processed_data/atac_cluster_peaks/c_means_peaks/cluster_1.bed"
tf <- read.table("~/R_packages/RunATAC/inst/exdata/ctcf.pwm") %>% as.matrix()





regions <- read_bed(regions) 
#regions <- regions[seqnames(regions) == "chr19"]
regions_tf <-  motif_gr(regions, pwm = tf)

ins <- range_summary(bigwig = atac_sig, gr = regions_tf, range = 150, log = FALSE)
bias <- range_summary(bigwig = atac_bias, gr = regions_tf, range = 150, log = TRUE)

ins_full <- range_summary(bigwig = atac_sig, gr = regions_tf, range = 150,
                          log = FALSE, aggregate = FALSE)
bias_full <- range_summary(bigwig = atac_bias, gr = regions_tf, range = 150,
                           log = TRUE, aggregate = FALSE)


plot(ins$position, ins$signal, type='l')
plot(bias$position, bias$signal, type='l')

plot(bias$position, ins$signal / bias$signal, type='l')


load(file = "~/polo_iPSC/resources/pwm_matrix_list.Rda")
# Get mofif PWM for Oct4-Sox2
os_pwm <- pwm_matrix_list$MA0142.1

regions <- ("~/polo_iPSC/ATACseq/processed_data/atac_cluster_peaks/c_means_peaks/cluster_1.bed")
regions <- read_bed(regions)
reg_os <-  motif_gr(regions, pwm = os_pwm, min.score = "70%")

ins <- range_summary(bigwig = atac_sig, gr = reg_os, range = 150, log = FALSE, aggregate = FALSE)
bias <- range_summary(bigwig = atac_bias, gr = reg_os, range = 150, log = TRUE, aggregate = FALSE)
nuc_sig <- "http://cpebrazor.ivec.org/public/listerlab/sam/polo_mm_iPSC/atac/atac_d12_combined_replicates.nucleoatac_signal.bigwig"
nuc <- range_summary(bigwig = nuc_sig, gr = reg_os, range = 500, aggregate = FALSE)
o_sig <- range_summary(bigwig = "http://cpebrazor.ivec.org/public/listerlab/sam/polo_mm_iPSC/chip_seq/ips_oct_subtract.bw",
                       gr = reg_os, range = 500, log = FALSE, aggregate = FALSE)
s_sig <- range_summary(bigwig = "http://cpebrazor.ivec.org/public/listerlab/sam/polo_mm_iPSC/chip_seq/ips_sox_subtract.bw",
                       gr = reg_os, range = 500, log = FALSE, aggregate = FALSE)

mc <- range_summary(bigwig = "http://cpebrazor.ivec.org/public/listerlab/sam/polo_mm_iPSC/methylCseq/mmiPS_6__p1GFPpos.CG.level.unstranded.bigwig",
                     gr = reg_os, range = 500, log = FALSE, aggregate = FALSE)

# Normalise insertions by bias signal
ins_over_bias <- ins  / bias


## Get the nuc occupancy class
occ <- range_summary(bigwig = "http://cpebrazor.ivec.org/public/listerlab/sam/polo_mm_iPSC/atac/atac_d12_combined_replicates.occ.bigwig",
                     gr = reg_os, range = 10, log = FALSE, aggregate = FALSE)

occ_means <- rowMeans(occ)

set_occ_class <- function(x){
        
        dat <- NULL
        
        if (x < 0.1) {
                dat <- "A"
        } else if (x >= 0.1 & x < 0.3) {
                dat <- "B" 
        } else if (x >= 0.3 & x < 0.5) {
                dat <- "C"
        } else if (x > 0.5) { 
                dat <- "D"
        }
        
        return(dat)
}

occ_class <- lapply(occ_means, set_occ_class) %>% unlist()


facet_plot_occ_class <- function(dat, occ_class, y_lab=""){
        
        # set the occ class
        dat$class <- occ_class
        
        dat_melt <- melt(dat, id.vars = c('class'))
        
        dat_mean <- dat_melt %>%
                group_by(variable, class) %>%
                summarise(yy=mean(value))
        
        dat_mean$variable <- as.character(dat_mean$variable) %>% as.numeric()
        dash <- 0
        ggplot(dat_mean, aes(x = variable, y = yy, group=class)) +
                geom_line() +
                facet_grid(class~.) +
                ylab(y_lab) +
                xlab("Position relative to motif") +
                theme_bw() +
                theme(plot.background = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.text = element_text(color = 'black'),
                      axis.line = element_line(),
                      #text = element_text(size=8),
                      axis.line.x = element_line(color = 'black'),
                      axis.line.y = element_line(color = 'black'))
}


gg_nuc <- facet_plot_occ_class(nuc, occ_class = occ_class, y_lab = "Nucleosome signal")
gg_nuc
gg_oct <- facet_plot_occ_class(o_sig, occ_class = occ_class, y_lab = "Oct4 ChIP-seq (CPM)")
gg_ins <- facet_plot_occ_class(ins_over_bias, occ_class = occ_class, y_lab = "ATAC-seq insertions / insertion bias")

library(cowplot)
pdf("~/Desktop/nucleosome_plots.pdf", width = 7.5, height = 5)
plot_grid(gg_nuc, 
          gg_ins, 
          gg_oct,
          nrow=1, ncol=3)
dev.off()


