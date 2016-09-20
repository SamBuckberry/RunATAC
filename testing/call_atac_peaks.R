source("https://bioconductor.org/biocLite.R")
library(BiocInstaller)
useDevel()
biocLite(c("GenomicRanges", "GenomicAlignments", "rtracklayer", "motifRG"))
biocLite("BayesPeak")
library(chromstaR)
library(devtools)
library(RunATAC)
library(BayesPeak)
chrom_info <- fetchExtendedChromInfoFromUCSC(genome = "mm10")

# Load the tn5 insertion positions
atac_gr <- read_atac_pos(bam_file = "~/Desktop/atac_MEF_chr16_test.bam")
length(atac_gr)

#read_width <- 25
#atac_gr <- atac_gr + read_width

# df <- data.frame(atac_gr)
# colnames(df) <- c("chr", "start", "end", "strand")
# 
# rd <-as(atac_gr, "RangedData") 
# pks <- bayespeak(treatment = df, bin.size = 100, chr = "chr16",
#                  start = 1, end = 98207768)
# 
# sp <- summarise.peaks(pks)

bin_size <- 200
assembly <- "mm10"
blacklist <- read_bed("~/polo_iPSC/resources/mm10_blacklist_regions.bed")

call_atac_peaks <- function(atac_gr...){
        
        bin_gr <- chromstaR::binReads(file = atac_gr, assembly = assembly,
                                      binsizes = bin_size, blacklist = blacklist)
        bin_gr <- bin_gr[seqnames(bin_gr) %in% c("chr16")]
        
        model <- chromstaR::callPeaksUnivariate(binned.data = bin_gr,
                                                #read.cutoff.quantile = 0.99, 
                                                keep.posteriors = TRUE,
                                                prefit.on.chr = "chr16")
        
        plotHistogram(model) + ggtitle("ATAC-seq")
        
        segments <- model$segments
        peaks <- segments[segments$state ==  'modified']
        length(peaks); summary(width(peaks))
        
        
        # Change the posterior cutoff and get number of peaks
        model <- changePostCutoff(model, post.cutoff=0.999)
        segments <- model$segments
        peaks <- segments[segments$state ==  'modified' ]
        length(peaks); summary(width(peaks))
        
        peaks <- peaks[width(peaks) < 2000]
        peaks[sample(1:length(peaks), 3)]
                               
}

