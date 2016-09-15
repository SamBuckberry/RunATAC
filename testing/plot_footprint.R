
# Load the PWM to plot footprints

# Load the file with regions of interest
rgns <- read_bed("~/polo_iPSC/ATACseq/processed_data/atac_cluster_peaks/atac_cluster_1_chip_sub_cluster_1.bed")

# Read in the BAM and get insertion signal for regions
inserts <- read_atac_pos(bam_file = "~/Desktop/test.bam", yieldSize = 1e6, which=rgns)
length(inserts)
