library(RunATAC)
library(magrittr)


ins_bw <- "inst/extdata/chrIV.ins.bigwig"
peaks <- read_bed("inst/extdata/chrIV_peaks.narrowPeak")
reb1 <- getMatrixByName(JASPAR2016, name=c("REB1"))
pwm <- reb1@profileMatrix
genome <- Scerevisiae

ips_pks <- "~/polo_iPSC/ATACseq/processed_data/atac_peaks/merged_replicate_peaks/atac_iPSC_combined_replicates_peaks.narrowPeak"
ips_pks <- read_bed(ips_pks)
peaks <- ips_pks[seqnames(ips_pks) %in% "chr1"]

ins_bw <- "/Volumes/Datasets/atac_iPSC_combined_replicates.ins.bigwig"
nuc_bw <- "/Volumes/Datasets/atac_iPSC_combined_replicates.occ.bigwig"

os <- getMatrixByName(JASPAR2016, name=c("Pou5f1::Sox2"))
pwm <- os@profileMatrix

genome <- Mmusculus
min.score="75%"

bias_bw <- "/Volumes/Datasets/mm10_bias.Scores.bigwig"

fp_centipede <- function(ins_bw, peaks, pwm, genome,
                         min_base_cov=0.1, min.score="75%"){
        
        # Get the motif match positions
        motif_pos <- motif_gr(gr = peaks, pwm = pwm, genome = genome,
                              min.score=min.score)
        
        # Get the pileup data of Tn5 insertions
        dat <- range_summary(bigwig = ins_bw, gr = motif_pos, range = 100)
        
        # Filter regions with < min fraction of bases covered
        covered_bases <- dat > 0
        pass_cov <- rowSums(covered_bases) >= floor(min_base_cov * ncol(dat))
        dat <- dat[pass_cov, ]
        
        # Get the nucleosome scores
        #nuc <- range_summary(bigwig = nuc_bw, gr = motif_pos, range = 20) %>% rowMeans()
        #nuc_score <- -nuc + 1
        
        # Get bias scores
        #bias <- range_summary(bigwig = bias_bw, gr = motif_pos)
        
        #ins_over_bias <- dat / bias
        
        # Conservation scores
        cons_track <- "http://cpebrazor.ivec.org/public/listerlab/sam/polo_mm_iPSC/resources/mm10.60way.phyloP60way.bw"
        cons <- range_summary(bigwig = cons_track, gr = motif_pos, range = ceiling(ncol(pwm)/2))
        cons <- rowMeans(cons)
        
        # Setup the priors
        priors <- data.frame(X=1,
                             PWMscore=log(score(motif_pos[pass_cov])),
                             Cons=cons[pass_cov]) %>% 
                as.matrix()
        
        # Fit the CENTIPEDE model
        centFit <- fitCentipede(Xlist = list(ATAC=as.matrix(cbind(dat, dat[ ,ncol(dat):1]))),
                                Y = priors)
        
        plotProfile(centFit$LambdaParList[[1]],Mlen=1)
        
        # Format the results
        df <- data.frame(chr=seqnames(motif_pos), start=start(motif_pos), end=end(motif_pos),
                         PostPr=centFit$PostPr, LogRatios=centFit$LogRatios)
        
        results <- list(centFit=df, countMat=dat, motifPos=motif_pos)
        
        return(results)
}


fp1 <- fp_centipede(ins_bw, peaks, pwm, genome)


# iPSC footprinting for Oct/Sox
library(BSgenome.Mmusculus.UCSC.mm10)
ips_pks <- "~/polo_iPSC/ATACseq/processed_data/atac_peaks/merged_replicate_peaks/atac_iPSC_combined_replicates_peaks.narrowPeak"
ips_pks <- read_bed(ips_pks)
ips_pks <- ips_pks[seqnames(ips_pks) %in% "chr1"]

ips_ins <- "/Volumes/Datasets/atac_iPSC_combined_replicates.ins.bigwig"

os <- getMatrixByName(JASPAR2016, name=c("Pou5f1::Sox2"))
os <- os@profileMatrix
        
        
os_ips_fp <- fp_centipede(ins_bw = ips_ins, peaks = ips_pks, pwm = os,
                          genome = Mmusculus, min.score="75%")


ins_bw <- "/Volumes/Datasets/atac_d9_combined_replicates.ins.bigwig"
peaks <- ips_pks
pwm <- os
sox_fl <- "~/polo_iPSC/ChIPseq/processed_data/macs_peaks_replicates/d9_sox_peaks.narrowPeak"
oct_fl <- "~/polo_iPSC/ChIPseq/processed_data/macs_peaks_replicates/d9_oct_peaks.narrowPeak"
sox_peaks <- read_bed(sox_fl)
oct_peaks <- read_bed(oct_fl)
sox_peaks <- sox_peaks[overlapsAny(sox_peaks, oct_peaks)]
oct_peaks <- oct_peaks[overlapsAny(oct_peaks, sox_peaks)]
chip_peaks <- c(sox_peaks, oct_peaks) %>% reduce()

library(caret)
test_fp <- function(peaks, ins_bw, chip_peaks, pwm){

        fp <- fp_centipede(ins_bw = ins_bw, 
                           peaks = peaks, pwm = pwm,
                           genome = Mmusculus, min.score="75%")
        
        # Get the motifs bound by TF 
        chip_hits <- overlapsAny(fp$motifPos, chip_peaks)
        
        # FP predicted binding
        fp_hits <- overlapsAny(fp$motifPos, fp$motifPos[fp$centFit$PostPr > 0.99])
        
        comp <- data.frame(chip=chip_hits, atac=fp_hits)
        
        # Create the matrix for PPV calculation
        A <- sum(comp$atac == TRUE & comp$chip == TRUE)
        B <- sum(comp$atac == TRUE & comp$chip == FALSE)
        C <- sum(comp$atac == FALSE & comp$chip == TRUE)
        D <- sum(comp$atac == FALSE & comp$chip == FALSE)
        
        dat <- matrix(c(A, B, C, D), byrow = TRUE, nrow = 2)
        
        res <- c(posP=posPredValue(dat),
                 negP=negPredValue(dat),
                 sens=(A/(A+C)),
                 spec=(D/(D+B)), 
                 n_motif=length(fp$motifPos),
                 chip_hits=sum(chip_hits),
                 fp_hits=sum(fp_hits),
                 overlap=A)
        return(res)
        
}

test_1 <- test_fp(peaks, ins_bw, chip_peaks, pwm)



m2 <- "http://cpebrazor.ivec.org/public/listerlab/sam/polo_mm_iPSC/atac/IPS_05673_M2_ATAC.ins.bw"

test_m2 <- test_fp(peaks, ins_bw=m2, nuc_bw, chip_peaks, pwm) # No PWM scores

bam_fl <- system.file("extdata", "chrIV.bam", package = "RunATAC")
frags <- read_atac_frags(bam_file = bam_fl)

ins <- read_atac_insertions(bam_fl)
ins

peak_fl <- system.file("extdata", "chrIV_peaks.narrowPeak", package = "RunATAC")
peaks <- read_bed(peak_fl)
peaks

library(JASPAR2016)
library(TFBSTools)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

reb1 <- getMatrixByName(JASPAR2016, name=c("REB1"))
reb1 <- reb1@profileMatrix


motif_pos <- motif_gr(gr = peaks, pwm = reb1, genome = Scerevisiae)

motif_pos


ins_bw <- system.file("extdata", "chrIV.ins.bigwig", package = "RunATAC")
dat <- range_summary(bigwig = ins_bw, gr = motif_pos)


# Calculate the bias scores
ins_pwm <- calc_ins_pwm(ins_gr = ins, genome = Scerevisiae)

range=100
regions_gr <- GenomicRanges::resize(motif_pos, width = range+(range+1), 
                                   fix = 'center')

bias <- calc_ins_bias(ins_pwm = ins_pwm, regions_gr = regions_gr,
                      genome = Scerevisiae)
bias <- do.call(rbind, bias)
plot(colMeans(bias))

# Correct for bias
dat <- dat / bias
plot(colSums(dat))

anno_dat <- data.frame(X=1, PWMscore=score(motif_pos)) %>% as.matrix()

library(CENTIPEDE)
data(NRSFcuts, package='CENTIPEDE')
data("NRSF_Anno")
<<<<<<< HEAD
centFit <- fitCentipede(Xlist = list(DNase=as.matrix(cbind(dat,dat))),
                        Y = anno_dat)
=======
centFit <- fitCentipede(Xlist = list(DNase=as.matrix(cbind(dat))),
                        Y = anno_dat)
centFit2 <- fitCentipede(Xlist = list(DNase=as.matrix(cbind(dat, dat))),
                        Y = anno_dat)
table(centFit$PostPr > 0.99)
table(centFit2$PostPr > 0.99)
>>>>>>> origin/master

plotProfile(centFit$LambdaParList[[1]],Mlen=1)

plot(centFit$LambdaParList$DNase[1:201, ])


peak_fl <- "~/polo_iPSC/ATACseq/processed_data/ATAC_peaks/merged_replicate_peaks/atac_iPSC_combined_replicates_peaks.narrowPeak"
peaks <- read_bed(peak_fl)
peaks

os <- getMatrixByName(JASPAR2016, name=c("Pou5f1::Sox2"))
os <- os@profileMatrix

library(BSgenome.Mmusculus.UCSC.mm10)
motif_pos <- motif_gr(gr = peaks, pwm = os, genome = Mmusculus)
motif_pos

ins_bw <- "/Volumes/Datasets/atac_iPSC_combined_replicates.ins.bigwig"
dat <- range_summary(bigwig = ins_bw, gr = motif_pos)

# Calculate the bias scores
ins_pwm <- calc_ins_pwm(ins_gr = ins, genome = Mmusculus)

range=100
regions_gr <- GenomicRanges::resize(motif_pos, width = range+(range+1), 
                                    fix = 'center')

bias <- calc_ins_bias(ins_pwm = ins_pwm, regions_gr = regions_gr,
                      genome = Mmusculus)
bias <- do.call(rbind, bias)
plot(colMeans(bias))

# Correct for bias
dat <- dat / bias
plot(colSums(dat))

anno_dat <- data.frame(X=1, PWMscore=score(motif_pos)) %>% as.matrix()

library(CENTIPEDE)
data(NRSFcuts, package='CENTIPEDE')
data("NRSF_Anno")
centFit <- fitCentipede(Xlist = list(DNase=as.matrix(cbind(dat,dat))),
                        Y = anno_dat)

plotProfile(centFit$LambdaParList[[1]],Mlen=1)

plot(centFit$LambdaParList$DNase[1:201, ])



