library(RunATAC)
library(CENTIPEDE)

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

bias <- calc_ins_bias(ins_pwm = ins_pwm, regions_gr = regions_gr, genome = Scerevisiae)
bias <- do.call(rbind, bias)
plot(colMeans(bias))

# Correct for bias
dat <- dat / bias

anno_dat <- data.frame(X=1, PWMscore=score(motif_pos)) %>% as.matrix()

library(CENTIPEDE)
data(NRSFcuts, package='CENTIPEDE')
data("NRSF_Anno")
centFit <- fitCentipede(Xlist = list(DNase=as.matrix(cbind(dat,dat))),
                        Y = anno_dat, )

plotProfile(centFit$LambdaParList[[1]],Mlen=1)

plot(centFit$LambdaParList$DNase[1:201, ])


centFit <- fitCentipede(Xlist = list(DNase=as.matrix(NRSFcuts)),
                        Y=cbind(rep(1, dim(NRSF_Anno)[1]), NRSF_Anno[,5], NRSF_Anno[,6]))
