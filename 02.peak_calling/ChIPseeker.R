##############################################################################
# Script information                                                      
# Title: ChIPseeker
# Author: Erping Long
# Date: 2022-02-10
# Description: None
##############################################################################

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#read the results of peak calling

data <- read.csv("peaks_by_cell_types.csv")

#change the result to GRanges object

peaks <- with(data, GRanges(seqnames = seqnames, IRanges(start, end), strand))

#annotate peaks

peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

write.csv(peakAnno,"peaks_annotation_by_ChIPseeker.csv")

