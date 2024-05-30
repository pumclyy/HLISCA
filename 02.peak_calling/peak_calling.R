##############################################################################
# Script information                                                      
# Title: Peak calling
# Author: Erping Long
# Date: 2022-02-10
# Description: None
##############################################################################

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(readr)
library(harmony)
library(GenomicRanges)
library(SeuratData)

####ATAC analysis

lung_4 <- readRDS(file = "lung_1.rds") 

DefaultAssay(lung_4) <- "ATAC"

# call peaks using MACS2

peaks <- CallPeaks(lung_4, group.by = "final_cluster")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

write.csv(peaks,"peaks_by_cell_types.csv")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(lung_4),
  features = peaks,
  cells = colnames(lung_4)
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# create a new assay using the MACS2 peak set and add it to the Seurat object
lung_4[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c(":", "-"),
  fragments = Fragments(lung_4),
  annotation = annotations
)

DefaultAssay(lung_4) <- "peaks"

lung_4 <- FindTopFeatures(lung_4, min.cutoff = 5)
lung_4 <- RunTFIDF(lung_4)
lung_4 <- RunSVD(lung_4)

saveRDS(lung_4, file = "lung.rds")
