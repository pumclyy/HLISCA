##############################################################################
# Script information                                                      
# Title: cCRE-gene linkage
# Author: Erping Long, Yuyan Li
# Date: 2023-2-17
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
library(BSgenome.Hsapiens.UCSC.hg38)

lung <- readRDS(file = "lung.rds") 

####Linking peaks to genes

DefaultAssay(lung) <- "peaks"

# first compute the GC content for each peak
lung <- RegionStats(lung, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
lung <- LinkPeaks(
  object = lung,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 1e+06,
  min.cells = 10,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  verbose = TRUE
  )

write.csv(Links(lung[["peaks"]]),"./linkage_peaks_by_cell_types_SCT_3.csv")

saveRDS(lung, file = "lung_linkage_1Mb.rds")

# 2Mb

lung <- RegionStats(lung, genome = BSgenome.Hsapiens.UCSC.hg38)

lung <- LinkPeaks(
  object = lung,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 2e+06
)

saveRDS(lung, "lung_linkage_2Mb.rds")

write.csv(Links(lung[["peaks"]]), "links_2Mb.csv")

# 5Mb

lung <- RegionStats(lung, genome = BSgenome.Hsapiens.UCSC.hg38)

lung <- LinkPeaks(
  object = lung,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 5e+06
)

saveRDS(lung, "lung_linkage_5Mb.rds")

write.csv(Links(lung[["peaks"]]), "links_5Mb.csv")

