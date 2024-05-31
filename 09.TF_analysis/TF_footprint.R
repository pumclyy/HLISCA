##############################################################################
# Script information                                                      
# Title: TF footprinting
# Author: Erping Long
# Date: 2023-5-17
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
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

lung_peaks <- readRDS(file = "lung.rds")

DefaultAssay(lung_peaks) <- "peaks"

#### TF Footprinting 
lung_peaks <- Footprint(
  object = lung_peaks,
  motif.name = c("EGR1","EGR2","EGR3","EGR4","EHF"),
  genome = BSgenome.Hsapiens.UCSC.hg38)

p2 <- PlotFootprint(lung_peaks, features = c("EGR1","EGR2","EGR3","EGR4","EHF"))

p2 + patchwork::plot_layout(ncol = 1)

ggsave("footprint_EGR1_EGR2_EGR3_EGR4_EHF.pdf", 
       p2+ patchwork::plot_layout(ncol = 1), width = 10, height = 30)


