##############################################################################
# Script information                                                      
# Title: Merging data from different samples
# Author: Erping Long
# Date: 2021-03-15
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

lung_M1 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/M1/M1_lung_4.rds")
lung_M3 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/M3/M3_lung_4.rds")
lung_M4 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/M4/M4_lung_4.rds")
lung_M5 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/M5/M5_lung_4.rds")
lung_M6 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/M6/M6_lung_4.rds")
lung_M7 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/M7/M7_lung_4.rds")
lung_M8 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/M8/M8_lung_4.rds")
lung_M9 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/M9/M9_lung_4.rds")

lung_NCI_13 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/NCI_13/NCI_13_lung_4.rds")
lung_NCI_28 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/NCI_28/NCI_28_lung_4.rds")
lung_NCI_55 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/NCI_55/NCI_55_lung_4.rds")
lung_NCI_115 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/NCI_115/NCI_115_lung_4.rds")
lung_NCI_118 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/NCI_118/NCI_118_lung_4.rds")
lung_S1 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/S1/S1_lung_4.rds")
lung_S2 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/S2/S2_lung_4.rds")
lung_N17_20 <- readRDS(file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/N17_20/N17_20_lung_4.rds")

lung_M1@meta.data$orig.ident = c("M1")
lung_M3@meta.data$orig.ident = c("M3")
lung_M4@meta.data$orig.ident = c("M4")
lung_M5@meta.data$orig.ident = c("M5")
lung_M6@meta.data$orig.ident = c("M6")
lung_M7@meta.data$orig.ident = c("M7")
lung_M8@meta.data$orig.ident = c("M8")
lung_M9@meta.data$orig.ident = c("M9")
lung_NCI_13@meta.data$orig.ident = c("NCI_13")
lung_NCI_28@meta.data$orig.ident = c("NCI_28")
lung_NCI_55@meta.data$orig.ident = c("NCI_55")
lung_NCI_115@meta.data$orig.ident = c("NCI_115")
lung_NCI_118@meta.data$orig.ident = c("NCI_118")
lung_S1@meta.data$orig.ident = c("S1")
lung_S2@meta.data$orig.ident = c("S2")
lung_N17_20@meta.data$orig.ident = c("N17_20")

lung_4 <- merge(lung_M1,
                  y = c(lung_M3, lung_M4, lung_M5, lung_M6, lung_M7, lung_M8, lung_M9, lung_NCI_13, lung_NCI_28, lung_NCI_55, lung_NCI_115, lung_NCI_118, lung_S1, lung_S2, lung_N17_20),
                  add.cell.ids = c("M1", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "NCI_13", "NCI_28", "NCI_55", "NCI_115", "NCI_118", "S1", "S2", "N17_20"), 
                  project = "lung_Merge")

dim(lung_4)

table(lung_4@meta.data$orig.ident)

saveRDS(lung_4, file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/lung_merge_all.rds")
