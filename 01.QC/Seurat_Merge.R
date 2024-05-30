##############################################################################
# Script information                                                      
# Title: Merging data from different samples
# Author: Erping Long
# Date: 2023-01-06
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

lung_M1 <- readRDS(file = "M1_lung_4.rds")
lung_M3 <- readRDS(file = "M3_lung_4.rds")
lung_M4 <- readRDS(file = "M4_lung_4.rds")
lung_M5 <- readRDS(file = "M5_lung_4.rds")
lung_M6 <- readRDS(file = "M6_lung_4.rds")
lung_M7 <- readRDS(file = "M7_lung_4.rds")
lung_M8 <- readRDS(file = "M8_lung_4.rds")
lung_M9 <- readRDS(file = "M9_lung_4.rds")

lung_NCI_13 <- readRDS(file = "NCI_13_lung_4.rds")
lung_NCI_28 <- readRDS(file = "NCI_28_lung_4.rds")
lung_NCI_55 <- readRDS(file = "NCI_55_lung_4.rds")
lung_NCI_115 <- readRDS(file = "NCI_115_lung_4.rds")
lung_NCI_118 <- readRDS(file = "NCI_118_lung_4.rds")
lung_S1 <- readRDS(file = "S1_lung_4.rds")
lung_S2 <- readRDS(file = "S2_lung_4.rds")
lung_N17_20 <- readRDS(file = "N17_20_lung_4.rds")

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
