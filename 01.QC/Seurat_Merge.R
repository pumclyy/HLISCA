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

lung_MN1 <- readRDS(file = "MN1_lung_4.rds")
lung_MS1 <- readRDS(file = "MS1_lung_4.rds")
lung_MS2 <- readRDS(file = "MS2_lung_4.rds")
lung_MS3 <- readRDS(file = "MS3_lung_4.rds")
lung_MS4 <- readRDS(file = "MS4_lung_4.rds")
lung_MN2 <- readRDS(file = "MN2_lung_4.rds")
lung_MN3 <- readRDS(file = "MN3_lung_4.rds")
lung_MN4 <- readRDS(file = "MN4_lung_4.rds")

lung_FN1 <- readRDS(file = "FN1_lung_4.rds")
lung_FN2 <- readRDS(file = "FN2_lung_4.rds")
lung_FS1 <- readRDS(file = "FS1_lung_4.rds")
lung_FN3 <- readRDS(file = "FN3_lung_4.rds")
lung_FS2 <- readRDS(file = "FS2_lung_4.rds")
lung_FS3 <- readRDS(file = "FS3_lung_4.rds")
lung_FS4 <- readRDS(file = "FS4_lung_4.rds")
lung_FN4 <- readRDS(file = "FN4_lung_4.rds")

lung_MN1@meta.data$orig.ident = c("MN1")
lung_MS1@meta.data$orig.ident = c("MS1")
lung_MS2@meta.data$orig.ident = c("MS2")
lung_MS3@meta.data$orig.ident = c("MS3")
lung_MS4@meta.data$orig.ident = c("MS4")
lung_MN2@meta.data$orig.ident = c("MN2")
lung_MN3@meta.data$orig.ident = c("MN3")
lung_MN4@meta.data$orig.ident = c("MN4")
lung_FN1@meta.data$orig.ident = c("FN1")
lung_FN2@meta.data$orig.ident = c("FN2")
lung_FS1@meta.data$orig.ident = c("FS1")
lung_FN3@meta.data$orig.ident = c("FN3")
lung_FS2@meta.data$orig.ident = c("FS2")
lung_FS3@meta.data$orig.ident = c("FS3")
lung_FS4@meta.data$orig.ident = c("FS4")
lung_FN4@meta.data$orig.ident = c("FN4")

lung_4 <- merge(lung_MN1,
                  y = c(lung_MS1, lung_MS2, lung_MS3, lung_MS4, lung_MN2, lung_MN3, lung_MN4, lung_FN1, lung_FN2, lung_FS1, lung_FN3, lung_FS2, lung_FS3, lung_FS4, lung_FN4),
                  add.cell.ids = c("MN1", "MS1", "MS2", "MS3", "MS4", "MN2", "MN3", "MN4", "FN1", "FN2", "FS1", "FN3", "FS2", "FS3", "FS4", "FN4"), 
                  project = "lung_Merge")

dim(lung_4)

table(lung_4@meta.data$orig.ident)

saveRDS(lung_4, file = "/data/Choi_lung/SHARE-seq/Seurat_SHARE/lung_merge_all.rds")
