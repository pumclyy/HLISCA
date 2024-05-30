##############################################################################
# Script information                                                      
# Title: Spliting data
# Author: Yuyan Li
# Date: 2023-05-05
# Description: None
##############################################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(celldex)
library(cowplot)
library(tidyverse)
library(CellChat)

lung <- readRDS("16_multiome.rds")

split_list <- SplitObject(lung, split.by = "orig.ident") 


list_ns <- split_list[c("M1", "M7", "M8", "M9", "NCI_13", "NCI_28", "NCI_115", "N17_20")]

list_s <- split_list[c("M3", "M4", "M5", "M6", "NCI_55", "NCI_118", "S1", "S2")]


lung_ns <- merge(list_ns$M1, 
                 y = c(list_ns$M7, list_ns$M8, list_ns$M9, list_ns$NCI_13, list_ns$NCI_28, list_ns$NCI_115, list_ns$N17_20))
lung_s <- merge(list_s$M3, 
                y = c(list_s$M4, list_s$M5, list_s$M6, list_s$NCI_55, list_s$NCI_118, list_s$S1, list_s$S2))

