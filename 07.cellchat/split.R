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


list_ns <- split_list[c("MN1", "MN2", "MN3", "MN4", "FN1", "FN2", "FN3", "FN4")]

list_s <- split_list[c("MS1", "MS2", "MS3", "MS4", "FS1", "FS2", "FS3", "FS4")]


lung_ns <- merge(list_ns$MN1, 
                 y = c(list_ns$MN2, list_ns$MN3, list_ns$MN4, list_ns$FN1, list_ns$FN2, list_ns$FN3, list_ns$FN4))
lung_s <- merge(list_s$MS1, 
                y = c(list_s$MS2, list_s$MS3, list_s$MS4, list_s$FS1, list_s$FS2, list_s$FS3, list_s$FS4))

