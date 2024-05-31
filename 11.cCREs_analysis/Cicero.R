##############################################################################
# Script information                                                      
# Title: Identifying cCRE-module
# Author: Erping Long
# Date: 2023-03-19
# Description: None
##############################################################################

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
library(monocle3)

lung_4 <- readRDS(file = "lung.rds")

DefaultAssay(lung_4) <- "peaks"

# convert to CellDataSet format and make the cicero object
# wnn.umap has been changed to WNN.UMAP

lung_4.cds <- as.cell_data_set(lung_4)
lung_4.cicero <- make_cicero_cds(lung_4.cds, reduced_coordinates = reducedDims(lung_4.cds)$WNN.UMAP) 

#use a customized hg38 file

hg38_genome <- read.csv("hg38_sequence_length.csv")

## Usually run with sample_num = 100 ##
conns <- run_cicero(lung_4.cicero, hg38_genome, sample_num = 100)
write.csv(conns,"./conns.csv")

ccans <- generate_ccans(conns)

write.csv(ccans,"./ccans.csv")

links <- ConnectionsToLinks(conns = conns, ccans = ccans)

write.csv(links,"./links.csv")
