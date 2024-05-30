##############################################################################
# Script information                                                      
# Title: Clustering
# Author: Erping Long
# Date: 2022-01-06
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

pbmc_4 <- readRDS(file = "lung_merge_all.rds")

dim(pbmc_4)

table(pbmc_4@meta.data$orig.ident)

####Harmony on RNA data

DefaultAssay(pbmc_4) <- "RNA"

pbmc_4 <- SCTransform(pbmc_4, verbose = TRUE)

##Set default assay to SCT

DefaultAssay(pbmc_4) <- "SCT"

#dimensional reduction and clustering

pbmc_4 <- ScaleData(pbmc_4, verbose = TRUE)
pbmc_4 <- RunPCA(pbmc_4, npcs = 50, verbose = TRUE)
pbmc_4 <- RunUMAP(pbmc_4, reduction = "pca", dims = 1:50, reduction.name = 'umap.RNA', reduction.key = 'rnaUMAP_')
pbmc_4 <- FindNeighbors(pbmc_4, reduction = "pca", dims = 1:50)
pbmc_4 <- FindClusters(pbmc_4, resolution = 0.5)

p1 <- DimPlot(pbmc_4, reduction = "umap.RNA", group.by = "orig.ident")
p2 <- DimPlot(pbmc_4, reduction = "umap.RNA", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_SCT.pdf", p1+p2, width = 24, height = 12)

#correcting batch effect

pbmc_4 <- RunHarmony(object = pbmc_4, group.by.vars = 'orig.ident',
reduction = 'pca',
assay.use = 'SCT',
reduction.save = "harmony_SCT",
project.dim = FALSE)

pbmc_4 <- RunUMAP(pbmc_4, dims = 1:50, reduction = 'harmony_SCT', reduction.name = "umap.SCT_harm", reduction.key = "sctUMAPharm_")
pbmc_4 <- FindNeighbors(pbmc_4, reduction = "harmony_SCT", dims = 1:50)
pbmc_4 <- FindClusters(pbmc_4, resolution = 0.5)

p3 <- DimPlot(pbmc_4, reduction = "umap.SCT_harm", group.by = "orig.ident")
p4 <- DimPlot(pbmc_4, reduction = "umap.SCT_harm", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_SCT_Har.pdf", p3+p4, width = 24, height = 12)

####Harmony on ATAC data

#dimensional reduction and clustering

DefaultAssay(pbmc_4) <- "ATAC"

pbmc_4 <- RunTFIDF(pbmc_4)
pbmc_4 <- FindTopFeatures(pbmc_4, min.cutoff = 'q0')
pbmc_4 <- RunSVD(pbmc_4)
pbmc_4 <- RunUMAP(pbmc_4, reduction = 'lsi', dims = 2:50, reduction.name = "umap.ATAC", reduction.key = "atacUMAP_")

pbmc_4 <- FindNeighbors(pbmc_4, reduction = "lsi", dims = 2:50)
pbmc_4 <- FindClusters(pbmc_4, resolution = 0.5)

p5 <- DimPlot(pbmc_4, reduction = "umap.ATAC", group.by = "orig.ident")
p6 <- DimPlot(pbmc_4, reduction = "umap.ATAC", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_ATAC.pdf", p5+p6, width = 24, height = 12)

#correcting batch effect

pbmc_4 <- RunHarmony(object = pbmc_4, group.by.vars = 'orig.ident',
reduction = 'lsi',
assay.use = 'ATAC',
reduction.save = "harmony_ATAC",
project.dim = FALSE)

pbmc_4 <- RunUMAP(pbmc_4, dims = 2:50, reduction = 'harmony_ATAC', reduction.name = "umap.ATAC_harm", reduction.key = "atacUMAPharm_")
pbmc_4 <- FindNeighbors(pbmc_4, reduction = "harmony_ATAC", dims = 2:50)
pbmc_4 <- FindClusters(pbmc_4, resolution = 0.5)

p7 <- DimPlot(pbmc_4, reduction = "umap.ATAC_harm", group.by = "orig.ident")
p8 <- DimPlot(pbmc_4, reduction = "umap.ATAC_harm", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_ATAC_Har.pdf", p7+p8, width = 24, height = 12)


####WNN on two modalities after harmony

pbmc_4 <- FindMultiModalNeighbors(pbmc_4, reduction.list = list("harmony_SCT", "harmony_ATAC"), dims.list = list(1:50, 2:50))
pbmc_4 <- RunUMAP(pbmc_4, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

#final clustering

pbmc_4 <- FindClusters(pbmc_4, graph.name = "wsnn", algorithm = 3, resolution = 0.5)

p9 <- DimPlot(pbmc_4, reduction = "wnn.umap", group.by = "orig.ident")
p10 <- DimPlot(pbmc_4, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_WNN_Har.pdf", p9+p10, width = 24, height = 12)

saveRDS(pbmc_4, file = "pbmc_integrate_all.rds")



