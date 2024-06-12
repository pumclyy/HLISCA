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

lung_4 <- readRDS(file = "lung_merge_all.rds")

dim(lung_4)

table(lung_4@meta.data$orig.ident)

####Harmony on RNA data

DefaultAssay(lung_4) <- "RNA"

lung_4 <- SCTransform(lung_4, verbose = TRUE)

##Set default assay to SCT

DefaultAssay(lung_4) <- "SCT"

#dimensional reduction and clustering

lung_4 <- ScaleData(lung_4, verbose = TRUE)
lung_4 <- RunPCA(lung_4, npcs = 50, verbose = TRUE)
lung_4 <- RunUMAP(lung_4, reduction = "pca", dims = 1:50, reduction.name = 'umap.RNA', reduction.key = 'rnaUMAP_')
lung_4 <- FindNeighbors(lung_4, reduction = "pca", dims = 1:50)
lung_4 <- FindClusters(lung_4, resolution = 0.5)

p1 <- DimPlot(lung_4, reduction = "umap.RNA", group.by = "orig.ident")
p2 <- DimPlot(lung_4, reduction = "umap.RNA", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_SCT.pdf", p1+p2, width = 24, height = 12)

#correcting batch effect

lung_4 <- RunHarmony(object = lung_4, group.by.vars = 'orig.ident',
reduction = 'pca',
assay.use = 'SCT',
reduction.save = "harmony_SCT",
project.dim = FALSE)

lung_4 <- RunUMAP(lung_4, dims = 1:50, reduction = 'harmony_SCT', reduction.name = "umap.SCT_harm", reduction.key = "sctUMAPharm_")
lung_4 <- FindNeighbors(lung_4, reduction = "harmony_SCT", dims = 1:50)
lung_4 <- FindClusters(lung_4, resolution = 0.5)

p3 <- DimPlot(lung_4, reduction = "umap.SCT_harm", group.by = "orig.ident")
p4 <- DimPlot(lung_4, reduction = "umap.SCT_harm", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_SCT_Har.pdf", p3+p4, width = 24, height = 12)

####Harmony on ATAC data

#dimensional reduction and clustering

DefaultAssay(lung_4) <- "ATAC"

lung_4 <- RunTFIDF(lung_4)
lung_4 <- FindTopFeatures(lung_4, min.cutoff = 'q0')
lung_4 <- RunSVD(lung_4)
lung_4 <- RunUMAP(lung_4, reduction = 'lsi', dims = 2:50, reduction.name = "umap.ATAC", reduction.key = "atacUMAP_")

lung_4 <- FindNeighbors(lung_4, reduction = "lsi", dims = 2:50)
lung_4 <- FindClusters(lung_4, resolution = 0.5)

p5 <- DimPlot(lung_4, reduction = "umap.ATAC", group.by = "orig.ident")
p6 <- DimPlot(lung_4, reduction = "umap.ATAC", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_ATAC.pdf", p5+p6, width = 24, height = 12)

#correcting batch effect

lung_4 <- RunHarmony(object = lung_4, group.by.vars = 'orig.ident',
reduction = 'lsi',
assay.use = 'ATAC',
reduction.save = "harmony_ATAC",
project.dim = FALSE)

lung_4 <- RunUMAP(lung_4, dims = 2:50, reduction = 'harmony_ATAC', reduction.name = "umap.ATAC_harm", reduction.key = "atacUMAPharm_")
lung_4 <- FindNeighbors(lung_4, reduction = "harmony_ATAC", dims = 2:50)
lung_4 <- FindClusters(lung_4, resolution = 0.5)

p7 <- DimPlot(lung_4, reduction = "umap.ATAC_harm", group.by = "orig.ident")
p8 <- DimPlot(lung_4, reduction = "umap.ATAC_harm", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_ATAC_Har.pdf", p7+p8, width = 24, height = 12)


####WNN on two modalities after harmony

lung_4 <- FindMultiModalNeighbors(lung_4, reduction.list = list("harmony_SCT", "harmony_ATAC"), dims.list = list(1:50, 2:50))
lung_4 <- RunUMAP(lung_4, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

#final clustering

lung_4 <- FindClusters(lung_4, graph.name = "wsnn", algorithm = 3, resolution = 0.5)

p9 <- DimPlot(lung_4, reduction = "wnn.umap", group.by = "orig.ident")
p10 <- DimPlot(lung_4, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_WNN_Har.pdf", p9+p10, width = 24, height = 12)

saveRDS(lung_4, file = "lung_integrate_all.rds")



