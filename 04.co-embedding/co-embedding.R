##############################################################################
# Script information                                                      
# Title: Co-embedding
# Author: Yuyan Li
# Date: 2024-04-11
# Description: None
##############################################################################

library(SeuratData)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)

lung_4 <- readRDS("lung.rds")

lung.rna <- DietSeurat(lung_4, assay = "RNA")
lung.atac <- DietSeurat(lung_4, assays = "ATAC")

# Perform standard analysis of each modality independently RNA analysis

lung.rna <- NormalizeData(lung.rna)
lung.rna <- FindVariableFeatures(lung.rna)
lung.rna <- ScaleData(lung.rna)
lung.rna <- RunPCA(lung.rna)
lung.rna <- RunUMAP(lung.rna, dims = 1:30)

# ATAC analysis add gene annotation information

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(lung.atac) <- annotations

# We exclude the first dimension as this is typically correlated with sequencing depth

lung.atac <- RunTFIDF(lung.atac)
lung.atac <- FindTopFeatures(lung.atac, min.cutoff = "q0")
lung.atac <- RunSVD(lung.atac)
lung.atac <- RunUMAP(lung.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# quantify gene activity

gene.activities <- GeneActivity(lung.atac, features = VariableFeatures(lung.rna))

# add gene activities as a new assay

lung.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities

DefaultAssay(lung.atac) <- "ACTIVITY"
lung.atac <- NormalizeData(lung.atac)
lung.atac <- ScaleData(lung.atac, features = rownames(lung.atac))

# Identify anchors

transfer.anchors <- FindTransferAnchors(reference = lung.rna, query = lung.atac, features = VariableFeatures(object = lung.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

genes.use <- VariableFeatures(lung.rna)
refdata <- GetAssayData(lung.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = lung.atac[["lsi"]],
                           dims = 2:30)

lung.atac[["RNA"]] <- imputation

lung.rna$or <- rep("rna", 117911)
lung.atac$or <- rep("atac", 117911)
coembed <- merge(x = lung.rna, y = lung.atac, add.cell.ids = c("rna", "atac"))

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("or", "final_cluster"))

ggsave("coembeding.pdf", width = 24, height = 12)
