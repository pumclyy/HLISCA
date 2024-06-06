##############################################################################
# Script information                                                      
# Title: DAR
# Author: Bolun Li
# Date: 2024-06-01
# Description: None
##############################################################################

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(GenomicRanges)
library(SeuratData)
library(ChIPseeker)

####ATAC analysis
seurat_object <- readRDS("lung.rds")
seurat_object = RenameIdents(seurat_object, '0' = 'NK', '1' = 'Lymphatic', '2' = 'AT2', '3' = 'AT1', '4' = 'T',
                          '5' = 'Club', '6' = 'Ciliated', '7' = 'Artery', '8' = 'Vein', '9' = 'Macrophage', '10' = 'Goblet',
                          '11' = 'Fibroblast', '12' = 'Monocyte', '13' = 'Basal', '14' = 'SMC', '15' = 'Capillary',
                          '16' = 'Mesothelial', '17' = 'B', '18' = 'Dendritic', '19' = 'AT1_AT2', '20' = 'NK_T',
                          '21' = 'MyoFib', '22' = 'AT2_pro')

seurat_object$CellType <- Idents(seurat_object)
seurat_object$CellType <- factor(seurat_object$CellType, levels = c('AT1', 'AT2', 'AT1_AT2','AT2_pro',
                                                              'Club','Ciliated','Goblet','Basal',
                                                              'Artery','Vein','Capillary', "Lymphatic",
                                                              'SMC','MyoFib','Fibroblast','Mesothelial',
                                                              'Macrophage','Monocyte','Dendritic','NK','NK_T', 'B', 'T'))

DefaultAssay(seurat_object) <- "ATAC"

#dividing each cell type into smoker and non_smoker clusters

Idents(seurat_object) <- "orig.ident"
seurat_object = RenameIdents(seurat_object, 'M1' = 'non_smoker', 'M3' = 'smoker', 'M4' = 'smoker', 'M5' = 'smoker', 'M6' = 'smoker',
                          'M7' = 'non_smoker', 'M8' = 'non_smoker', 'M9' = 'non_smoker', 'N17_20' = 'non_smoker', 'NCI_115' = 'non_smoker', 
                          'NCI_118' = 'smoker', 'NCI_13' = 'non_smoker', 'NCI_28' = 'non_smoker', 'NCI_55' = 'smoker', 
                          'S1' = 'smoker', 'S2' = 'smoker')
seurat_object$Smoking <- Idents(seurat_object)
seurat_object$Celltype_ss <- paste(seurat_object$CellType, seurat_object$Smoking, sep = "_")
table(seurat_object$Celltype_ss)

###calling peaks using divided 46 clusters

# call peaks using MACS2
peaks <- CallPeaks(seurat_object, group.by = "Celltype_ss",
                   macs2.path = "conda/envs/PeakCalling_analysis/bin/macs2",
                   outdir ="peaks")
peaks

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

saveRDS(peaks, file = "peaks.rds")

# quantify counts in each peak
peaks <- readRDS(file = "peaks.rds")
macs2_counts <- FeatureMatrix(
  fragments = Fragments(seurat_object),
  features = peaks,
  cells = colnames(seurat_object)
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# create a new assay using the MACS2 peak set and add it to the Seurat object
seurat_object[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c(":", "-"),
  fragments = Fragments(seurat_object),
  annotation = annotations
)

DefaultAssay(seurat_object) <- "peaks"

seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 5)
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- RunSVD(seurat_object)

saveRDS(seurat_object, file = "peak_called_by_cell_types_smoking.rds")
seurat_object <- readRDS("peak_called_by_cell_types_smoking.rds")


# run FindMarkers by parallel
library(future)
availableCores()
plan("multisession", workers = 16)

# Signac DAR
celltypes <- levels(seurat_object$CellType)
DefaultAssay(seurat_object)
table(Idents(seurat_object))
Idents(seurat_object) <- "Celltype_ss"
DAR.list <- list()
for (i in celltypes){
  cluster1 <- paste0(i,"_smoker")
  cluster2 <- paste0(i,"_non_smoker")
  marker_i <- FindMarkers(seurat_object, ident.1 = cluster1, ident.2 = cluster2, 
                          test.use = 'LR',latent.vars = 'nCount_peaks')
  DAR.list[[i]] <- marker_i
  filename <- paste0(i,".smoking.DAR.LR.csv")
  write.csv(marker_i,filename)
}

# select significant DAR with adjusted p value < 0.05
DAR.list.sig <- lapply(DAR.list, function(x) {
  x$peak <- rownames(x)
  x <- subset(x, p_val_adj < 0.05 )
})


# add annotation info for DAR
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(org.Hs.eg.db)
peaks <- seurat_object@assays$peaks@ranges
peakAnno <- annotatePeak(peaks,
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

write.csv(peakAnno,"project/AT2_co_CS/peaks_annotation_by_ChIPseeker.csv")
peaks_annotation_by_ChIPseeker <- read.csv("project/AT2_co_CS/peaks_annotation_by_ChIPseeker.csv")
peaks_annotation_by_ChIPseeker <- peaks_annotation_by_ChIPseeker[,-1]
rownames(peaks_annotation_by_ChIPseeker) <- paste(peaks_annotation_by_ChIPseeker$seqnames, 
                                                  peaks_annotation_by_ChIPseeker$start,
                                                  peaks_annotation_by_ChIPseeker$end, sep = "-")

DAR.list.new <-list()
celltype <- names(DAR.list.sig)
for (i in 1:23) {
  if (nrow(DAR.list.sig[[i]]) > 0) {
    markers <- DAR.list.sig[[i]]
    markers$cluster <- celltype[i]
    markers_annot <- cbind(markers, peaks_annotation_by_ChIPseeker[markers$peak,])
    DAR.list.new[[i]] <- markers_annot
  }else{
    DAR.list.new[[i]] <- NULL
  }
}

DAR.new <- Reduce(rbind, DAR.list.new)
write.table(DAR.new, file = "project/AT2_co_CS/DAR_LR_smoking_w_CT.txt",
            sep = "\t",quote = FALSE)
