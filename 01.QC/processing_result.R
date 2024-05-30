##############################################################################
# Script information                                                      
# Title: processing QC result
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

# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5(paste0("/data/Choi_lung/SHARE-seq/cellranger-arc/",ID,"/outs/filtered_feature_bc_matrix.h5"))

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
lung <- CreateSeuratObject(counts = rna_counts)
lung[["percent.mt"]] <- PercentageFeatureSet(lung, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- paste0("/data/Choi_lung/SHARE-seq/cellranger-arc/",ID,"/outs/atac_fragments.tsv.gz")

chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
lung[["ATAC"]] <- chrom_assay


#load the result of droplet and doublet
DropletQC_scrublet <- read.csv("scrublet.csv")

lung$nuclear_fraction <- DropletQC_scrublet$nuclear_fraction

lung$scrublet <- DropletQC_scrublet$scrublet


#filter by the result of droplet and doublet
lung_nf <- subset(
  x = lung,
  subset = nuclear_fraction > 0.25
    )

lung_sc <- subset(
  x = lung_nf,
  scrublet == FALSE
    )


#filter by "nCount_ATAC", "nCount_RNA" and "percent.mt"
plot1 <- VlnPlot(lung_sc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

lung_2 <- subset(
  x = lung_sc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 1e3
    )

lung_3 <- subset(
  x = lung_2,
  subset = nCount_RNA < 25000 &
    nCount_RNA > 1000
)

lung_4 <- subset(
  x = lung_3,
  subset = percent.mt < 10
)

# RNA analysis and SCTransform
DefaultAssay(lung_4) <- "RNA"
lung_4 <- SCTransform(lung_4, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

saveRDS(lung_4, file = "lung.rds")


