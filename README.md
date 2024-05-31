# Integrated atlas of human lung across chromatin accessibilty and transcriptome

This repository contains the code to analyze the single cell multiome data of human lung as part of the manuscript "Context-aware single-cell multiome approach identified cell-type specific lung cancer susceptibility genes".

---

## data

The data we used is available on [Zenodo](https://)

---

## pipeline

### Quality control and filtering

We used DropletQC(v.0.9) to remove “empty” droplets containing ambient RNA from the gene expression matrices. The resulting expression matrices were processed individually in R (v.4.1.3) using Seurat (v.4.0.6) and Gencode v.27 for gene identification. Filtered gene–barcode matrices were normalized with the ‘SCTransform’ function of Seurat, and the top 2,000 variable genes were identified. Gene expression matrices were scaled and centered using the ‘ScaleData’ function. See details in [01.QC](https://github.com/pumclyy/16_multiome/tree/main/01.QC).

### Peak Calling

The snATAC peak calling and annotation were performed following the Signac pipeline. Specifically, peaks were called using MACS2 with default parameters after combining the reads of all the cells in each cell type to determine the genomic regions enriched for Tn5 accessibility from snATAC fragments, resulting in 330,453 peaks in total. Peaks were then annotated according to distance to protein-coding genes using ChIPseeker. The  gene activity score was calculated via the "GeneActivity" function of Signac. See details in [02.peak_calling](https://github.com/pumclyy/16_multiome/tree/main/02.peak_calling).

### Clustering

Using the normalized gene expression data, we performed principal component analysis
(PCA) with 50 PCs to compute and store. A uniform manifold approximation and projection (UMAP)-based approach was applied for expression matrices with the first 50 PCs and for chromatin accessibility matrices with the 2nd through 50th PCs (the first PC was excluded as this is typically correlated with sequencing depth). Both expression and chromatin accessibility matrices ere corrected for batch effect using Harmony. A Weighted Nearest Neighbor (WNN) method was applied to integrate the weighted combination of RNA and ATAC-seq modalities. The ‘FindClusters’ function was applied for clustering using smart local moving (SLM) algorithm for modularity optimization at a resolution of 0.5. See details in [03.clustering](https://github.com/pumclyy/16_multiome/tree/main/03.clustering)

### Co-embedding analysis

We performed co-embedding analysis between snRNA and snATAC modalities following the Seurat pipeline. See details in [04.co-embedding](https://github.com/pumclyy/16_multiome/tree/main/04.co-embedding)

### Analysis of smoking-responsive genes.

Smoking-responsive genes for each cell type were inferred by pseudobulk differential gene expression analysis using DESeq2 (v1.41.1). Gender was incorporated as the covariate into the model. See details in [05.DEG](https://github.com/pumclyy/16_multiome/tree/main/05.DEG)

### Analysis of smoking-responsive cCREs

### cellchat

CellChat (version 1.6.0) was used to infer ligand–receptor interactions based on scRNA-seq data. The intercellular communication analyses were performed separately using cells from ever-smokers and never-smokers for comparison at different levels. Then the integrated Human Lung Cell Atlas (HLCA) core dataset was used to validate the trend of MHC-I and MHC-II communication. See details in [07.cellchat](https://github.com/pumclyy/16_multiome/tree/main/07.cellchat)

### Trait relevance score calculation

We inferred the lung-cancer-associated score for each cell using SCAVENGE pipeline, based on snATAC-seq and lung cancer GWAS data. See details in [08.SCAVENGE](https://github.com/pumclyy/16_multiome/tree/main/08.SCAVENGE)

### Allelic effects of predicted TF binding, cell-type-specific TF abundance assessment and TF footprinting

Prediction of variant effects on TF binding sites was performed with the motifbreakR package. "Abundant TF" was identified by expression level and expression percent of TF. TF footprint analysis was performed for each allelic-binding TF using the ‘Footprint’ function in Signac by restricting to the peak regions. See details in [09.TF_analysis](https://github.com/pumclyy/16_multiome/tree/main/09.TF_analysis)

### Reporter assays for variant allelic function test

Allelic transcriptional activity of CCVs from the 15_5p15.33 locus was assessed as part of MPRA. See details in [10.MPRA](https://github.com/pumclyy/16_multiome/tree/main/10.MPRA)