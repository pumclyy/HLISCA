

# Integrated atlas of human lung across chromatin accessibilty and transcriptome

This repository contains the code to analyze the single cell multiome data of human lung as part of the manuscript "Context-aware single-cell multiome approach identified cell-type specific lung cancer susceptibility genes".

---

## data

The data we used is available on [Zenodo](https://)

---

## pipeline

### Quality control and filtering

We used DropletQC(v.0.9) to remove “empty” droplets containing ambient RNA from the gene expression matrices. The resulting expression matrices were processed individually in R (v.4.1.3) using Seurat (v.4.0.6) and Gencode v.27 for gene identification. Filtered gene–barcode matrices were normalized with the ‘SCTransform’ function of Seurat, and the top 2,000 variable genes were identified. Gene expression matrices were scaled and centered using the ‘ScaleData’ function. See details in [01.QC](https://)

### Peak Calling

