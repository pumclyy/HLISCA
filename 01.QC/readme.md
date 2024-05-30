# Quality control and filtering

We used DropletQC and scrublet to remove droplets and doublets respectively. First we performed DropletQC to indentify information on "empty" droplets. Similarly, we performed scrublet to get information of doublets. **(DropletQC.R and scrublet.py)**

Then, the resulting expression matrices were processed individually in R (v.4.1.3) using Seurat (v.4.0.6) and Gencode v.27 for gene identification. We integrated the information on "empty" droplets and doublets into Seurat object. We excluded cells with less than 500 or more than 25000 nCount_RNA (number of RNA read counts), cells with less than 1000 and more than 70,000 nCount_ATAC (number of ATAC read counts), and cells with more than 10% of counts corresponding to mitochondrial genes. And we removed doublets and "empty" droplets. After that, we normalized gene-barcode matrices with the ‘SCTransform’ function of Seurat.    **(processing_result.R)**

Finally, we merged the result of each sample together. **(Seurat_Merge.R)**
