# Co-embedding analysis

The co-embedding analysis involved three main steps. First, a set of anchors were identified between RNA (reference) and gene activity score (query) via the “FindTransferAnchors” function of Seurat v.4.0.6. Second, based on the identified anchors we imputed RNA expression into the snATAC-seq cells via the “TransferData” functions, and then merged the two datasets. Third, the PCA and UMAP were perfomed on the merged dataset to visualize the co-embedding dimensions. Top 2,000 variable genes were used in the co-embedding analysis. **(co-embedding.R)**
