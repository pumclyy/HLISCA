# Trait relevance score calculation

We inferred the lung-cancer-associated score for each cell using SCAVENGE pipeline, based on snATAC-seq and lung cancer GWAS data.

For GWAS input, we used flat probability scores evenly divided across all the CCVs in each locus.

The peak-by-cell matrix of snATAC-seq data were processed using ArchR package. We processed snATAC-seq data using ArchR pipeline to obtain a peak-by-cell matrix. The peak-by-cell count matrix and corresponding meta data were extracted and stored in a RangedSummarizedExperimentobject. **(ArchR.R)**

gchromVAR was performed to calculate the original colocalization scores between GWAS and snATAC-seq. Seed cells proportion was set as 5%, and network propagation was applied to calculate the lung-cancer-associated score for each cell. **(SCAVENGE.R)**