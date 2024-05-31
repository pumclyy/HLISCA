# cCRE-module, cCRE-cCRE, and cCRE-gene correlation

The co-accessible cCRE modules of two or more cCREs were identified by Cicero with
Louvain community detection algorithm and co-accessible score cutoff of 0.32 (automatically defined). A more stringent co-accessible score cutoff of 0.5 was used to define “directly co-accessible” cCREs. **(Cicero.R)**

For cCRE-gene correlation, we identified cCREs that may regulate a given gene by computing the correlation between gene expression and accessibility at nearby cCREs. Specifically,we performed the cCRE-gene Pearson correlation analysis by running the ‘LinkPeaks’ function with ‘distance’ of 1e+6 (+/- 1Mb of TSS), ‘min.cells’ of 10, ‘pvalue_cutoff’ of 0.05, and ‘score_cutoff’ of 0.05. Additionally, cCRE-gene correlation analysis was performed using larger distance settings of +/- 2Mb and +/-5Mb to provide insights to potential longer-range correlation. **(Linkage.R)**
