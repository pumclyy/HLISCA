# Peak calling, annotation and gene activity score

Peaks were called using MACS2 with default parameters after combining the reads of all the cells in each cell type to determine the genomic regions enriched for Tn5 accessibility from snATAC fragments. Term frequency inverse document frequency (TF-IDF) normalization on a matrix was performed using ‘RunTFIDF’ function. **(peak_calling.R)**

Peaks were then annotated according to distance to protein-coding genes using ChIPseeker, with “tssRegion” setting of -3000 to 3000. **(ChIPseeker.R)**

The gene activity score was calculated via the "GeneActivity" function of Signac. Then we calculated the correlation between percentages of cells expressing a marker gene and those with a corresponding gene activity across the cell types, among 41 canonical cell-type marker genes. **(GeneActivity.R)**
