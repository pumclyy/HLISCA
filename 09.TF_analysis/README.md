# Allelic effects of predicted TF binding, cell-type-specific TF abundance assessment and TF footprinting

Prediction of variant effects on TF binding sites was performed with the motifbreakR
package and a comprehensive collection of human TF-binding site models (HOCOMOCO, v11). We selected the information content algorithm and used a threshold of e-4 as the maximum p value for a TF-binding site match in motifbreakR. The allelic-binding effect was defined by the difference between alternative allele score and reference allele score larger than 0.7. **(motifbreakR.R)**

"Abundant TF" was identified by expression level and expression percent of TF. “Abundant TF” in a given cell type was defined as: 1) TF is expressed in > 50% of the cells in that cell type and 2) TF expression level in the same cell type is above 75 percentile of the values from all predicted TF-cell type pairs.

TF footprint analysis was performed for each allelic-binding TF using the ‘Footprint’ function in Signac by restricting to the peak regions. **(TF_footprint.R)**
