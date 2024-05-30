##############################################################################
# Script information                                                      
# Title: Gene activity score
# Author: Yuyan Li
# Date: 2024-04-12
# Description: None
##############################################################################

library(Seurat)
library(Signac)
library(harmony)
library(ggplot2)

pbmc_4 <- readRDS("16_multiome.rds")

#Calculate the Gene activity score

expscore <- GeneActivity(pbmc_4, assay = "peaks", features = rownames(pbmc_4))

pbmc_4@assays$score <- CreateAssayObject(counts = expscore)

#read the markers list

markers <- read.csv("markers.csv")

#calculate the correlation

for (i in markers) {
  p1 <- DotPlot(pbmc_4, features = i, assay = "SCT")
  p2 <- DotPlot(pbmc_4, features = i, assay = "score")
  
  cor.test(p1$data$pct.exp, p2$data$pct.exp)
}
