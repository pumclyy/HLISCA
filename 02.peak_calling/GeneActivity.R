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

lung <- readRDS("16_multiome.rds")

#Calculate the Gene activity score

expscore <- GeneActivity(lung, assay = "peaks", features = rownames(lung))

lung@assays$score <- CreateAssayObject(counts = expscore)

#read the markers list

markers <- read.csv("markers.csv")

#calculate the correlation

for (i in markers) {
  p1 <- DotPlot(lung, features = i, assay = "SCT")
  p2 <- DotPlot(lung, features = i, assay = "score")
  
  cor.test(p1$data$pct.exp, p2$data$pct.exp)
}
