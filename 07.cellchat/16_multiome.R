#引用
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(celldex)
library(cowplot)
library(tidyverse)
library(CellChat)
#读入文件
pbmc <- readRDS("/data/liyy/16_multiome.rds")
#生成cellchat对象并存储
pbmc$type <- pbmc$final_cluster #由于final_cluster中含有0，令加一行作为分类依据
levels(pbmc$type) <- c(1:23) #将因子替换为1到23
DefaultAssay(pbmc) <- "SCT" #更改默认矩阵,用于cellchat分析
cellchat <- createCellChat(object = pbmc, meta = pbmc@meta.data, group.by = "type") #生成cellchat对象
saveRDS(cellchat, file = "/data/liyy/16_multiome/cellchat_16multiome.rds") #存储cellchat对象，方便之后调用
cellchat <- readRDS("/data/liyy/16_multiome/cellchat_16multiome.rds") #读取数据
#设置配体受体数据库
cellchatDB <- CellChatDB.human #选择人类的细胞通讯数据库
cellchatDB.use <- cellchatDB #选择整个数据库用于比对
cellchat@DB <- cellchatDB.use #将用于比对的数据库存入对象中
#预处理数据
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#计算细胞通讯概率
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#可视化
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weight/strength")
#由于细胞数量太多，整体难以看出细节，因此对单个细胞可视化
#细胞数量太多，一幅图不能全部呈现，分四次呈现
mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:6) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
par(mfrow = c(2,3), xpd=TRUE)
for (i in 7:12) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
par(mfrow = c(2,3), xpd=TRUE)
for (i in 13:18) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
par(mfrow = c(2,3), xpd=TRUE)
for (i in 19:23) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#分析各信号通路在不同细胞间的概率
cellchat@netP$pathways #显示所有的信号通路
pathways.show <- c("CD46") #选取想要分析的信号通路
vertex.receiver = seq(1,12)
#以三种方式可视化
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#EPCAM--3,4,6,7,11,14,20,23   PTPRC--1,5,10,13,18,19,21
netVisual_bubble(cellchat, sources.use = c(3,4,6,7,11,14,20,23), targets.use = c(1,5,10,13,18,19,21), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,5,10,13,18,19,21), remove.isolate = FALSE)
