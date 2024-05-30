##############################################################################
# Script information                                                      
# Title: Cellchat
# Author: Yuyan Li
# Date: 2023-05-10
# Description: None
##############################################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(celldex)
library(cowplot)
library(tidyverse)
library(CellChat)

#create cellchat object

cellchat_ns <- createCellChat(lung_ns, meta = lung_ns@meta.data, group.by = "type")
cellchat_s <- createCellChat(lung_s, meta = lung_s@meta.data, group.by = "type")
cellchat_ns@DB <- CellChatDB.human
cellchat_s@DB <- CellChatDB.human

#calculate the communication probability

cellchat_ns <- subsetData(cellchat_ns)
cellchat_ns <- identifyOverExpressedGenes(cellchat_ns)
cellchat_ns <- identifyOverExpressedInteractions(cellchat_ns)
cellchat_ns <- computeCommunProb(cellchat_ns)
cellchat_ns <- computeCommunProbPathway(cellchat_ns)
cellchat_ns <- aggregateNet(cellchat_ns)

cellchat_s <- subsetData(cellchat_s)
cellchat_s <- identifyOverExpressedGenes(cellchat_s)
cellchat_s <- identifyOverExpressedInteractions(cellchat_s)
cellchat_s <- computeCommunProb(cellchat_s)
cellchat_s <- computeCommunProbPathway(cellchat_s)
cellchat_s <- aggregateNet(cellchat_s)

cellchat_ns <- netAnalysis_computeCentrality(object = cellchat_ns, slot.name = "netP")
cellchat_s <- netAnalysis_computeCentrality(object = cellchat_s, slot.name = "netP")
object.list <- list(ns = cellchat_ns, s = cellchat_s)
cellchat_com <- mergeCellChat(object.list, add.names = names(object.list))

#save the probability

ns.netp <- subsetCommunication(cellchat_ns, slot.name = 'netP')
write.csv(ns.netp, "pathway_ns.csv")

s.netp <- subsetCommunication(cellchat_s, slot.name = 'netP')
write.csv(s.netp, "pathway_s.csv")

ns.net <- subsetCommunication(cellchat_ns, slot.name = 'net')
write.csv(ns.net, "L_R_ns.csv")

s.net <- subsetCommunication(cellchat_s, slot.name = 'net')
write.csv(s.net, "L_R_s.csv")
