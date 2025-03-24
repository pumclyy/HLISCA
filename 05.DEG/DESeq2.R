##############################################################################
# Script information                                                      
# Title: DESeq2
# Author: Yuyan Li
# Date: 2023-09-15
# Description: None
##############################################################################

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(readr)
library(harmony)
library(GenomicRanges)
library(SeuratData)
library(DESeq2)
library(SingleCellExperiment)
library(Matrix.utils)
library(magrittr)
library(purrr)
library(scater)
library(cowplot)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(ashr)
library(MAST)

lung <- readRDS("lung.rds")


#generate pseudobulk data


##creat single cell experiment object
counts <- lung@assays$RNA@counts
metadata <- lung@meta.data

metadata$cluster_id <- factor(lung$final_cluster)
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

sce$sample_id <- as.factor(sce$orig.ident)
groups <- colData(sce)[, c("cluster_id", "sample_id")]


##generate sample level metadata

kids <- purrr::set_names(levels(sce$cluster_id))
nk <- length(kids)
sids <- purrr::set_names(levels(sce$sample_id))
ns <- length(sids)

n_cells <- as.numeric(table(sce$sample_id))
m <- match(sids, sce$sample_id)
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")


##aggregate the counts per sample_id and cluster_id to generate psudobulk data 

groups <- colData(sce)[, c("cluster_id", "sample_id")]


pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "[MNS].*")))
class(pb)
str(pb)

options(width = 100)
table(sce$cluster_id, sce$sample_id)


##differential expression analysis

get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()



gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

ei$group_id <- factor(c("N", "S", "S", "S", "S", "N", "N", "N", "N", "N", "S", "N", "N", "S", "S", "S"))

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 

metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id) 


clusters <- levels(as.factor(metadata$cluster_id))

cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
rownames(cluster_metadata) <- cluster_metadata$sample_id

counts <- pb[[clusters[1]]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

all(rownames(cluster_metadata) == colnames(cluster_counts))

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id + sample_id)
dds <- DESeq(dds)
plotDispEsts(dds)

contrast <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[1])

res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

write.csv(res_tbl,
          paste0(clusters[1], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


