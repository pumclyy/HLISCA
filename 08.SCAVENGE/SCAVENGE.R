##############################################################################
# Script information                                                      
# Title: Trait relevance score calculation
# Author: Erping Long
# Date: 2023-5-17
# Description: None
##############################################################################

library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)


trait_file <- "Overall_PP_SHARE_Flat.bed"

pbmc5krda <- "NCI_M1_SE_gvar.rda"

NCI_M1 <- load(pbmc5krda)

###gchromVAR analysis

SE_pbmc5k <- addGCBias(SE_gvar, genome = BSgenome.Hsapiens.UCSC.hg38)
SE_pbmc5k_bg <- getBackgroundPeaks(SE_pbmc5k, niterations=200)
trait_import <- importBedScore(rowRanges(SE_pbmc5k), trait_file, colidx=5)
SE_pbmc5k_DEV <- computeWeightedDeviations(SE_pbmc5k, trait_import, background_peaks = SE_pbmc5k_bg)

###Reformat results
z_score_mat <- data.frame(colData(SE_pbmc5k), z_score=t(assays(SE_pbmc5k_DEV)[["z"]]) %>% c)
head(z_score_mat)

### Generate the seed cell index (using the top 5% if too many cells are eligible)
seed_idx <- seedindex(z_score_mat$z_score, 0.05)

###calculate scale factor
scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)

###Construct m-knn graph
peak_by_cell_mat <- assay(SE_pbmc5k)
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)

###Calculate lsi-mat
lsi_mat <- do_lsi(tfidf_mat, dims=30)

###Calculate m-knn graph
mutualknn30 <- getmutualknn(lsi_mat, 30)

###Network propagation
np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)

###Trait relevant score (TRS) with scaled and normalized
omit_idx <- np_score==0
mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
TRS <- TRS * scale_factor
mono_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)

write.csv(mono_mat,"results.csv")

###UMAP plots of cell type annotation and cell-to-cell graph 
p <- ggplot(data=mono_mat, aes(x, y, color=color)) + geom_point(size=1, na.rm = TRUE) + 
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p

###Comparsion before and after SCAVENGE analysis

###scatter plot
p1 <- ggplot(data=mono_mat, aes(x, y, color=z_score)) + geom_point(size=1, na.rm = TRUE, alpha = 0.6) + 
scale_color_gradientn(colors = viridis) + scale_alpha()+
    pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p1

###bar plot
pp1 <- ggplot(data=mono_mat,  aes(x=color, y=z_score))  +
    geom_boxplot(aes(fill=color, color=color), outlier.shape=NA) + 
    guides(fill=FALSE) + pretty_plot(fontsize = 10) +
    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp1

### Trait relevant cell determination from permutation test 
mono_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, seed_idx=mono_mat$seed_idx, topseed_npscore=mono_mat$np_score, permutation_times=1000, true_cell_significance=0.05, rda_output=F, mycores=8, rw_gamma=0.05)
mono_mat2 <- data.frame(mono_mat, mono_permu)


### Look at the distribution of statistically significant phenotypically enriched and depleted cells
##Enriched cells

mono_mat2 %>%
    group_by(color) %>% 
        summarise(enriched_cell=sum(true_cell_top_idx)) %>% 
            ggplot(aes(x=color, y=enriched_cell, fill=color)) + geom_bar(stat="identity") + theme_classic()

##Depleted cells

mono_mat2$rev_true_cell_top_idx <- !mono_mat2$true_cell_top_idx
mono_mat2 %>%
    group_by(color) %>% 
        summarise(depleted_cell=sum(rev_true_cell_top_idx)) %>% 
            ggplot(aes(x=color, y=depleted_cell, fill=color)) + geom_bar(stat="identity") + theme_classic()

