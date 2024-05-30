#######sbatch --cpus-per-task=4 --mem=200g --time=6-00:00:00 run_ChIPseeker.sh
#######sinteractive --cpus-per-task=4 --mem=400g
##module load R/4.1

#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")

#BiocManager::install("ChIPseeker")

library(ChIPseeker)
#library(EnsDb.Hsapiens.v86)

#files <- getSampleFiles()
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(org.Hs.eg.db)

peakAnno <- annotatePeak("/data/Choi_lung/SHARE-seq/Seurat_SHARE/peaks_for_annotation.bed",
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

write.csv(peakAnno,"./peaks_annotation_by_ChIPseeker.csv")

