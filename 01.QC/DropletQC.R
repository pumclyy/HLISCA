## sinteractive --mem=200g --cpus-per-task=16

## module load R/4.0

## install.packages("devtools")
## devtools::install_github("powellgenomicslab/DropletQC", build_vignettes = TRUE)


library(Rsamtools)
library(GenomicRanges)
library(DropletQC)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
ID <- args[1]

# Create a reference to the BAM file
bam_file <- paste0("/data/Choi_lung/SHARE-seq/cellranger-arc/",ID,"/outs/gex_possorted_bam.bam")
umi_file <- read.csv(paste0("/data/Choi_lung/SHARE-seq/cellranger-arc/",ID,"/outs/per_barcode_metrics.csv"))

nf2 <- nuclear_fraction_tags(
   bam = bam_file,
   barcodes = paste0("/data/Choi_lung/SHARE-seq/cellranger-arc/",ID,"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),
   tiles = 1,
   cores = 1,
   verbose = FALSE
)
head(nf2)

write.csv(nf2,paste0("./",ID,"_NF.csv"))

#nf2_umi_file <- read.csv("/gpfs/gsfs12/users/Choi_lung/SHARE-seq/DropletQC/M2_NF2_UMI.csv")


# Get data frame with the nuclear fraction in the first column and umi counts in
# the second
#gbm <- filter(nf2, sample=="GBM")
gbm.nf.umi <- data.frame(nf=nf2_umi_file$nuclear_fraction,
                         umi=nf2_umi_file$gex_umis_count)

# Run identify_empty_drops
gbm.ed <- identify_empty_drops(nf_umi=nf2)
head(gbm.ed)

#table(gbm.ed$cell_status)

gbm.ed <- identify_empty_drops(nf_umi=gbm.nf.umi, include_plot = TRUE)

#gbm.ed <- identify_empty_drops(nf_umi=gbm.nf.umi,
 #                              nf_rescue = 0.02,umi_rescue = 1000,
  #                             include_plot = TRUE)


#gbm.ed$cell_type <- gbm$cell_type
#head(gbm.ed)

#pdf("DropletQC_plot.pdf")
#plot(gbm.ed)
#dev.off()


