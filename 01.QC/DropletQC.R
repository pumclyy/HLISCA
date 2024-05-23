##############################################################################
# Script information                                                      
# Title: Droplet removing
# Author: Erping Long
# Date: 2021-03-15
# Description: None
##############################################################################

#import packages
library(Rsamtools)
library(GenomicRanges)
library(DropletQC)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)

bam_file <- gex_possorted_bam.bam
umi_file <- read.csv("per_barcode_metrics.csv")

# Calculating the nuclear fraction score

nf <- nuclear_fraction_tags(
   bam = bam_file,
   barcodes = "barcodes.tsv.gz",
   tiles = 1,
   cores = 1,
   verbose = FALSE
)

write.csv(nf , "NF.csv")


gbm.nf.umi <- data.frame(nf=nf$nuclear_fraction,
                         umi=umi_file$gex_umis_count)

# Identifying empty droplets

gbm.ed <- identify_empty_drops(nf_umi=gbm.nf.umi, include_plot = TRUE)

