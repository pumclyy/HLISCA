##############################################################################
# Script information                                                      
# Title: MotifbreakR
# Author: Erping Long, Yuyan Li
# Date: 2023-5-17
# Description: None
##############################################################################

library(motifbreakR)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(stringr)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(MotifDb)
library(ggplot2)

snp <- read.csv("snp_list.csv")

snps.fromlist <- snps.from.rsid(rsid = snp,dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,search.genome = BSgenome.Hsapiens.UCSC.hg38)

data("hocomoco")
results <- motifbreakR(snpList = snps.fromlist, filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
					   BPPARAM = BiocParallel::bpparam())

results <- calculatePvalue(results)
results

write.csv(results,"motifBreakR_results.csv")



