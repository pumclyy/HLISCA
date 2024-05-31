#!/usr/bin/env Rscript
library(motifbreakR)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(stringr)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(MotifDb)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)

#gwasfile <- args[1]
#qtlfile <- args[2]
#prefix <- args[3]
#tissue <- args[4]
#args <- commandArgs(TRUE)


snp <- args[1]

snp.plot = c("rs3769823")
snp.plot

#rawlist <- scan(snpfile, what = character())
#rawlist

#snps.frombed <- snps.from.file(file = args[1],dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,search.genome = BSgenome.Hsapiens.UCSC.hg19,format = "bed")

snps.fromlist <- snps.from.rsid(rsid = snp.plot,dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,search.genome = BSgenome.Hsapiens.UCSC.hg38)
snps.fromlist

data("hocomoco")
results <- motifbreakR(snpList = snps.fromlist, filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
					   BPPARAM = BiocParallel::bpparam())

results <- calculatePvalue(results)
results

pdf("rs3769823_motif_IRF8_20230811.pdf", height = 10)
plotMB(results = results, rsid = snp.plot, effect = "strong")
dev.off()


#write.csv(results,paste0(prefix,"_motifBreakR_results.csv"))



#results2 <- mcols(results) %>% as.data.frame() %>% tbl_df() %>% mutate(seqnames= as.character(seqnames(results))) %>% select(seqnames,REF:effect)%>% mutate(dscore=scoreAlt-scoreRef,dpct=pctAlt-pctRef,dscore2=abs(dscore)) %>% filter(effect=="strong") %>%  arrange(dscore) %>% mutate(mutation=paste0(seqnames,":",snpPos,":",REF,":", ALT))
#results2 <- unique(results2)
#write_csv(results2,paste0(args[2],"_strong.csv"))

#data1 <- results2 %>% arrange(desc(dscore2)) %>% filter(dscore<0) %>% group_by(seqnames,REF, ALT, snpPos)  %>% slice(1) %>% ungroup() %>% arrange(desc(dscore2))
#data2 <- results2 %>% arrange(desc(dscore2)) %>% filter(dscore>0) %>% group_by(seqnames,REF, ALT, snpPos)  %>% slice(1) %>% ungroup() %>% arrange(desc(dscore2))
#datatmp <- full_join(data1 %>% select(mutation,providerId,dscore), data2 %>% select(mutation,providerId,dscore),by=c("mutation"="mutation")) %>% mutate(effect=ifelse((!is.na(dscore.x) & abs(dscore.x)>=abs(dscore.y))|is.na(dscore.y),"Disrupte motif","Creat motif"),motif=ifelse((!is.na(dscore.x) & abs(dscore.x)>=abs(dscore.y))|is.na(dscore.y),providerId.x,providerId.y),dscore=ifelse((!is.na(dscore.x) & abs(dscore.x)>=abs(dscore.y))|is.na(dscore.y),abs(dscore.x),abs(dscore.y))) %>% mutate(motif=str_replace_all(motif,"_HUMAN","")) 
#motiftmp <- datatmp %>% group_by(motif) %>% tally() %>% arrange(desc(n)) %>% filter(n>1) %>% .[[1]]
#datatmp <- datatmp %>% mutate(motif2=ifelse(motif %in% motiftmp,motif,"Others"))
#write_csv(datatmp,paste0(args[2],"_top.csv"))


#test
pca.snps.file <- system.file("extdata", "pca.enhancer.snps", package = "motifbreakR")
pca.snps <- as.character(read.table(pca.snps.file)[,1])

variants <- snps.from.rsid(rsid = pca.snps,
                           dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                           search.genome = BSgenome.Hsapiens.UCSC.hg19)
motifbreakr.results <- motifbreakR(snpList = variants, pwmList = MotifDb, threshold = 0.9)
plotMB(results = motifbreakr.results, rsid = "rs7837328", effect = "strong")



