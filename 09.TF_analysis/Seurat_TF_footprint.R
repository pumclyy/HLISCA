#######sbatch --cpus-per-task=4 --mem=200g --time=6-00:00:00 run_TF_Footprint.sh
#######sinteractive --cpus-per-task=4 --mem=400g
##module load R/4.1

##devtools::install_local('/gpfs/gsfs12/users/longe2/Software/harmony-master/')
##Harmony bugs, addressed in https://github.com/immunogenomics/harmony/issues/159

##setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
##install.packages("Signac")
##BiocManager::install("EnsDb.Hsapiens.v86")
##devtools::install_github('satijalab/seurat-data')
##BiocManager::install("JASPAR2020")
##BiocManager::install("TFBSTools")
##BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
##TFBSTool need root installation, sudo apt install libgsl-dev
##BiocManager::install("motifmatchr")

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

sessionInfo()

pbmc_peaks <- readRDS(file = "/NASdata/longep/pbmc_peaks_for_TF_footprint.rds")

DefaultAssay(pbmc_peaks) <- "peaks"

library(JASPAR2020)
library(TFBSTools)
#library(DirichletMultinomial)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)

#### TF Footprinting 
pbmc_peaks <- Footprint(
  object = pbmc_peaks,
  motif.name = c("EGR1","EGR2","EGR3","EGR4","EHF"),
  genome = BSgenome.Hsapiens.UCSC.hg38)

p2 <- PlotFootprint(pbmc_peaks, features = c("EGR1","EGR2","EGR3","EGR4","EHF"))

p2 + patchwork::plot_layout(ncol = 1)

ggsave("/NASdata/longep/TF_Footprint/footprint_EGR1_EGR2_EGR3_EGR4_EHF.pdf", 
       p2+ patchwork::plot_layout(ncol = 1), width = 10, height = 30)


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
pbmc_4 <- AddMotifs(
  object = pbmc_4,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)



##https://satijalab.org/seurat/reference/dietseurat


saveRDS(pbmc_peaks, file = "/NASdata/longep/pbmc_peaks_for_TF_footprint.rds")

#### TF Footprinting 
pbmc_peaks <- Footprint(
  object = pbmc_peaks,
  motif.name = c("ZNF148"),
  genome = BSgenome.Hsapiens.UCSC.hg38)

p3 <- PlotFootprint(pbmc_peaks, features = c("ZNF148"))

ggsave("/NASdata/longep/TF_Footprint/footprint_ZNF148.pdf", p3, width = 10, height = 10)

#### fragments path updates
pbmc_peaks@assays$peaks@fragments[[1]]@path <- "/NASdata/longep/cellranger_outs/M1/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[2]]@path <- "/NASdata/longep/cellranger_outs/M3/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[3]]@path <- "/NASdata/longep/cellranger_outs/M4/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[4]]@path <- "/NASdata/longep/cellranger_outs/M5/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[5]]@path <- "/NASdata/longep/cellranger_outs/M6/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[6]]@path <- "/NASdata/longep/cellranger_outs/M7/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[7]]@path <- "/NASdata/longep/cellranger_outs/M8/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[8]]@path <- "/NASdata/longep/cellranger_outs/M9/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[9]]@path <- "/NASdata/longep/cellranger_outs/NCI_13/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[10]]@path <- "/NASdata/longep/cellranger_outs/NCI_28/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[11]]@path <- "/NASdata/longep/cellranger_outs/NCI_55/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[12]]@path <- "/NASdata/longep/cellranger_outs/NCI_66/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[13]]@path <- "/NASdata/longep/cellranger_outs/NCI_118/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[14]]@path <- "/NASdata/longep/cellranger_outs/S1/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[15]]@path <- "/NASdata/longep/cellranger_outs/S2/outs/atac_fragments.tsv.gz"
pbmc_peaks@assays$peaks@fragments[[16]]@path <- "/NASdata/longep/cellranger_outs/N17_20/outs/atac_fragments.tsv.gz"


