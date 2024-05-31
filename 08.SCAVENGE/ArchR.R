##############################################################################
# Script information                                                      
# Title: Processing snATAC-seq data
# Author: Erping Long
# Date: 2023-5-14
# Description: None
##############################################################################

library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)

addArchRThreads(threads = 16)
addArchRGenome("hg38")

#creating Arrow Files

ArrowFiles <- createArrowFiles(
  inputFiles = "atac_fragments.tsv.gz",
  sampleNames = "MS1",
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
)

#inferring doublets

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

#creating an ArchRProject

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE
)

#removing doublets

proj <- filterDoublets(ArchRProj = proj)

#dimentional reduction and clustering

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

plot1 <- ggAlignPlots(p2, type = "h")

plotPDF(p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj <- saveArchRProject(ArchRProj = proj)

projHeme4 <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters") 

#calling peaks

projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "Clusters",
    pathToMacs2 = "/python3.9/site-packages/MACS3"
)


projHeme5 <- addPeakMatrix(projHeme4)

getAvailableMatrices(ArchRProj = projHeme5)

proj_PeakMatrix <- getMatrixFromProject(
  ArchRProj = projHeme5,
  useMatrix = "PeakMatrix",
)

dim(proj_PeakMatrix)
class(proj_PeakMatrix)

save(proj_PeakMatrix, file="peak_by_cell_matrix-summarizedexperiment.rda")

#obtaining peak-by-cell matrix of snATAC-seq

peakbycellmat <- assay(proj_PeakMatrix) 

#storing the matrix in a RangedSummarizedExperiment object

SE_gvar <- SummarizedExperiment(assays = list(counts = peakbycellmat),
                           rowRanges = rowRanges(proj_PeakMatrix), 
                           colData = DataFrame(names = colnames(peakbycellmat)))

assayNames(SE_gvar) <- "counts"

SE_gvar <- addGCBias(SE_gvar, genome = BSgenome.Hsapiens.UCSC.hg38)
SE_gvar_bg <- getBackgroundPeaks(SE_gvar, niterations=200)
save(SE_gvar, SE_gvar_bg, file="NCI_M1_SE_gvar.rda") 
