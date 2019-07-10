## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning=FALSE, 
  message = FALSE,
  fig.width = 7,
  fig.align = "center"
  )

## ----setup---------------------------------------------------------------
library(clustably)
library(Seurat)

## ----load_pbmc-----------------------------------------------------------

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc

## ----default_clustering--------------------------------------------------
set.seed(42)
n <- length(Cells(pbmc))
sample.size <- floor(.75*n)
pbmc.subset = subset(pbmc, cells = sample(Cells(pbmc), sample.size))

doClustering <- function(obj) {
  obj <- NormalizeData(obj) 
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = F) 
  obj <- RunUMAP(obj, reduction="pca", dims=1:10, verbose = F) 
  obj <- FindNeighbors(obj, dims=1:10, verbose = F) 
  obj <- FindClusters(obj, verbose=F)
}

pbmc <- doClustering(pbmc)
pbmc.subset <- doClustering(pbmc.subset)

## ----compare_dims--------------------------------------------------------
CombinePlots(list(DimConfidencePlot(pbmc), DimConfidencePlot(pbmc.subset)))

## ----plot_naive_alignment, fig.width=4-----------------------------------
library(ggplot2)
tab <- table(Idents(pbmc)[Cells(pbmc.subset)],Idents(pbmc.subset), dnn=c("full","subset")) 
df <- as.data.frame(scale(tab, center=F))
ggplot(df, aes(full,subset,fill=Freq)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "OrRd", direction=1)

## ----alignment-----------------------------------------------------------
pbmc.subset.aligned <- align2Idents.Seurat(pbmc.subset, pbmc)

CombinePlots(list(
    PlotCrossClassification(pbmc.subset.aligned, pbmc,  names=c("subset", "orig")),
    PlotCrossClassification(pbmc.subset.aligned, pbmc,  names=c("subset", "orig"), scale = T)
  )
)
  

## ----jackknife-----------------------------------------------------------
label <- function(obj) {
  obj <- NormalizeData(obj, verbose=F) 
  obj <- FindVariableFeatures(obj, verbose = F) 
  obj <- ScaleData(obj, verbose = F) 
  obj <- RunPCA(obj, verbose = F) 
  obj <- FindNeighbors(obj, verbose = F)
  obj <- FindClusters(obj, verbose = F) 
  Idents(obj)
}

pbmc.jack <- jackknifeClustering(pbmc, labelingFunc=label, n=100, verbose=T)
pbmc.jack <- RunUMAP(pbmc.jack, reduction = "pca", dims = 1:10)

## ---- fig.width=7, fig.height=7------------------------------------------
Idents(pbmc.jack) <- "jackknifeConsensus"
DimConfidencePlot(pbmc.jack, confidence="jackknifeConfidence")

## ----crosslabeling, include=F--------------------------------------------
# This figure is too big.  We can also plot the frequency with which the labels are crossed. 
DimCrossClassificationPlot(pbmc.jack)

## ------------------------------------------------------------------------
p1 <- PlotCrossClassificationFrequency(pbmc.jack)
p2 <- PlotCellClassifications(pbmc.jack)
CombinePlots(list(p1,p2))

## ------------------------------------------------------------------------
p1 <- PlotConfidenceDistributions(pbmc.jack)
p2 <- PlotConfidenceDistributions(pbmc.jack, violin=T)
CombinePlots(list(p1,p2))

## ------------------------------------------------------------------------
label.opt <- function(obj) {
  obj <- NormalizeData(obj, verbose=F) 
  obj <- FindVariableFeatures(obj, verbose = F) 
  obj <- ScaleData(obj, verbose = F) 
  obj <- RunPCA(obj, verbose = F) 
  obj <- FindNeighbors(obj, verbose = F)
  obj <- FindClusters(obj, resolution = 0.5, verbose = F) 
  Idents(obj)
}

pbmc.jack.opt <- jackknifeClustering(pbmc, labelingFunc=label.opt, n=100, verbose=T)
pbmc.jack.opt <- RunUMAP(pbmc.jack.opt, reduction = "pca", dims = 1:10)

## ------------------------------------------------------------------------
CombinePlots(list(
  PlotCrossClassificationFrequency(pbmc.jack, plotTitle="Default"),
  PlotCrossClassificationFrequency(pbmc.jack.opt, plotTitle="Optimized")))

## ------------------------------------------------------------------------
CombinePlots(list(
  PlotCellClassifications(pbmc.jack, plotTitle="Default"),
  PlotCellClassifications(pbmc.jack.opt, plotTitle="Optimized")))

## ------------------------------------------------------------------------
CombinePlots(list(
  PlotConfidenceDistributions(pbmc.jack, plotTitle="Default"),
  PlotConfidenceDistributions(pbmc.jack.opt, plotTitle="Optimized")))

