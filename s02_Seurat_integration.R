library(Seurat)
require(Matrix)
library(dplyr) 
require(stringr)
require(Signac)
library(rtracklayer)
library(GenomicRanges)

library(EnsDb.Hsapiens.v86)

# RNA
data_dir = "Chen/RNA/"
features <- read.delim(file.path(data_dir, "features.tsv.gz"), header = FALSE)$V1
barcodes <- read.delim(file.path(data_dir, "barcodes.tsv.gz"), header = FALSE)$V1
matrix <- readMM(file.path(data_dir, "matrix.mtx.gz"))
rownames(matrix) <- features
colnames(matrix) <- barcodes
rm(barcodes, features)

# Compute single RNA
rna <- CreateSeuratObject(counts = matrix, project = "GSE126074", min.cells = 1, min.features = 1)
rm(matrix)
DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 4000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20, reduction = "pca")
rna <- FindClusters(rna, resolution = 0.5)
rna <- RunUMAP(rna, reduction = "pca", dims = 1:15)
p1 <- DimPlot(
  rna,
  reduction = "umap",     
  group.by  = "seurat_clusters",  
  label     = TRUE,      
  pt.size   = 0.5         
) + NoLegend() 
p1


# ATAC
data_dir = "Chen/ATAC"
peaks <- read.delim(file.path(data_dir, "peaks.tsv.gz"), header = FALSE)$V1
barcodes <- read.delim(file.path(data_dir, "barcodes.tsv.gz"), header = FALSE)$V1
matrix <- readMM(file.path(data_dir, "matrix.mtx.gz"))
rownames(matrix) <- peaks
colnames(matrix) <- barcodes
rm(barcodes, peaks)


chrom_assay <- CreateChromatinAssay(
  counts = Matrix(matrix, sparse = TRUE),
  sep = c(":", "-"),
  fragments = NULL
)

# Compute single ATAC
rna[["ATAC"]] <- chrom_assay
DefaultAssay(rna) <- "ATAC"

rna <- RunTFIDF(rna)
rna <- FindTopFeatures(rna, min.cutoff = 'q0')
rna <- RunSVD(rna,n = 50)
rna <- FindNeighbors(rna, dims = 1:20, reduction = "lsi")
rna <- FindClusters(rna, resolution = 0.5)
rna <- RunUMAP(rna, reduction = "lsi", dims = 1:15)
p2 <- DimPlot(
  rna,
  reduction = "umap",     
  group.by  = "seurat_clusters",  
  label     = TRUE,      
  pt.size   = 0.5         
) + NoLegend() 
p2

# Integration
rna <- FindMultiModalNeighbors(rna, reduction.list = list("pca", "lsi"), dims.list = list(1:15, 1:15))

rna <- RunUMAP(rna, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
rna <- FindClusters(rna, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
p3 <- DimPlot(
  rna,
  reduction = "wnn.umap",     
  group.by  = "seurat_clusters",  
  label     = TRUE,      
  pt.size   = 0.5         
) + NoLegend() 
p3

writeMM(rna@graphs$wsnn,"Seurat_connectivities.mtx")
writeMM(rna@graphs$wknn,"Seurat_distance.mtx")

