# Load libraries
knitr::opts_chunk$set(echo = TRUE)
library(RColorBrewer)
library(Seurat)
library(DoubletFinder)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(SingleR)
library(celldex)
library(speckle)
library(xlsx)
library(scales)
library(hdf5r)
ref.se <- celldex::ImmGenData()

# Preprocess each sample
setwd("/Volumes/LaCie/qian_dock8/")
sample_list <- c("DOCK8_H962R", "DOCK8_N1987Y", "DOCK8_WT", "DOCK8_stat3flox")
so_list <- list()
for (isample in sample_list) {
  set.seed(123) 
  # Load raw 10X h5 matrix (change the path according to the sample_filtered_feature_bc_matrix.h5 location for each sample
  myfile <- paste0("/path/", isample,
                   "/Per_sample_outputs/", isample,
                   "_001/Count outputs/sample_filtered_feature_bc_matrix.h5")
  h5_obj <- Read10X_h5(filename = myfile, use.names = TRUE, unique.features = TRUE)
  # Create Seurat object
  Seurat_obj <- CreateSeuratObject(counts = h5_obj, project = isample, 
                                   min.cells = 0, min.features = 200)
  Seurat_obj <- RenameCells(Seurat_obj, add.cell.id = isample)
  # Calculate QC metrics
  Seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-")
  VlnPlot(Seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  # Filter cells
  Seurat_obj <- subset(Seurat_obj,
    subset = nCount_RNA > 500 & nCount_RNA < 30000 &
             nFeature_RNA > 500 & nFeature_RNA < 5000 &
             percent.mt < 10)
  # Normalize and PCA
  Seurat_obj <- NormalizeData(Seurat_obj) %>%
                FindVariableFeatures(nfeatures = 3000) %>%
                ScaleData() %>%
                RunPCA(verbose = FALSE)
  ElbowPlot(Seurat_obj, ndims = 50)
  # Clustering and UMAP
  ndims <- 40
  Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:ndims) %>%
                FindClusters(resolution = c(0.2, 0.5, 0.7)) %>%
                RunUMAP(dims = 1:ndims)
  DimPlot(Seurat_obj, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)
  # Run DoubletFinder
  sweep.res.list <- paramSweep(Seurat_obj, PCs = 1:ndims, sct = FALSE, num.cores = 8)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    pull(pK) %>%
    as.numeric()
  annotations <- Seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * ncol(Seurat_obj))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  Seurat_obj <- doubletFinder(Seurat_obj, PCs = 1:ndims, pN = 0.25,
                              pK = pK, nExp = nExp_poi.adj)
  # Keep only singlets
  Seurat_obj <- subset(Seurat_obj, subset = DF.classifications == "Singlet")
  so_list[[isample]] <- Seurat_obj
}
# Save individual objects
saveRDS(so_list, file = "/path/rds/so_list.rds")

# Merge datasets
options(future.globals.maxSize = 30 * 1024^3)
so_list <- readRDS("/path/rds/so_list.rds")
combined <- merge(x = so_list[[1]], y = so_list[2:4])
print(combined)

#set dimension
ndims = 40

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 3000)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
ElbowPlot(combined, ndims= 50)
combined <- FindNeighbors(combined, dims = ndims, reduction = "pca")
combined <- FindClusters(combined, resolution =0.2, cluster.name = "unintegrated_clusters")
combined <- RunUMAP(combined, dims = 1:ndims, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(combined, reduction = "umap.unintegrated", group.by = c("orig.ident"))

# Perform integration

options(future.globals.maxSize = 30 * 1024^3) 
obj <- IntegrateLayers(
    object = combined, method = HarmonyIntegration,   
  new.reduction = "harmony", verbose = FALSE, features = all_genes)

dims= 50
set.seed(1234) 
integrated <- FindNeighbors(integrated, reduction= "harmony", dims=1:dims)
ElbowPlot(integrated, ndims = 50)
integrated <- FindClusters(integrated, resolution= 0.5, cluster.name= "harmony_0.5")
integrated <- FindClusters(integrated, resolution= 0.6, cluster.name= "harmony_0.6")
integrated <- FindClusters(integrated, resolution= 0.7, cluster.name= "harmony_0.7")
integrated <- FindClusters(integrated, resolution= 0.65, cluster.name= "harmony_0.65")
integrated <- FindClusters(integrated, resolution= 0.8, cluster.name= "harmony_0.8")
integrated <- FindClusters(integrated, resolution= 0.9, cluster.name= "harmony_0.9") 
integrated <- RunUMAP(integrated, reduction= "harmony", reduction.name= "umap_harmony", dims = 1:dims)

integrated@meta.data$sample<-integrated$orig.ident

DimPlot(
  integrated,
  reduction = "umap_harmony",
  group.by = "harmony_0.6",
  combine = FALSE, label.size = 2)

p2 <- DimPlot(
  integrated,
  reduction = "umap_harmony",
  group.by = c("sample"),
  combine = FALSE, label.size = 2)


integrated <- JoinLayers(integrated)


## Save the integrated Seurat object with processing and clustering
saveRDS(integrated, file = "/path/rds/seurat.integrated_processed_harmony_50dims.rds")
integrated = readRDS("/path/rds/seurat.integrated_processed_harmony_50dims.rds")

# Find allMarkers
library(dplyr)
Idents(obj) <- integrated$harmony_0.6
all.markers_0.6 <-FindAllMarkers(integrated,logfc.threshold = 0.1,min.pct = 0.25, thresh.use = 0.25,
                               only.pos = TRUE)
features_06<-all.markers_0.6 %>%
  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>% 
  slice_max(n = 5, order_by = avg_log2FC)

p6<-DotPlot(obj, features = unique(features_06$gene)) + RotatedAxis()
p6

#check top 25 markers
features25_06<-all.markers_0.6 %>%
  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>% 
  slice_max(n = 25, order_by = avg_log2FC)

#check top 50 markers
features50_06<-all.markers_0.6 %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)


