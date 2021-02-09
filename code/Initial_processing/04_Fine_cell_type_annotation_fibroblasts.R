#!/usr/bin/env Rscript



#### script for Fine cell type annotation of fibroblast subset

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

seu<-readRDS('data/lungs_all/data_lungs_all.rds')

DimPlot(seu, reduction = "umap", label = T)
# Subset to fibroblasts
seu <- subset(seu, cell_type_main == "Fibroblasts")
DimPlot(seu, reduction = "umap", label = T)
seu@meta.data$orig.Seurat_clusters <- seu@meta.data$seurat_clusters
seu <- RunPCA(object = seu)
ElbowPlot(seu, ndims = 50)
ndim = 15
seu <- FindNeighbors(seu, dims = 1:ndim)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(object = seu, dims = 1:ndim)

DimPlot(seu, reduction = "umap", label = T)
DimPlot(seu, reduction = "umap", label = T, group.by = "orig.Seurat_clusters")
DimPlot(seu, reduction = "umap", label = T, group.by = "cell_type_fine")

# getting previous fibroblast cluster assignments
lung_fibroblast_annotation <- read.csv("annotation_required/lung_fibroblast_annotation.csv", 
                                       row.names = 1)
head(lung_fibroblast_annotation)
all(rownames(lung_fibroblast_annotation) == rownames(seu@meta.data))
seu@meta.data$seurat_clusters <- as.factor(lung_fibroblast_annotation$seurat_clusters)
DimPlot(seu, reduction = "umap", label = T, group.by = "seurat_clusters")
seu <- SetIdent(seu, value = "seurat_clusters")
DimPlot(seu, reduction = "umap", label = T)

# Cell cycle assignment
DefaultAssay(seu) <- "RNA"
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, 
                        g2m.features = cc.genes.updated.2019$g2m.genes, 
                        set.ident = F)
DefaultAssay(seu) <- "integrated"
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 
                           0.1, "cycling", "non-cycling")

# Differential gene expression (DGE) based on
# clusters
seu@active.assay <- "RNA"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25)
write.csv(markers, "data/lungs_all/markers_lungs_all_fibro.csv", 
          row.names = F)

# Differential gene expression (DGE) based on
# pathogenic and non-pathogenic clusters
seu <- SetIdent(seu, value = "cell_type_fine")
markers <- FindMarkers(seu, ident.1 = c("Pathological FB", 
                                        "Intermediate pathological FB"), only.pos = FALSE, 
                       min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "data/lungs_all/markers_lungs_pathological_fibro.csv", 
          row.names = T)
seu <- SetIdent(seu, value = "seurat_clusters")

# Fine cell type assignment after manual annotation
# based on DGE
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(0, 1), "Intermediate pathological FB", "NA")
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(6), "Other FB", seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(4), "Pathological FB", seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(2), "Adventitial FB", seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(3, 5), "Alveolar FB", seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(7), "Airway smooth muscle", seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(8), "Pericyte", seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(9), "Vascular smooth muscle", seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$seurat_clusters %in% 
                               c(10), "Mesothelial FB", seu$cell_type_fine)

DimPlot(seu, group.by = "cell_type_fine")

# Save object and labels
saveRDS(seu, "data/lungs_all/data_lungs_all_fibro.rds")
celltype <- seu@meta.data %>% select("barcode", "cell_cycle", 
                                     "cell_type_fine")
write.csv(celltype, "data/lungs_all/data_lungs_all_fibro_celltype.csv")