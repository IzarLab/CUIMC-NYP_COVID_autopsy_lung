#!/usr/bin/env Rscript

#### Global analysis: Harmony UMAPs
#### Author: Jana Biermann, PhD

library(Seurat)
library(dplyr)
library(ggplot2)

# Color palette
colors <- c('#006E82', '#AA0A3C', '#8214A0', '#00A0FA', '#FA5078', '#005AC8', '#CC79A7', '#FAE6BE', '#0072B2', '#A0FA82', '#F0F032', '#0AB45A', '#FA7850', '#14D2DC', '#FA78FA')

# Read in file and harmony coordinates (generated in Tcell_analysis.py)
seu <- readRDS('data/lungs_all/data_lungs_all_v4.rds')
harmony <- as.matrix(read.csv('data/lungs_all/Harmony_UMAP.csv', header = F))
rownames(harmony) <- colnames(seu@assays$RNA@data)

# Add harmony coordinates
seu@reductions$harmony <- CreateDimReducObject(embeddings = harmony, key = 'harmony_')

pdf(file = 'data/lungs_all/plots_harmony_umaps.pdf')
DimPlot(seu, reduction = 'harmony', group.by = 'initial_clustering', shuffle = T, raster = T, pt.size = 0.001)
DimPlot(seu, reduction = 'harmony', group.by = 'initial_clustering', shuffle = T, raster = T, pt.size = 0.001) + NoLegend()
DimPlot(seu, reduction = 'harmony', group.by = 'initial_clustering', label = T, repel = T, shuffle = T, raster = T, pt.size = 0.001) + NoLegend()

DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_main', shuffle = T, raster = T, pt.size = 0.001)
DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_main', shuffle = T, raster = T, pt.size = 0.001) + NoLegend()
DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_main', label = T, repel = T, shuffle = T, raster = T, pt.size = 0.001) + NoLegend()

DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_intermediate', shuffle = T, raster = T, pt.size = 0.001)
DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_intermediate', shuffle = T, raster = T, pt.size = 0.001) + NoLegend()
DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_intermediate', label = T, repel = T, shuffle = T, raster = T, pt.size = 0.001) + NoLegend()

DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_fine', shuffle = T, raster = T, pt.size = 0.001)
DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_fine', shuffle = T, raster = T, pt.size = 0.001) + NoLegend()
DimPlot(seu, reduction = 'harmony', group.by = 'cell_type_fine', label = T, repel = T, shuffle = T, raster = T, pt.size = 0.001) + NoLegend()

DimPlot(seu, reduction = 'harmony', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors, pt.size = 0.001)
DimPlot(seu, reduction = 'harmony', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors, pt.size = 0.001) + NoLegend()
dev.off()