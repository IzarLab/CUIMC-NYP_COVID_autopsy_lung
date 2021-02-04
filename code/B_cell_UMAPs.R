#!/usr/bin/env Rscript

#### B cell analysis: UMAPs
library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)

# Color palette
colors <- c('#006E82', '#AA0A3C', '#8214A0', '#00A0FA', '#FA5078', '#005AC8', '#CC79A7', '#FAE6BE', '#0072B2', '#A0FA82', '#F0F032', '#0AB45A', '#FA7850', '#14D2DC', '#FA78FA')

# Read in file
patient <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to B cells
patient <- subset(patient, cell_type_fine %in% c('B cells', 'Activated B cells', 'Plasma cells'))
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- RunUMAP(object = patient, dims = 1:30)

# Save UMAPs
ifelse(!dir.exists(file.path(paste0('data/lungs_all/b_cells'))), dir.create(file.path(paste0('data/lungs_all/b_cells'))), FALSE)
pdf('data/lungs_all/b_cells/umap_bcell_cov_ctr.pdf')

# UMAPs
DimPlot(patient, reduction = 'umap', label = F, group.by = 'cell_type_fine', raster = T, cols = c('#AF1900', '#E19600', '#193264'))
DimPlot(patient, reduction = 'umap', label = F, group.by = 'cell_type_fine', raster = T, cols = c('#AF1900', '#E19600', '#193264')) + NoLegend()
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', raster = T, cols = colors) + ggtitle('Disease status')
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', raster = T, cols = colors) + ggtitle('Disease status') + NoLegend()

# Dotplot
DotPlot(patient, assay = 'RNA', features = c('MS4A1', 'IL2RA', 'SDC1', 'IRF4'), group.by = 'cell_type_fine', dot.scale = 10) + scale_color_viridis() + RotatedAxis() + coord_flip()

dev.off()

# Save object
saveRDS(patient, 'data/lungs_all/b_cells/data_bcells_cov_ctr.rds')