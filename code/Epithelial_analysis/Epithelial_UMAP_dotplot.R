#!/usr/bin/env Rscript

#### Epithelial cell analysis: UMAPs and dot plots
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(viridis)

# Color palettes
colors <- c('#006E82', '#AA0A3C', '#8214A0', '#00A0FA', '#FA5078', '#005AC8', '#CC79A7', '#FAE6BE', '#0072B2', '#A0FA82', 
            '#F0F032', '#0AB45A', '#FA7850', '#14D2DC', '#FA78FA')
epi_palette <- c('#E69F00', '#CC79A7', '#F0E442', '#56B4E9', '#D55E00', '#0072B2', '#0AB45A', '#A0FA82', '#8214A0')

# Read in file
seu <- readRDS('data/lungs_all/data_lungs_all_epithelial.rds')

# Rename labels
seu$group <- ifelse(seu$group == 'cov', 'COVID-19', 'Control')

# Save UMAPS and dot plots
pdf(file = 'data/lungs_all/plots_epithelial_analysis.pdf')

# UMAPs
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, cols = epi_palette)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, cols = epi_palette) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors) + NoLegend()

# Dot plots
DotPlot(seu, assay = 'RNA', features = c('basal_cell1', 'ciliated_cell1', 'club_cell1', 'goblet_cell1', 'mucous_cell1', 'AT11', 'AT21'), 
        group.by = 'cell_type_fine', dot.scale = 8) + scale_color_viridis() + RotatedAxis() + coord_flip()

DotPlot(seu, assay = 'RNA', features = c('MKI67', 'TOP2A', 'B2M', 'HLA-DRA', 'IRF8'), group.by = 'cell_type_fine', dot.scale = 8) + 
  scale_color_viridis() + RotatedAxis() + coord_flip()

DotPlot(seu, assay = 'RNA', features = c('COL1A1', 'COL1A2', 'COL14A1', 'COL6A3'), group.by = 'cell_type_fine', dot.scale = 8) + 
  scale_color_viridis() + RotatedAxis() + coord_flip()

dev.off()
