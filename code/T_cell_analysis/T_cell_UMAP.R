#!/usr/bin/env Rscript

#### T cell analysis: UMAPs
library(ggplot2)
library(Seurat)
library(dplyr)

# Color palettes
colorsgroup <- c('#006E82', '#AA0A3C')
colors <- c('#96B400', '#FA5078', '#14D2DC', '#8214A0', '#A0FA82')

# Load seurat object
seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to T/NK cells
seu <- subset(seu, cell_type_main %in% c('T cells'))

# Rerun Seurat workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- RunUMAP(object = seu, dims = 1:30)

# Save UMAPS and diffusion maps
ifelse(!dir.exists(file.path('data/lungs_all/t_cells')), dir.create(file.path('data/lungs_all/t_cells')), FALSE)
pdf('data/lungs_all/t_cells/t_cells_umaps.pdf')

# UMAPs
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', repel = T, raster = T, cols = colors)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'group', repel = T, raster = T, cols = colorsgroup)

# Feature plots
FeaturePlot(seu, c('CD4', 'CD8A', 'FOXP3', 'NCAM1'), order = T, min.cutoff = 'q05', max.cutoff = 'q95', raster = T)
FeaturePlot(seu, c('GZMA', 'GZMB', 'TOP2A', 'MKI67'), order = T, min.cutoff = 'q05', max.cutoff = 'q95', raster = T)

# Dot plot
DotPlot(seu, assay = 'RNA', features = c('TRAC', 'CD4', 'CD8A', 'FOXP3', 'NCAM1'), group.by = 'cell_type_fine', 
        dot.scale = 8) + scale_color_viridis() + RotatedAxis() + coord_flip()

dev.off()