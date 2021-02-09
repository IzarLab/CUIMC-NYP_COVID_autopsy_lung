#!/usr/bin/env Rscript

#### Global analysis: Sex differences
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(viridis)

# Read in file
seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Save UMAPS and dot plots
pdf(file = 'data/lungs_all/plots_global_sexdifferences.pdf')

# UMAP
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sex', shuffle = T, raster = T, cols = c('#009E73', '#D55E00'), pt.size = 0.001)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sex', shuffle = T, raster = T, cols = c('#009E73', '#D55E00'), pt.size = 0.001) + NoLegend()

# Dot plots
DotPlot(seu, assay = 'RNA', features = c('SARS-COV-2', 'ACE2', 'BSG', 'NPR1', 'TMPRSS2', 'FURIN', 'CTSL'), group.by = 'cell_type_intermediate', 
        dot.scale = 8, split.by = 'sex') + RotatedAxis() + ggtitle('All patients')

DotPlot(subset(seu, group == 'COVID-19' & sex == 'F'), assay = 'RNA', features = c('SARS-COV-2', 'ACE2', 'BSG', 'NPR1', 'TMPRSS2', 'FURIN', 'CTSL'), 
        group.by = 'cell_type_intermediate', dot.scale = 8) + scale_color_viridis() + RotatedAxis() + ggtitle('Female COVID-19 patients')
DotPlot(subset(seu, group == 'COVID-19' & sex == 'M'), assay = 'RNA', features = c('SARS-COV-2', 'ACE2', 'BSG', 'NPR1', 'TMPRSS2', 'FURIN', 'CTSL'), 
        group.by = 'cell_type_intermediate', dot.scale = 8) + scale_color_viridis() + RotatedAxis() + ggtitle('Male COVID-19 patients')

DotPlot(subset(seu, group == 'Control' & sex == 'F'), assay = 'RNA', features = c('SARS-COV-2', 'ACE2', 'BSG', 'NPR1', 'TMPRSS2', 'FURIN', 'CTSL'), 
        group.by = 'cell_type_intermediate', dot.scale = 8) + scale_color_viridis() + RotatedAxis() + ggtitle('Female control patients')
DotPlot(subset(seu, group == 'Control' & sex == 'M'), assay = 'RNA', features = c('SARS-COV-2', 'ACE2', 'BSG', 'NPR1', 'TMPRSS2', 'FURIN', 'CTSL'), 
        group.by = 'cell_type_intermediate', dot.scale = 8) + scale_color_viridis() + RotatedAxis() + ggtitle('Male control patients')

dev.off()