#!/usr/bin/env Rscript

### title: Analysis of T cells in lung samples (7 COVID-19 and 1 control samples). Code for Figure 2ij.
### author: Jana Biermann
### date: 10/05/20

print(Sys.time())

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gplots)
library(viridis)


pat <- 'lungs_all'

#### UMAP analysis (Figure 2ij) ####
# Read in file from script #3
patient <- readRDS(paste0('data/', pat, '/data_', pat, '_with_covid_sigs.rds'))
patient@meta.data$group <- ifelse(patient@meta.data$group == 'cov', 'COVID-19', 'Control')

# Subset to T and NK cells as identified in overallclassification and celltype_bped_fine (SingleR)
patient <- subset(patient, overallclassification %in% c('CD8+ T-cells', 'CD4+ T-cells', 'NK cells'))
patient <- subset(patient, celltype_bped_fine %in% c('CD8+ T-cells', 'CD8+ Tcm', 'CD8+ Tem', 'CD4+ T-cells', 'CD4+ Tcm', 'CD4+ Tem', 'NK cells', 'Tregs'))

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- FindNeighbors(patient, dims = 1:20)
patient <- FindClusters(patient, resolution = 0.8)
patient <- RunUMAP(object = patient, dims = 1:20)

# Exclude outliers
patient <- subset(patient, integrated_snn_res.0.8 != c(12))

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- FindNeighbors(patient, dims = 1:20)
patient <- FindClusters(patient, resolution = 0.8)
patient <- RunUMAP(object = patient, dims = 1:20)

# Add Tregs to 'overallclassification'
patient@meta.data$overallclassification <- ifelse(patient@meta.data$celltype_bped_fine == 'Tregs', 'Tregs', patient@meta.data$overallclassification)

# Save Figure 2ij
pdf(file = paste0('data/', pat, '/t_cells/Figure_2ij.pdf'))
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', repel = T, cols = rev(viridis(3)[1:2]))
DimPlot(patient, reduction = 'umap', label = F, group.by = 'overallclassification', repel = T, cols = rev(viridis(4)))
dev.off()


##### Save T cell image #####
save.image(paste0('data/', pat, '/t_cells/image_t_cells.RData'))


print(Sys.time())
