#!/usr/bin/env Rscript

#### Fibroblast analysis: Diffusion maps
library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(gplots)
library(viridis)
library(scales)

# Color palette
colors <- c('#006E82', '#AA0A3C', '#F0E442', '#00A0FA', '#FA5078', '#005AC8', '#CC79A7', '#FAE6BE', '#0072B2', '#A0FA82', '#F0F032', '#0AB45A', '#FA7850', '#14D2DC', '#FA78FA')

# Read in file
patient <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to fibroblasts
patient <- subset(patient, cell_type_fine %in% c('Adventitial FB', 'Alveolar FB', 'Intermediate pathological FB', 'Other FB', 'Pathological FB'))

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- RunUMAP(object = patient, dims = 1:30)

# Diffusion component analysis
es <- as.ExpressionSet(as.data.frame(t(patient@assays$integrated@data)))
es@phenoData@data <- patient@meta.data

# Make diffusion map
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

# Save diffusion maps
ifelse(!dir.exists(file.path('data/lungs_all/diffusion')), dir.create(file.path('data/lungs_all/diffusion')), FALSE)
ifelse(!dir.exists(file.path('data/lungs_all/diffusion/fibroblasts')), dir.create(file.path('data/lungs_all/diffusion/fibroblasts')), FALSE)
pdf('data/lungs_all/diffusion/fibroblasts/!dm_fibroblasts_cov_ctr.pdf')

# Diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(c('#8214A0', '#00A0FA', '#0AB45A', '#14D2DC', '#FA7850'))
plot(dm, col = as.factor(es@phenoData@data$cell_type_fine), main = 'Cell type (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(-0.05, -0.12), legend = levels(as.factor(es@phenoData@data$cell_type_fine)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$cell_type_fine))), bty = 'n', cex = 0.85)

palette(colors)
plot(dm, col = as.factor(es@phenoData@data$group), main = 'Disease status (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(0, 0), legend = levels(as.factor(es@phenoData@data$group)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$group))), bty = 'n')
dev.off()

# Save image
save.image('data/lungs_all/diffusion/fibroblasts/image_fibroblasts.RData')
saveRDS(patient, 'data/lungs_all/diffusion/fibroblasts/data_fibroblasts_cov_ctr.rds')