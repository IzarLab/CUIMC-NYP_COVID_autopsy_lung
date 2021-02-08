#!/usr/bin/env Rscript

#### Myeloid analysis: UMAPs and diffusion maps; Identification of AMs
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
colors <- c('#006E82', '#AA0A3C', '#8214A0', '#00A0FA', '#FA5078', '#005AC8', '#CC79A7', '#FAE6BE', '#0072B2', '#A0FA82', '#F0F032', '#0AB45A', '#FA7850', '#14D2DC', '#FA78FA')

# Read in file
patient <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to monocytes and macrophages
patient <- subset(patient, cell_type_fine %in% c('Monocytes', 'Macrophages', 'Transitioning MDM', 'Monocyte-derived macrophages'))
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- RunUMAP(object = patient, dims = 1:30)

# Identify macrophage lineage 
# Gene signature to identify alveolar macrophages (AM) based on differential gene expression was obtained from Travaglini et al. (Travaglini et al., 2020)
sigs <- read.csv('data/myeloid_signatures_cov.csv', na.strings = c(''))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  patient <- AddModuleScore(object = patient, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

patient$cell_type_fine <- ifelse(patient$alveolar_mac_DEG1 > 0.15, 'Alveolar macrophages', 
                                 ifelse(patient$cell_type_fine %in% c('Macrophages', 'Monocyte-derived macrophages'), 'Monocyte-derived macrophages',
                                        ifelse(patient$cell_type_fine == 'Monocytes', patient$cell_type_fine, 'Transitioning MDM')))


# Diffusion component analysis
es <- as.ExpressionSet(as.data.frame(t(patient@assays$integrated@data)))
es@phenoData@data <- patient@meta.data

# Make a diffusion map
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

# Save UMAPS and diffusion maps
ifelse(!dir.exists(file.path('data/lungs_all/diffusion')), dir.create(file.path('data/lungs_all/diffusion')), FALSE)
ifelse(!dir.exists(file.path('data/lungs_all/diffusion/myeloid')), dir.create(file.path('data/lungs_all/diffusion/myeloid')), FALSE)
pdf('data/lungs_all/diffusion/myeloid/dm_myeloid_cov_ctr.pdf')

# UMAPs
DimPlot(patient, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, cols = c('#E73F74', '#F2B701', '#11A579', '#3969AC'))
DimPlot(patient, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, cols = c('#E73F74', '#F2B701', '#11A579', '#3969AC')) + NoLegend()
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors) + ggtitle('Disease status')
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors) + ggtitle('Disease status') + NoLegend()

# Feature plot
FeaturePlot(patient, features = c('alveolar_mac_DEG1', 'Monocytes1', 'macrophage1'), min.cutoff = 'q05', max.cutoff = 'q95', order = T, raster = T)

# Diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(c('#E73F74', '#F2B701', '#11A579', '#3969AC'))
plot(dm, col = as.factor(es@phenoData@data$cell_type_fine), main = 'cell_type_fine (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(-0.05, -0.15), legend = levels(as.factor(es@phenoData@data$cell_type_fine)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$cell_type_fine))), bty = 'n', cex = 0.8)

par(mar = c(5.1, 4.1, 4.1, 6.5), xpd = TRUE)
palette(viridis(100))
plot(dm, col = es@phenoData@data$alveolar_mac_DEG1, main = 'alveolar_mac_DEG signature (COVID-19 and control cells)', pch = 20)
colorlegend(es@phenoData@data$alveolar_mac_DEG1, viridis(length(es@phenoData@data$alveolar_mac_DEG1)), posx = c(0.88, 0.9), posy = c(0, 0.65))

par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colors)
plot(dm, col = as.factor(es@phenoData@data$group), main = 'Disease status (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$group)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$group))), bty = 'n')
dev.off()

# Save myeloid cell image and label
save.image('data/lungs_all/diffusion/myeloid/image_myeloid.RData')
saveRDS(patient, 'data/lungs_all/diffusion/myeloid/data_myeloid_cov_ctr.rds')
celltype <- patient@meta.data %>% select('barcode', 'cell_type_fine')
write.csv(celltype, 'data/lungs_all/data_lungs_all_macs_celltype.csv')