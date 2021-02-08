#!/usr/bin/env Rscript

#### Tuft cell analysis: Reintegration; UMAPs, dot plots and diffusion maps; DGE and hypergeometric test
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
library(hypeR)

# Color palette
colors <- c('#006E82', '#AA0A3C', '#F0E442', '#00A0FA', '#FA5078', '#005AC8', '#CC79A7', '#FAE6BE', '#0072B2', '#A0FA82', '#F0F032', '#0AB45A', '#FA7850', '#14D2DC', '#FA78FA')

# Read in file
patient <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to basal, club and goblet cells and reintegrate
patient <- subset(patient, cell_type_fine %in% c('Airway basal', 'Airway club', 'Airway goblet'))
DefaultAssay(patient) <- 'RNA'
obj.list <- SplitObject(patient, split.by = 'patient')

# keep only patients with >50 cells
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 50) {
    obj.list[i] <- NA
  }
}
# remove NAs
obj.list <- obj.list[!is.na(obj.list)]

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, k.filter = 50)

# Integrate data sets
patient <- IntegrateData(anchorset = anchors, dims = 1:30)

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- RunUMAP(object = patient, dims = 1:30)

# Run gene signatures
sigs <- read.csv('data/epithelial_signatures.csv', na.strings = c('', 'NA'))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  patient <- AddModuleScore(object = patient, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Diffusion component analysis
es <- as.ExpressionSet(as.data.frame(t(patient@assays$integrated@data)))
es@phenoData@data <- patient@meta.data

# Make diffusion map
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

# Add DCs to Seurat object
patient <- AddMetaData(patient, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

# Identify tuft-like cells based on DC1 cut-off
patient$cell_type_fine <- ifelse(patient$DC1 > 0.015, 'Tuft-like', patient$cell_type_fine)

# Save UMAPS and diffusion maps
ifelse(!dir.exists(file.path('data/lungs_all/diffusion')), dir.create(file.path('data/lungs_all/diffusion')), FALSE)
ifelse(!dir.exists(file.path('data/lungs_all/diffusion/tuft_cells')), dir.create(file.path('data/lungs_all/diffusion/tuft_cells')), FALSE)
pdf('data/lungs_all/diffusion/tuft_cells/dm_tuft_int_cov_ctr.pdf')

# UMAPs
DimPlot(patient, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, cols = c('#E69F00', '#F0E442', '#56B4E9', '#FA5078'))
DimPlot(patient, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, cols = c('#E69F00', '#F0E442', '#56B4E9', '#FA5078')) + NoLegend()
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors) + ggtitle('Disease status')
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors) + ggtitle('Disease status') + NoLegend()

# Feature and dot plots with gene signatures
FeaturePlot(patient, features = c('club_cell1', 'basal_cell1', 'goblet_cell1', 'Tuft_montoro_long1'), 
            min.cutoff = 'q05', max.cutoff = 'q95', order = T, raster = T)
DotPlot(patient, features = c('club_cell1', 'basal_cell1', 'goblet_cell1', 'Tuft_montoro_long1'), 
        group.by = 'cell_type_fine', dot.scale = 7) + scale_color_viridis() + RotatedAxis() + coord_flip()

# Diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(c('#E69F00', '#F0E442', '#56B4E9', '#FA5078'))
plot(dm, col = as.factor(es@phenoData@data$cell_type_fine), main = 'Cell type (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$cell_type_fine)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$cell_type_fine))), bty = 'n')

palette(colors)
plot(dm, col = as.factor(es@phenoData@data$group), main = 'Disease status (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$group)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$group))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 6.5), xpd = TRUE)
palette(viridis(100))
plot(dm, col = es@phenoData@data$Tuft_montoro_long1, main = 'Tuft_montoro_long1 signature (COVID-19 and control cells)', pch = 20)
colorlegend(es@phenoData@data$Tuft_montoro_long1, viridis(length(es@phenoData@data$Tuft_montoro_long1)), posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = es@phenoData@data$club_cell1, main = 'club_cell1 signature (COVID-19 and control cells)', pch = 20)
colorlegend(es@phenoData@data$club_cell1, viridis(length(es@phenoData@data$club_cell1)), posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = es@phenoData@data$basal_cell1, main = 'basal_cell1 signature (COVID-19 and control cells)', pch = 20)
colorlegend(es@phenoData@data$basal_cell1, viridis(length(es@phenoData@data$basal_cell1)), posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = es@phenoData@data$goblet_cell1, main = 'goblet_cell1 signature (COVID-19 and control cells)', pch = 20)
colorlegend(es@phenoData@data$goblet_cell1, viridis(length(es@phenoData@data$goblet_cell1)), posx = c(0.88, 0.9), posy = c(0, 0.65))
dev.off()

# Save tuft cell image and label
save.image('data/lungs_all/diffusion/tuft_cells/image_tuft_int_cells.RData')
saveRDS(patient, 'data/lungs_all/diffusion/tuft_cells/data_tuft_int_cov_ctr.rds')
celltype <- patient@meta.data %>% select('barcode', 'cell_type_fine')
write.csv(celltype, 'data/lungs_all/data_lungs_all_tuft_celltype.csv')

# Preparation for differential gene expression
DefaultAssay(patient) <- 'RNA'
patient <- NormalizeData(patient)
patient <- ScaleData(patient)
DefaultAssay(patient) <- 'integrated'

# Identify markers
Idents(patient) <- patient$cell_type_fine
markers <- FindAllMarkers(patient, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
markers %>% group_by(cluster) %>% write.csv('data/lungs_all/diffusion/tuft_cells/markers_tuft.csv', row.names = F)
markers <- markers %>% filter(p_val_adj < 0.05)

# hypeR: Calculate enrichment of gene signatures
sigs2 <- read.csv('data/epithelial_signatures2.csv', na.strings = c('', NA))
genesets <- as.list(sigs2)
background <- rownames(patient@assays$RNA@data)

pdf('data/lungs_all/diffusion/tuft_cells/tuft_pathways.pdf')
tuft <- markers %>% dplyr::filter(cluster == 'Tuft-like') %>% magrittr::use_series(gene)
hyp_dots(hypeR(as.character(tuft), genesets = genesets, background = background), title = 'Epithelial signatures on tuft-like cells', abrv = 80)
dev.off()