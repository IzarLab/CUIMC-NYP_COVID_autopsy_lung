#!/usr/bin/env Rscript

#### AT cell analysis: Reintegration; UMAPs, diffusion maps, dot, feature and violin plots; DGE and hypergeometric test; DATP cell state
#### fractions
#### Author: Jana Biermann, PhD

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
colors <- c('#006E82', '#AA0A3C', '#F0E442', '#00A0FA', '#FA5078', '#005AC8', '#CC79A7', '#FAE6BE', '#0072B2', '#A0FA82', '#F0F032', '#0AB45A', 
            '#FA7850', '#14D2DC', '#FA78FA')

# Read in file
patient <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to AT1 and AT2 cells
patient <- subset(patient, cell_type_fine %in% c('AT1', 'AT2'))
DefaultAssay(patient) <- 'RNA'
obj.list <- SplitObject(patient, split.by = 'patient')

# Keep only patient with > 50 cells (all patients have > 50 cells in data set)
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 50) {
    obj.list[i] <- NA
  }
}

# Remove NAs
obj.list <- obj.list[!is.na(obj.list)]

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, k.filter = 60)

# Integrate data sets
patient <- IntegrateData(anchorset = anchors, dims = 1:30)

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- RunUMAP(object = patient, dims = 1:30)

# Identify AT1/2 cell subtypes Signatures to identify DATP cells, primed and cycling AT2 cells were obtained from Choi et al. (Choi et al.,
# 2020)
sigs <- read.csv('data/epithelial_signatures.csv', na.strings = c('', 'NA'))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  patient <- AddModuleScore(object = patient, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# DATP initial classification based on 'DATP_short1' signature
patient$cell_type_datp <- ifelse(patient$DATP_short1 > 0.7, 'DATP', patient$cell_type_fine)

# DATP final classification based on our DATP signature
patient$cell_type_ourdatp <- ifelse(patient$our_DATP_sig1 > 0.4, 'DATP', patient$cell_type_fine)

# Subset to cell states and rerun CellCycleScoring
datp <- subset(patient, cell_type_ourdatp == 'DATP') %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
DefaultAssay(datp) <- 'RNA'
datp <- CellCycleScoring(datp, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(datp) <- 'integrated'

at <- subset(patient, cell_type_ourdatp %in% c('AT1', 'AT2')) %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
DefaultAssay(at) <- 'RNA'
at <- CellCycleScoring(at, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(at) <- 'integrated'

# Diffusion component analysis
es <- as.ExpressionSet(as.data.frame(t(patient@assays$integrated@data)))
es@phenoData@data <- patient@meta.data

# Make diffusion map
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

# Add DCs to Seurat object
patient <- AddMetaData(patient, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

# Save UMAPS and diffusion maps
ifelse(!dir.exists(file.path('data/lungs_all/diffusion')), dir.create(file.path('data/lungs_all/diffusion')), FALSE)
ifelse(!dir.exists(file.path('data/lungs_all/diffusion/AT')), dir.create(file.path('data/lungs_all/diffusion/AT')), FALSE)
pdf('data/lungs_all/diffusion/AT/dm_at_int_cov_ctr.pdf')

# UMAPs
DimPlot(patient, reduction = 'umap', label = F, group.by = 'cell_type_ourdatp', shuffle = T, raster = T, cols = c('#0072B2', '#0AB45A', '#FA5078'))
DimPlot(patient, reduction = 'umap', label = F, group.by = 'cell_type_ourdatp', shuffle = T, raster = T, cols = c('#0072B2', '#0AB45A', '#FA5078')) + 
  NoLegend()
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors) + ggtitle('Disease status')
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', shuffle = T, raster = T, cols = colors) + ggtitle('Disease status') + NoLegend()

# Feature plots
FeaturePlot(patient, features = c('UP_in_DATPs1', 'DATP_short1', 'DATP_long1', 'our_DATP_sig1'), min.cutoff = 'q05', max.cutoff = 'q95', order = T, 
            raster = T)
FeaturePlot(patient, features = c('MKI67', 'TOP2A', 'primed_AT21', 'cycling_AT21'), min.cutoff = 'q05', max.cutoff = 'q95', order = T, raster = T)
FeaturePlot(patient, features = c('CLDN4', 'KRT8', 'CDKN1A'), min.cutoff = 'q05', max.cutoff = 'q95', order = T, raster = T)
FeatureScatter(patient, feature1 = 'G2M.Score', feature2 = 'S.Score', group.by = 'cell_type_ourdatp', cols = c('#0072B2', '#0AB45A', '#FA5078'), 
               pt.size = 0.2)
FeatureScatter(patient, feature1 = 'G2M.Score', feature2 = 'S.Score', group.by = 'cell_type_fine', cols = c('#0072B2', '#0AB45A'), pt.size = 0.2)
FeatureScatter(datp, feature1 = 'G2M.Score', feature2 = 'S.Score', group.by = 'cell_type_ourdatp', cols = c('#FA5078'), pt.size = 0.2)
FeaturePlot(datp, features = c('MKI67', 'TOP2A', 'G2M.Score', 'S.Score'), min.cutoff = 'q05', max.cutoff = 'q95', order = T, raster = T)
FeaturePlot(at, features = c('MKI67', 'TOP2A', 'G2M.Score', 'S.Score'), min.cutoff = 'q05', max.cutoff = 'q95', order = T, raster = T)

# Violin and dot plots
VlnPlot(object = subset(patient, cell_type_fine == 'AT1'), features = c('CAV1', 'ETV5', 'SFTPB', 'our_DATP_sig1'), ncol = 4, pt.size = 0, group.by = 'cell_type_fine', 
        cols = colors, split.by = 'group')
VlnPlot(object = subset(patient, cell_type_fine == 'AT2'), features = c('CAV1', 'ETV5', 'SFTPB', 'our_DATP_sig1'), ncol = 4, pt.size = 0, group.by = 'cell_type_fine', 
        cols = colors, split.by = 'group')
VlnPlot(object = patient, features = c('CAV1', 'ETV5', 'SFTPB', 'our_DATP_sig1'), ncol = 4, pt.size = 0, cols = colors, group.by = 'group')
DotPlot(patient, features = c('AT11', 'AT21', 'DATP_short1', 'our_DATP_sig1'), group.by = 'cell_type_ourdatp', dot.scale = 7) + scale_color_viridis() + 
  RotatedAxis() + coord_flip()

# Diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(c('#0072B2', '#0AB45A'))
plot(dm, col = as.factor(es@phenoData@data$cell_type_fine), main = 'Cell type (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(0, 0), legend = levels(as.factor(es@phenoData@data$cell_type_fine)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$cell_type_fine))), 
       bty = 'n')

palette(colors)
plot(dm, col = as.factor(es@phenoData@data$group), main = 'Disease status (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$group)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$group))), 
       bty = 'n')

palette(c('#0072B2', '#0AB45A', '#FA5078'))
plot(dm, col = as.factor(es@phenoData@data$cell_type_ourdatp), main = 'cell_type_ourdatp (COVID-19 and control cells)', pch = 20)
legend('bottomright', inset = c(0, 0), legend = levels(as.factor(es@phenoData@data$cell_type_ourdatp)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$cell_type_ourdatp))), 
       bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
plot(dm, col = es@phenoData@data$DATP_short1, main = 'DATP_short1 COVID-19 and control cells', pch = 20)
colorlegend(es@phenoData@data$DATP_short1, viridis(length(es@phenoData@data$DATP_short1)), posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = es@phenoData@data$our_DATP_sig1, main = 'our_DATP_sig1 COVID-19 and control cells', pch = 20)
colorlegend(es@phenoData@data$our_DATP_sig1, viridis(length(es@phenoData@data$our_DATP_sig1)), posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = es@phenoData@data$AT11, main = 'AT1 signature COVID-19 and control cells', pch = 20)
colorlegend(es@phenoData@data$AT11, viridis(length(es@phenoData@data$AT11)), posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = es@phenoData@data$AT21, main = 'AT2 signature COVID-19 and control cells', pch = 20)
colorlegend(es@phenoData@data$AT21, viridis(length(es@phenoData@data$AT21)), posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = es@phenoData@data$primed_AT21, main = 'Primed AT2 in COVID-19 and control cells', pch = 20)
colorlegend(es@phenoData@data$primed_AT21, viridis(length(es@phenoData@data$primed_AT21)), posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = es@phenoData@data$cycling_AT21, main = 'Cycling AT2 COVID-19 and control cells', pch = 20)
colorlegend(es@phenoData@data$cycling_AT21, viridis(length(es@phenoData@data$cycling_AT21)), posx = c(0.88, 0.9), posy = c(0, 0.65))
dev.off()

# Save AT cell image
save.image('data/lungs_all/diffusion/AT/image_AT_cells_int.RData')
saveRDS(patient, 'data/lungs_all/diffusion/AT/data_AT_int_cov_ctr.rds')


### DATP signature generation Add scaled data for DEG
DefaultAssay(patient) <- 'RNA'
patient <- NormalizeData(patient)
patient <- ScaleData(patient)
DefaultAssay(patient) <- 'integrated'

# Identify markers in initial DATP signature
Idents(patient) <- patient$cell_type_datp
markers <- FindAllMarkers(patient, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
markers %>% group_by(cluster) %>% write.csv('data/lungs_all/diffusion/AT/markers_datp.csv', row.names = F)
markers <- markers %>% filter(p_val_adj < 0.05)  # add significant DATP markers to 'data/epithelial_signatures.csv' and run 'AddModuleScore'

# Pathway analysis for markers in our_DATP_sig
HALLMARK <- msigdb_gsets(species = 'Homo sapiens', category = 'H')

pdf('data/lungs_all/diffusion/AT/our_DATP_sig_pathways.pdf')
background <- rownames(patient@assays$RNA@data)
datp <- markers %>% dplyr::filter(cluster == 'DATP') %>% magrittr::use_series(gene)
hyp_dots(hypeR(as.character(datp), genesets = HALLMARK, background = background), title = 'Hallmarks in DATP cells ', abrv = 80)
dev.off()


### DATP quantification
df_tobesummed_fine = data.frame(orig.ident = patient$orig.ident, group = patient$group, cell_type_ourdatp = patient$cell_type_ourdatp)
df_summed_fine = df_tobesummed_fine %>% group_by(orig.ident, cell_type_ourdatp, group) %>% tally()
df_summed_fine = df_summed_fine %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

# Wilcoxon test
pdatp <- subset(df_summed_fine, cell_type_ourdatp == 'DATP')
wilcox.test(pdatp$freq ~ pdatp$group)$p.value

# Fractions of cell states
write.csv(prop.table(table(patient$cell_type_ourdatp, patient$group)), 'data/lungs_all/cell_type_fraction_datp_among_AT1AT2.csv', row.names = T)
write.csv(prop.table(table(subset(patient, group == 'COVID-19')$cell_type_ourdatp)), 'data/lungs_all/cell_type_fraction_datp_among_AT1AT2_cov.csv', 
          row.names = F)
write.csv(prop.table(table(subset(patient, group == 'Control')$cell_type_ourdatp)), 'data/lungs_all/cell_type_fraction_datp_among_AT1AT2_ctr.csv', 
          row.names = F)