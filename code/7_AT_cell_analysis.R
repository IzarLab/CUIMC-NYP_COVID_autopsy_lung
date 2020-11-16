#!/usr/bin/env Rscript

### title: Analysis of AT1 and AT2 cells in lung samples (7 COVID-19 and 1 control samples). Code for Figure 3f and Extended Data Figure 4a-k.
### author: Jana Biermann
### date: 10/05/20

print(Sys.time())

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(gplots)
library(viridis)

pat <- 'lungs_all'

#### UMAPs and violin plots for COVID-19 and control samples (Extended Data Figure 4abcd)
#### Read in file from script #3
patient <- readRDS(paste0('data/', pat, '/data_', pat, '_with_covid_sigs.rds'))
patient@meta.data$group <- ifelse(patient@meta.data$group == 'cov', 'COVID-19', 'Control')

# Subset to AT1 and AT2 cells as identified in overallclassification and overlap with 'Epithelial cells' in 'cell_type_main'
patient <- subset(patient, overallclassification %in% c('AT1', 'AT2'))
patient <- subset(patient, cell_type_main %in% c('Epithelial cells'))

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- FindNeighbors(patient, dims = 1:20)
patient <- FindClusters(patient, resolution = 1.8)
patient <- RunUMAP(object = patient, dims = 1:20)

# Exclude outliers
patient <- subset(patient, integrated_snn_res.1.8 != 21)

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- FindNeighbors(patient, dims = 1:20)
patient <- FindClusters(patient, resolution = 0.8)
patient <- RunUMAP(object = patient, dims = 1:20)

# Identify AT1/2 cell subtypes 
# Signatures to identify DATP cells, primed and cycling AT2 cells were obtained from Choi et al. (Choi et al., 2020)
sigs <- read.csv('data/AT_sigs.csv', na.strings = c(''))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  patient <- AddModuleScore(object = patient, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Save Extended Data Figure 4abcd
pdf(file = paste0('data/', pat, '/AT/Extended_Data_Figure_4abcd.pdf'))
DimPlot(patient, reduction = 'umap', label = F, group.by = 'overallclassification', repel = T, cols = viridis(3)[1:2])
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', repel = T, cols = viridis(3)[1:2]) + ggtitle('Disease status')
FeaturePlot(patient, features = 'primed_AT2_jca1', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Primed AT2 cells')
plots <- VlnPlot(object = subset(patient, overallclassification == 'AT2'), features = c('primed_AT2_jca1', 'cycling_AT2_jca1', 'UP_in_DATPs_jca1'), pt.size = 0, combine = FALSE, group.by = 'group')
plots <- lapply(X = plots, FUN = function(p) p + scale_fill_manual(values = viridis(3)[1:2]) + theme(axis.title.x = element_blank()))
CombinePlots(plots = plots, legend = 'none', ncol = 3)
dev.off()

# P-values for Extended Data Figure 4d
AT2 <- subset(patient, overallclassification == 'AT2')
paste('P-value for primed AT2 cell signature:', wilcox.test(AT2@meta.data$primed_AT2_jca1 ~ AT2@meta.data$group)$p.val)
paste('P-value for cycling AT2 cell signature:', wilcox.test(AT2@meta.data$cycling_AT2_jca1 ~ AT2@meta.data$group)$p.val)
paste('P-value for DATP AT2 cell signature:', wilcox.test(AT2@meta.data$UP_in_DATPs_jca1 ~ AT2@meta.data$group)$p.val)


#### Diffusion component analysis for COVID-19 and control samples (Extended Data Figure 4hijk) ####
es <- as.ExpressionSet(as.data.frame(t(patient@assays$integrated@data)))
es@phenoData@data <- patient@meta.data

# Make a diffusion map
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Add DCs to Seurat object
patient <- AddMetaData(patient, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

# Save Extended Data Figure 4hijk
pdf(paste0('data/', pat, '/AT/Extended_Data_Figure_4hijk.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(viridis(3))
plot(dm, col = as.factor(es@phenoData@data$overallclassification), main = 'Cell type COVID-19 and control cells', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$overallclassification)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$overallclassification))), bty = 'n')

plot(dm, col = as.factor(es@phenoData@data$group), main = 'Disease status COVID-19 and control cells', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$group)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$group))), bty = 'n')

palette('Okabe-Ito')
plot(dm, col = as.factor(es@phenoData@data$patient), main = 'Patient COVID-19 and control cells')
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(es@phenoData@data$patient)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$patient))), bty = 'n')

palette(viridis(3))
par(mar = c(5.1, 4.1, 4.1, 6.5), xpd = TRUE)
plot(dm, col = es@phenoData@data$UP_in_DATPs_jca1, main = 'DATP COVID-19 and control cells', pch = 20)
colorlegend(es@phenoData@data$UP_in_DATPs_jca1, viridis(length(es@phenoData@data$UP_in_DATPs_jca1)), posx = c(0.88, 0.9), posy = c(0, 0.65))
dev.off()


#### COVID-19 AT1/2 cell analysis #### 
#Subset the data set to include only COVID-19 samples
covid <- subset(patient, group %in% c('COVID-19'))

# Rerun Seurat workflow
covid <- ScaleData(covid)
covid <- RunPCA(covid)
covid <- FindNeighbors(covid, dims = 1:20)
covid <- FindClusters(covid, resolution = 0.8)
covid <- RunUMAP(covid, dims = 1:20)

# Rerun signatures to identify AT1/2 cell subtypes
# Signatures to identify DATP cells, primed and cycling AT2 cells were obtained from Choi et al. (Choi et al., 2020)
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  covid <- AddModuleScore(object = covid, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Diffusion component analysis for COVID-19 samples (Figure 3f and Extended Data Figure 4fg)
es_cov <- as.ExpressionSet(as.data.frame(t(covid@assays$integrated@data)))
es_cov@phenoData@data <- covid@meta.data

# Make a diffusion map
dm_cov <- DiffusionMap(es_cov, verbose = T, n_pcs = 20)

# Add DCs to Seurat object
covid <- AddMetaData(covid, dm_cov@eigenvectors, col.name = colnames(dm_cov@eigenvectors))

# Save Extended Data Figure 4fg
pdf(paste0('data/', pat, '/AT/Extended_Data_Figure_4fg.pdf'))
palette(viridis(3))
par(mar = c(5.1, 4.1, 4.1, 6.5), xpd = TRUE)
plot(dm_cov, col = es_cov@phenoData@data$primed_AT2_jca1, main = 'Primed AT2 COVID-19 cells', pch = 20)
colorlegend(es_cov@phenoData@data$primed_AT2_jca1, viridis(length(es_cov@phenoData@data$primed_AT2_jca1)), posx = c(0.88, 0.9), posy = c(0, 0.65))

palette(viridis(3))
par(mar = c(5.1, 4.1, 4.1, 6.5), xpd = TRUE)
plot(dm_cov, col = es_cov@phenoData@data$UP_in_DATPs_jca1, main = 'DATP COVID-19 cells', pch = 20)
colorlegend(es_cov@phenoData@data$UP_in_DATPs_jca1, viridis(length(es_cov@phenoData@data$UP_in_DATPs_jca1)), posx = c(0.88, 0.9), posy = c(0, 0.65))
dev.off()

# Save Figure 3f
pdf(paste0('data/', pat, '/AT/Figure_3f.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(viridis(3))
plot(dm_cov, col = as.factor(es_cov@phenoData@data$overallclassification), main = 'Cell type COVID-19 cells', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es_cov@phenoData@data$overallclassification)), pch = 16, 
       col = as.factor(levels(as.factor(es_cov@phenoData@data$overallclassification))), bty = 'n')
dev.off()


#### Control AT1/2 cell analysis ####
# Subset the data set to include only control samples
control <- subset(patient, group %in% c('Control'))

# Rerun Seurat workflow
control <- ScaleData(control)
control <- RunPCA(control)
control <- FindNeighbors(control, dims = 1:20)
control <- FindClusters(control, resolution = 0.8)
control <- RunUMAP(control, dims = 1:20)

# Rerun signatures to identify AT1/2 cell subtypes 
# Signatures to identify DATP cells, primed and cycling AT2 cells were obtained from Choi et al. (Choi et al., 2020)
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  control <- AddModuleScore(object = control, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Diffusion component analysis for control samples (Extended Data Figure 4e)
es_ctr <- as.ExpressionSet(as.data.frame(t(control@assays$integrated@data)))
es_ctr@phenoData@data <- control@meta.data

# Make a diffusion map
dm_ctr <- DiffusionMap(es_ctr, verbose = T, n_pcs = 20)

# Add DCs to Seurat object
control <- AddMetaData(control, dm_ctr@eigenvectors, col.name = colnames(dm_ctr@eigenvectors))

# Save Extended Data Figure 4e
pdf(paste0('data/', pat, '/AT/Extended_Data_Figure_4e.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(viridis(3))
plot(dm_ctr, col = as.factor(es_ctr@phenoData@data$overallclassification), main = 'Cell type control cells', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es_ctr@phenoData@data$overallclassification)), pch = 16, 
       col = as.factor(levels(as.factor(es_ctr@phenoData@data$overallclassification))), bty = 'n')
dev.off()


##### Save AT cell image #####
save.image(paste0('data/', pat, '/AT/image_AT_cells.RData'))


print(Sys.time())