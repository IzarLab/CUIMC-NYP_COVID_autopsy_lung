#!/usr/bin/env Rscript

### title: Analysis of myeloid cells in lung samples (7 COVID-19 and 1 control samples). Code for Figure 2abd, Extended Data Figure 2, Extended Data Table 2 and 3.
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

#### UMAP analysis (Figure 2a) ####
# Read in file from script #3
patient <- readRDS(paste0('data/', pat, '/data_', pat, '_with_covid_sigs.rds'))
patient@meta.data$group <- ifelse(patient@meta.data$group == 'cov', 'COVID-19', 'Control')

# Rename disease status
patient@meta.data$group <- ifelse(patient@meta.data$group == 'cov', 'COVID-19', 'Control')

# Subset to monocytes and macrophages
patient <- subset(patient, overallclassification %in% c('Monocytes', 'Macrophages'))

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- FindNeighbors(patient, dims = 1:20)
patient <- FindClusters(patient, resolution = 0.8)
patient <- RunUMAP(object = patient, dims = 1:20)

# Exclude outliers
patient <- subset(patient, integrated_snn_res.0.8 != 15)

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- FindNeighbors(patient, dims = 1:20)
patient <- FindClusters(patient, resolution = 0.8)
patient <- RunUMAP(object = patient, dims = 1:20)


# Identify macrophage lineage
# Gene signature to identify alveolar macrophages (AM) based on differential gene expression was obtained from Travaglini et al. (Travaglini et al., 2020)
sigs <- read.csv('data/myeloid_signatures.csv', na.strings = c(''))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  patient <- AddModuleScore(object = patient, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}
patient@meta.data$alveolar_mac <- ifelse(patient@meta.data$overallclassification == 'Monocytes', patient@meta.data$overallclassification, 
                                         ifelse(patient@meta.data$alveolar_mac_DEG1 > 0.2, 'Alveolar macrophages', 'Monocyte-derived macrophages'))

# Save UMAP (Figure 2a)
pdf(file = paste0('data/', pat, '/myeloid/Figure_2a.pdf'))
plot <- DimPlot(patient, reduction = 'umap', label = F, group.by = 'alveolar_mac', repel = T, cols = viridis(3))
AugmentPlot(plot, dpi = 300)
dev.off()


#### Diffusion component analysis (Figure 2b) ####
es <- as.ExpressionSet(as.data.frame(t(patient@assays$integrated@data)))
es@phenoData@data <- patient@meta.data

# Make a diffusion map
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Add DCs to Seurat object
patient <- AddMetaData(patient, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

# Save diffusion map (Figure 2b)
pdf(paste0('data/', pat, '/myeloid/Figure_2b.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(viridis(3))
plot(dm, col = as.factor(es@phenoData@data$alveolar_mac), main = 'Myeloid cell lineages', pch = 20)
legend('bottomright', inset = c(-0.05, -0.15), legend = levels(as.factor(es@phenoData@data$alveolar_mac)), pch = 16, 
       col = as.factor(levels(as.factor(es@phenoData@data$alveolar_mac))), bty = 'n')
dev.off()


#### Differential gene expression in alveolar macrophages (Extended Data Table 3 and Figure 2d) ####
# Subset to only AMs
am <- subset(patient, alveolar_mac == 'Alveolar macrophages')
DefaultAssay(am) <- 'RNA'
am <- NormalizeData(am)
am <- ScaleData(am)

# Identify DEGs between COVID-19 and control donors (Extended Data Table 3)
Idents(object = am) <- am@meta.data$group
markers <- FindAllMarkers(am, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
markers %>% group_by(cluster) %>% write.csv(paste0('data/', pat, '/myeloid/Extended_Data_Table_3.csv'), row.names = F)
markers <- markers %>% filter(p_val_adj < 0.05)
top50 <- markers[markers$cluster == 'COVID-19', ]$gene[1:25]
top50 <- c(top50, markers[markers$cluster == 'Control', ]$gene[1:25])

# Generate heatmap of top 50 DEGs (top 25 from each group) (Figure 2d)
pdf(paste0('data/', pat, '/myeloid/Figure_2d.pdf'))
DoHeatmap(am, features = top50, group.by = 'group', assay = 'RNA', size = 3,angle = 0) + ggtitle('Top 50 DEGs in alveolar macrophages')
dev.off()


#### Macrophage analysis (Extended Data Figure 2) #### 
# Subset to only macrophages
macs <- subset(patient, alveolar_mac %in% c('Monocyte-derived macrophages', 'Alveolar macrophages'))

# Rerun Seurat workflow
macs <- ScaleData(macs)
macs <- RunPCA(macs)
macs <- FindNeighbors(macs, dims = 1:20)
macs <- FindClusters(macs, resolution = 0.5)
macs <- RunUMAP(macs, dims = 1:20)


# Identify macrophage subtypes based on signatures
# Signatures to identify subtypes of macrophages were obtained from Mould et al. (Mould et al., 2019)
sigs2 <- read.csv('data/myeloid2.csv', na.strings = c(''))

# Basic function to convert mouse to human gene names 
# Source: https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/
convertMouseGeneList <- function(x) {
  require('biomaRt')
  human = useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
  mouse = useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
  genesV2 = getLDS(attributes = c('mgi_symbol'), filters = 'mgi_symbol', values = x, mart = mouse, attributesL = c('hgnc_symbol'), martL = human, uniqueRows = T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
sigs <- c()
sigs$c1_resi <- convertMouseGeneList(sigs2$c1_resi)
sigs$c2_resi_proliferation <- convertMouseGeneList(sigs2$c2_resi_proliferation)
sigs$c3_BM_proteasome <- convertMouseGeneList(sigs2$c3_BM_proteasome)
sigs$c4_BM_inflammaroty_signaling <- convertMouseGeneList(sigs2$c4_BM_inflammaroty_signaling)
sigs$c5_BM_growth <- convertMouseGeneList(sigs2$c5_BM_growth)
sigs <- as.data.frame(t(plyr::ldply(sigs, rbind)))
colnames(sigs) <- sigs[1, ]
sigs <- sigs[-1, ]

for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  macs <- AddModuleScore(object = macs, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}


# Extended Data Figure 2 a-f
pdf(file = paste0('data/', pat, '/myeloid/Extended_Data_Figure_2_abcdef.pdf'))
plot <- DimPlot(macs, reduction = 'umap', label = F, group.by = 'alveolar_mac', repel = T, cols = viridis(3))
AugmentPlot(plot, dpi = 300)
plot <- FeaturePlot(macs, features = 'c1_resi1', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('c1 Residential macrophages')
AugmentPlot(plot, dpi = 300)
plot <- FeaturePlot(macs, features = 'c2_resi_proliferation1', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('c2 Proliferation residential macrophages')
AugmentPlot(plot, dpi = 300)
plot <- FeaturePlot(macs, features = 'c3_BM_proteasome1', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('c3 Monocyte-derived macrophages proteasome')
AugmentPlot(plot, dpi = 300)
plot <- FeaturePlot(macs, features = 'c4_BM_inflammaroty_signaling1', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('c4 Monocyte-derived macrophages proinflammatory')
AugmentPlot(plot, dpi = 300)
plot <- FeaturePlot(macs, features = 'c5_BM_growth1', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('c5 Monocyte-derived macrophages growth factors')
AugmentPlot(plot, dpi = 300)
dev.off()


# Differential gene expression in macrophages (Extended Data Table 2 and Extended Data Figure 2g)
DefaultAssay(macs) <- 'RNA'
macs <- NormalizeData(macs)
macs <- ScaleData(macs)

# Identify markers based on macrophage lineage
Idents(macs) <- macs@meta.data$alveolar_mac
markers <- FindAllMarkers(macs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
markers %>% group_by(cluster) %>% write.csv(paste0('data/', pat, '/myeloid/Extended_Data_Table_2.csv'), row.names = F)
markers <- markers %>% filter(p_val_adj < 0.05)
top50 <- markers[markers$cluster == 'Alveolar macrophages', ]$gene[1:25]
top50 <- c(top50, markers[markers$cluster == 'Monocyte-derived macrophages', ]$gene[1:25])

# Save Extended Data Figure 2g
pdf(paste0('data/', pat, '/myeloid/Extended_Data_Figure_2g.pdf'))
DoHeatmap(macs, features = top50, group.by = 'alveolar_mac', assay = 'RNA', size = 3,angle = 0) + ggtitle('Top 50 DEGs between AM and monocyte-derived macrophages')
dev.off()


##### Save myeloid cell image #####
save.image(paste0('data/', pat, '/myeloid/image_myeloid.RData'))


print(Sys.time())
