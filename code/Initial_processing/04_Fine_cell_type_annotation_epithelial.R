#!/usr/bin/env Rscript

#### Fine cell type annotation of epithelial subset
library(dplyr)
library(Seurat)

seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to epithelial cells
seu <- subset(seu, cell_type_main == 'Epithelial cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.6)
seu <- RunUMAP(object = seu, dims = 1:30)

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.15, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/lungs_all/markers_lungs_all_epithelial.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(1, 2, 6), 'AT1', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(0, 3, 4, 10), 'AT2', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(11), 'ECM−high epithelial', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(12), 'Cycling epithelial', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(8), 'Airway mucous', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(7), 'Airway goblet', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(9), 'Airway club', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(5), 'Airway ciliated', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(13), 'Airway basal', seu$cell_type_fine)

# Save object and labels
saveRDS(seu, 'data/lungs_all/data_lungs_all_epithelial.rds')
celltype <- seu@meta.data %>% select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/lungs_all/data_lungs_all_epithelial_celltype.csv')