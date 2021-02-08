#!/usr/bin/env Rscript

#### Fine cell type annotation of endothelial subset
library(dplyr)
library(Seurat)

seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to endothelial cells
seu <- subset(seu, cell_type_main == 'Endothelial cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:20)

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.07 | seu$S.Score > 0.07, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/lungs_all/markers_lungs_all_endothelial.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0), 'Endothelial cells (general)', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(1), 'Capillary endothelial cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(2), 'Arterial endothelial cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(3), 'Systemic venous endothelial cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(4), 'Inflamed endothelial cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(5), 'Pulmonary venous endothelial cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(6), 'Endothelial cells (other)', seu$cell_type_fine)

# Save object and labels
saveRDS(seu, 'data/lungs_all/data_lungs_all_endothelial.rds')
celltype <- seu@meta.data %>% select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/lungs_all/data_lungs_all_endothelial_celltype.csv')