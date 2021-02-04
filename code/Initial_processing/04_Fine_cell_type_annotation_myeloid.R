#!/usr/bin/env Rscript

#### Fine cell type annotation of myeloid subset
library(dplyr)
library(Seurat)

seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to myeloid, APC-like and mast cells
seu <- subset(seu, cell_type_main == 'Myeloid' | integrated_snn_res.0.8 == 26 | integrated_snn_res.0.8 == 30)
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.6)
seu <- RunUMAP(object = seu, dims = 1:20)

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.13 | seu$S.Score > 0.13, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/lungs_all/markers_lungs_all_myeloid.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(0, 1, 5, 6, 7, 8, 11), 'Macrophages', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(2), 'Monocytes', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(3, 16), 'remove', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(4, 9, 14), 'Transitioning MDM', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(13), 'Monocyte-derived macrophages', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(12, 15), 'Dendritic cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(10), 'Mast cells', seu$cell_type_fine)

# Save object and labels
saveRDS(seu, 'data/lungs_all/data_lungs_all_myeloid_initial.rds')
celltype <- seu@meta.data %>% select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/lungs_all/data_lungs_all_myeloid_celltype.csv')