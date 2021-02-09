#!/usr/bin/env Rscript

#### Fine cell type annotation of T cell subset
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)

seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to T cells
seu <- subset(seu, cell_type_main == 'T cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.6)
seu <- RunUMAP(object = seu, dims = 1:20)

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.15 | seu$S.Score > 0.1, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/lungs_all/markers_lungs_all_tcells.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(9), 'Tregs', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(1, 2, 3, 6), 'CD4+ T cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(0, 8), 'CD8+ T cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(4, 11), 'NK cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(13, 5, 10, 12, 7), 'Cycling NK/T cells', seu$cell_type_fine)

# Save object and labels
saveRDS(seu, 'data/lungs_all/data_lungs_all_tcells.rds')
celltype <- seu@meta.data %>% select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/lungs_all/data_lungs_all_tcells_celltype.csv')